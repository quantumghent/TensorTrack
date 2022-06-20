classdef MatrixBlock < AbstractBlock
    %MATRIXBLOCK Summary of this class goes here
    %   Detailed explanation goes here
    
    %#ok<*INUSD>
    properties
        charge
        var
        rowsizes
        colsizes
        tdims
        rank
    end
    
    methods
        function b = MatrixBlock(codomain, domain)
            if nargin == 0, return; end
            
            rank = [length(codomain) length(domain)];
            
            trees = fusiontrees(codomain, domain);
            [c, ~, ic] = unique(trees.coupled);
            uncoupled = trees.uncoupled;
            
            spaces = [codomain flip(domain)];
            ch = charges(spaces).';
            d = degeneracies(spaces).';
            
            [lia, locb] = ismember(uncoupled, ch, 'rows');
            tdims = d(locb(lia), :);
%             tdims = computedegeneracies(codomain, domain, uncoupled);
            mdims = [prod(tdims(:, 1:rank(1)), 2) prod(tdims(:, (1:rank(2)) + rank(1)), 2)];
            splits = split(trees);
            fuses = fuse(trees);
            
            for i = length(c):-1:1
                ids = ic == i;
                
                [~, ia] = unique(splits(ids));
                rowdims = mdims(ids, 1);
                rowsizes = [0 cumsum(rowdims(ia)).'];
                
                [~, ia] = unique(fuses(ids));
                coldims = mdims(ids, 2);
                colsizes = [0 cumsum(coldims(ia)).'];
                
                b(i).charge = c(i);
                b(i).var = uninit([rowsizes(end) colsizes(end)]);
                b(i).rowsizes = rowsizes;
                b(i).colsizes = colsizes;
                b(i).tdims = tdims(ids, :);
                b(i).rank = rank;
            end
        end
        
        function Y = axpby(a, X, b, Y, p, map)
            % Compute ```Y = permute(X, p) .* a + Y .* b```.
            % This method is the computationally critical method of this class, thus has
            % special cases for scalar multiplication (a == 0), addition (nargin == 4), and
            % various optimisations when a == 1, b == 0 | b == 1.
            % 
            % Arguments
            % ---------
            % a : double
            %   scalar to multiply with X.
            %
            % X : MatrixBlock
            %   list of source blocks.
            %
            % b : double
            %   scalar to multiply with Y.
            %
            % Y : MatrixBlock
            %   list of destination blocks.
            %
            % p : int
            %   permutation vector for X.
            %
            % map : (sparse) double
            %   coefficient matrix for permuting X.
            %
            % Returns
            % -------
            % Y : MatrixBlock
            %   Result of computing Y = permute(X, p) .* a + Y .* b.
            
            %% Special case 1: scalar multiplication
            if a == 0
                Y = Y .* b;
                return
            end
            
            
            %% Special case 2: addition without permutation
            if nargin == 4 || (isempty(p) && isempty(map)) || ...
                    (all(p == 1:length(p)) && all(X(1).rank == Y(1).rank))
                % reduce to scalar multiplication
                if b == 0   % a ~= 0 -> case 1
                    Y = X .* a;
                    return
                end
                
                % reduce to simple addition or substraction
                if a == 1 && b == 1,    Y = X + Y;  return; end
                if a == 1 && b == -1,   Y = X - Y;  return; end
                if a == -1 && b == 1,   Y = Y - X;  return; end
                
                % general addition + scalar multiplication cases
                if a == 1   % b ~= 1
                    for i = 1:length(Y)
                        Y(i).var = Y(i).var .* b + X(i).var;
                    end
                    return
                end
                
                if b == 1   % a ~= 1
                    for i = 1:length(Y)
                        Y(i).var = Y(i).var + X(i).var .* a;
                    end
                    return
                end
                
                % a ~= [0 1], b ~= [0 1]
                for i = 1:length(Y)
                    Y(i).var = Y(i).var .* b + X(i).var .* a;
                end
                return
            end
            
            
            %% General case: addition with permutation
            % tensor indexing to matrix indexing
            rrx = rankrange(X(1).rank);
            rry = rankrange(Y(1).rank);
            p_eff(rry) = rrx(p);
            
            % extract small tensor blocks
            doA = a ~= 1;
            vars_in = cell(size(map, 1), 1);
            offset = 0;
            for i = 1:length(X)
                rowsz = X(i).rowsizes;
                colsz = X(i).colsizes;
                
                var_in = X(i).var;
                if doA, var_in = a .* var_in; end
                
                oldtdims = X(i).tdims(:, rrx);
                newmdims = [prod(X(i).tdims(:, p_eff(1:Y(1).rank(1))), 2) ...
                    prod(X(i).tdims(:, p_eff((1:Y(1).rank(2)) + Y(1).rank(1))), 2)];
                
                for k = 1:length(colsz) - 1
                    for j = 1:length(rowsz) - 1
                        ind = j + (k-1) * (length(rowsz) - 1);
                        offset = offset + 1;
                        
                        vars_in{offset} = reshape(permute(reshape(...
                            var_in(rowsz(j)+1:rowsz(j+1), colsz(k)+1:colsz(k+1)), ...
                            oldtdims(ind, :)), ...
                            p_eff), ...
                            newmdims(ind, :));
                    end
                end
            end
            
            % apply map
            vars_out = cell(size(map, 2), 1);
            [rows, cols, vals] = find(map);
            for i = 1:length(vals)
                if isempty(vars_out{cols(i)})
                    vars_out{cols(i)} = vals(i) * vars_in{rows(i)};
                else
                    vars_out{cols(i)} = vars_out{cols(i)} + vals(i) * vars_in{rows(i)};
                end
            end
            
            % inject small tensor blocks
            if b == 0
                offset = 0;
                for i = 1:length(Y)
                    rows = length(Y(i).rowsizes) - 1;
                    cols = length(Y(i).colsizes) - 1;
                    if rows < cols
                        m = cell(rows, 1);
                        for n = 1:rows
                            m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                        end
                        Y(i).var = cat(1, m{:});
                    else
                        m = cell(cols, 1);
                        for n = 1:cols
                            m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                        end
                        Y(i).var = cat(2, m{:});
                    end
                    offset = offset + rows * cols;
                end
                return
            end
            
            if b == 1
                offset = 0;
                for i = 1:length(Y)
                    rows = length(Y(i).rowsizes) - 1;
                    cols = length(Y(i).colsizes) - 1;
                    if rows < cols
                        m = cell(rows, 1);
                        for n = 1:rows
                            m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                        end
                        Y(i).var = Y(i).var + cat(1, m{:});
                    else
                        m = cell(cols, 1);
                        for n = 1:cols
                            m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                        end
                        Y(i).var = Y(i).var + cat(2, m{:});
                    end
                    offset = offset + rows * cols;
                end
                return
            end
            
            offset = 0;
            for i = 1:length(Y)
                rows = length(Y(i).rowsizes) - 1;
                cols = length(Y(i).colsizes) - 1;
                if rows < cols
                    m = cell(rows, 1);
                    for n = 1:rows
                        m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                    end
                    Y(i).var = b .* Y(i).var + cat(1, m{:});
                else
                    m = cell(cols, 1);
                    for n = 1:cols
                        m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                    end
                    Y(i).var = b .* Y(i).var + cat(2, m{:});
                end
                offset = offset + rows * cols;
            end
        end
        
        function C = mul(C, A, B, a, b)
            % Compute C = (A .* a) * (B .* b).
            %
            % Arguments
            % ---------
            % C : MatrixBlock
            %   location to store the result
            %
            % A : MatrixBlock
            %   first matrix factor
            %
            % B : MatrixBlock
            %   second matrix factor
            %
            % a : double = 1
            %   first scalar factor
            %
            % b : double = 1
            %   second scalar factor
            %
            % Returns
            % -------
            % C : MatrixBlock
            %   Result of computing C = (A .* a) * (B .* b)
            
            [indA, locA] = ismember([C.charge], [A.charge]);
            [indB, locB] = ismember([C.charge], [B.charge]);
            
            if nargin == 3 || (a == 1 && b == 1)
                for i = find(indA & indB)
                    C(i).var = A(locA(i)).var * B(locB(i)).var;
                end
                return
            end
            
            if a == 1 % b ~= 1
                for i = find(indA & indB)
                    C(i).var = (A(locA(i)).var .* a) * B(locB(i)).var;
                end
                return
            end
            
            if b == 1 % a ~= 1
                for i = find(indA & indB)
                    C(i).var = A(locA(i)).var * (B(locB(i)).var .* b);
                end
                return
            end
            
            % a ~= 1 && b ~= 1
            for i = find(indA & indB)
                C(i).var = (A(locA(i)).var .* a) * (B(locB(i)).var .* b);
            end
        end
        
        function tblocks = tensorblocks(b)
            % Extract the non-empty small tensor blocks.
            %
            % Arguments
            % ---------
            % b : MatrixBlock
            %   Big blocks to convert.
            %
            % Returns
            % -------
            % tblocks : cell
            %   list of non-zero small tensor blocks, in column-major order.
            
            nblocks = 0;
            for i = 1:length(b)
                nblocks = nblocks + size(b(i).tdims, 1);
            end
            tblocks = cell(nblocks, 1);
            
            offset = 0;
            p = rankrange(b(1).rank);
            for i = 1:length(b)
                rowsz = b(i).rowsizes;
                colsz = b(i).colsizes;
                for k = 1:length(colsz) - 1
                    for j = 1:length(rowsz) - 1
                        offset = offset + 1;
                        tblocks{offset} = permute(reshape(...
                            b(i).var(rowsz(j)+1:rowsz(j+1), colsz(k)+1:colsz(k+1)), ...
                            b(i).tdims(j + (k-1) * (length(rowsz)-1), p)), ...
                            p);
                    end
                end
            end
        end
        
        function [mblocks, mcharges] = matrixblocks(b)
            % Extract a list of coupled matrix blocks.
            %
            % Arguments
            % ---------
            % b : MatrixBlock
            %   list of input data.
            %
            % Returns
            % -------
            % mblocks : cell
            %   list of non-zero coupled matrix blocks, sorted according to its charge.
            %
            % mcharges : AbstractCharge
            %   list of coupled charges.
            
            mblocks = {b.var};
            if nargout > 1
                mcharges = [b.charge];
            end
        end
        
        function b = fill_matrix_data(b, vars, charges)
            if nargin < 3 || isempty(charges)
                assert(length(vars) == length(b));
                for i = 1:length(b)
                    b(i).var = vars{i};
                end
                return
            end
            
            [lia, locb] = ismember(charges, [b.charge]);
            assert(all(lia));
            for i = 1:length(vars)
                b(locb(i)).var = vars{i};
            end
        end
        
        function b = fill_matrix_fun(b, fun, charges)
            if nargin < 3 || isempty(charges)
                for i = 1:length(b)
                    b(i).var = fun(size(b(i).var), b(i).charge);
                end
            else
                [lia, locb] = ismember(charges, [b.charge]);
                for i = locb(lia)
                    b(i).var = fun(size(b(i).var), b(i).charge);
                end
            end
        end
        
        function Y = minus(X, Y)
            % Subtraction of X and Y.
            %
            % Usage
            % -----
            % Y = minus(X, Y)
            % Y = X - Y
            %
            % Arguments
            % ---------
            % X : MatrixBlock
            %   first list of input matrices.
            %
            % Y : MatrixBlock
            %   second list of input matrices.
            %
            % Returns
            % -------
            % Y : MatrixBlock
            %   list of output matrices.
            
            for i = 1:length(Y)
                Y(i).var = X(i).var - Y(i).var;
            end
        end
        
        function Y = plus(X, Y)
            % Addition of X and Y.
            %
            % Usage
            % -----
            % Y = plus(X, Y)
            % Y = X + Y
            %
            % Arguments
            % ---------
            % X : MatrixBlock
            %   first list of input matrices.
            %
            % Y : MatrixBlock
            %   second list of input matrices.
            %
            % Returns
            % -------
            % Y : MatrixBlock
            %   list of output matrices.
            
            for i = 1:length(Y)
                Y(i).var = Y(i).var + X(i).var;
            end
        end
        
        function Y = times(Y, a)
            % Scalar multiplication of Y and a.
            %
            % Usage
            % -----
            % Y = times(Y, a)
            % Y = Y .* a
            %
            % Arguments
            % ---------
            % Y : MatrixBlock
            %   list of input matrices.
            %
            % a : double
            %   scalar factor.
            %
            % Returns
            % -------
            % Y : MatrixBlock
            %   list of output matrices.
            
            if a == 1, return; end
            if a == -1, Y = -Y; return; end
            
            for i = 1:length(Y)
                Y(i).var = Y(i).var .* a;
            end
        end
        
        function Y = rdivide(Y, a)
            % Scalar division of Y and a.
            %
            % Usage
            % -----
            % Y = rdivide(Y, a)
            % Y = Y ./ a
            %
            % Arguments
            % ---------
            % Y : MatrixBlock
            %   list of input matrices.
            %
            % a : double
            %   scalar factor.
            %
            % Returns
            % -------
            % Y : MatrixBlock
            %   list of output matrices.
            
            if a == 1, return; end
            if a == -1, Y = -Y; return; end
            
            for i = 1:length(Y)
                Y(i).var = Y(i).var ./ a;
            end
        end
        
        function A = uplus(A)
            % Unary plus. Equivalent to making a copy.
            %
            % Usage
            % -----
            % A = uplus(A)
            % A = +A
            %
            % Arguments
            % ---------
            % A : MatrixBlock
            %   list of input matrices.
            %
            % Returns
            % -------
            % A : MatrixBlock
            %   list of output matrices.
            
        end
        
        function Y = uminus(Y)
            % Unary minus. Computes the additive inverse.
            %
            % Usage
            % -----
            % A = uminus(A)
            % A = -A
            %
            % Arguments
            % ---------
            % A : MatrixBlock
            %   list of input matrices.
            %
            % Returns
            % -------
            % A : MatrixBlock
            %   list of output matrices.
            
            for i = 1:length(Y)
                Y(i).var = -Y(i).var;
            end
        end
        
        function t = underlyingType(X)
            t = underlyingType(X(1).var);
        end
        
        function X = ctranspose(X)
            % Adjoint of a tensor.
            %
            % Usage
            % -----
            % X = ctranspose(X)
            % X = X'
            %
            % Arguments
            % ---------
            % X : :class:`MatrixBlock`
            %   list of input matrices.
            %
            % Returns
            % -------
            % X : :class:`MatrixBlock`
            %   list of adjoint output matrices.
            
            for i = 1:length(X)
                X(i).var = X(i).var';
                [X(i).rowsizes, X(i).colsizes] = swapvars(X(i).rowsizes, X(i).colsizes);
                X(i).tdims = fliplr(X(i).tdims);
                X(i).rank = fliplr(X(i).rank);
            end
        end
        
        function v = vectorize(X, type)
            switch type
                case 'complex'
                    blocks = cell(size(X));
                    for i = 1:length(X)
                        blocks{i} = reshape(X(i).var .* sqrt(qdim(X(i).charge)), [], 1);
                    end
                    v = vertcat(blocks{:});
                    
                case 'real'
                    blocks = cell(size(X));
                    for i = 1:length(X)
                        X(i).var = X(i).var .* sqrt(qdim(X(i).charge));
                        blocks{i} = [reshape(real(X(i).var), [], 1)
                            reshape(imag(X(i).var), [], 1)];
                    end
                    v = vertcat(blocks{:});
            end
        end
        
        function X = devectorize(v, X, type)
            switch type
                case 'complex'
                    ctr = 0;
                    for i = 1:length(X)
                        n = numel(X(i).var);
                        X(i).var = reshape(v(ctr + (1:n)), size(X(i).var)) ./ ...
                            sqrt(qdim(X(i).charge));
                        ctr = ctr + n;
                    end
                    
                case 'real'
                    ctr = 0;
                    for i = 1:length(X)
                        n = numel(X(i).var);
                        sz = size(X(i).var);
                        X(i).var = complex(reshape(v(ctr + (1:n)), sz), ...
                            reshape(v(ctr + (n + 1:2 * n)), sz)) ./ sqrt(qdim(X(i).charge));
                        ctr = ctr + 2 * n;
                    end
            end
        end
    end
end
