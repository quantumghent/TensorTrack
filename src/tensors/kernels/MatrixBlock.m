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
            if isempty(trees)
                warning('tensors:empty', ...
                    'No fusion channels available for the given spaces.');
%                 b = MatrixBlock.empty(0, 1);
                return
            end
            assert(~isempty(trees), 'tensors:empty', ...
                'no fusion channels available for the given spaces.');
            
            [c, ~, ic] = unique(trees.coupled);
            uncoupled = trees.uncoupled;
            
            spaces = [codomain flip(domain)];
            ch = charges(spaces).';
            d = degeneracies(spaces).';
            
            [lia, locb] = ismember(uncoupled, ch, 'rows');
            tdims = d(locb(lia), :);
            mdims = [prod(tdims(:, 1:rank(1)), 2) prod(tdims(:, (1:rank(2)) + rank(1)), 2)];
            splits = split(trees);
            fuses = fuse(trees);
            
            b.rank = rank;
            b.charge = reshape(c, 1, []);
            b.var       = cell(size(b.charge));
            b.rowsizes  = cell(size(b.charge));
            b.colsizes  = cell(size(b.charge));
            b.tdims     = cell(size(b.charge));
            
            for i = length(c):-1:1
                ids = ic == i;
                
                [~, ia] = unique(splits(ids));
                rowdims = mdims(ids, 1);
                rowsizes = [0 cumsum(rowdims(ia)).'];
                
                [~, ia] = unique(fuses(ids));
                coldims = mdims(ids, 2);
                colsizes = [0 cumsum(coldims(ia)).'];
                
                b.var{i} = uninit([rowsizes(end) colsizes(end)]);
                b.rowsizes{i} = rowsizes;
                b.colsizes{i} = colsizes;
                b.tdims{i} = tdims(ids, :);
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
            rx = X.rank;
            ry = Y.rank;
            
            if nargin == 4 || (isempty(p) && isempty(map)) || ...
                    (all(p == 1:length(p)) && isequal(rx, ry) && ...
                    isequal(map, speye(size(map))))
                % reduce to scalar multiplication
                if b == 0   % a ~= 0 -> case 1
                    Y = X .* a;
                    return
                end
                
                % reduce to simple addition or substraction
                if a == 1 && b == 1,    Y = X + Y;  return; end
                if a == 1 && b == -1,   Y = X - Y;  return; end
                if a == -1 && b == 1,   Y = Y - X;  return; end
                
                Yvar = Y.var;
                Xvar = X.var;
                
                % general addition + scalar multiplication cases
                if a == 1       % b ~= 1
                    for i = 1:length(Yvar)
                        Yvar{i} = Yvar{i} .* b + Xvar{i};
                    end
                elseif b == 1   % a ~= 1
                    for i = 1:length(Yvar)
                        Yvar{i} = Yvar{i} + Xvar{i} .* a;
                    end
                else            % a ~= [0 1], b ~= [0 1]
                    for i = 1:length(Yvar)
                        Yvar{i} = Yvar{i} .* b + Xvar{i} .* a;
                    end
                end
                
                Y.var = Yvar;
                return
            end
            
            
            %% General case: addition with permutation
            % tensor indexing to matrix indexing
            rrx = rankrange(rx);
            rry = rankrange(ry);
            p_eff(rry) = rrx(p);
            p_eff_1 = p_eff(1:ry(1));
            p_eff_2 = p_eff((1:ry(2)) + ry(1));
            
            % extract small tensor blocks
            doA = a ~= 1;
            vars_in = cell(size(map, 1), 1);
            offset = 0;
            
            Xrowsizes = X.rowsizes;
            Xcolsizes = X.colsizes;
            Xtdims = X.tdims;
            Xvar = X.var;
            for i = 1:length(Xvar)
                rowsz = Xrowsizes{i};
                colsz = Xcolsizes{i};
                
                var_in = Xvar{i};
                if doA, var_in = a .* var_in; end
                
                oldtdims = Xtdims{i}(:, rrx);
                newmdims = [prod(oldtdims(:, p_eff_1), 2) prod(oldtdims(:, p_eff_2), 2)];
                
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
            Yrowsizes = Y.rowsizes;
            Ycolsizes = Y.colsizes;
            Yvar = Y.var;
            if b == 0
                offset = 0;
                for i = 1:length(Yvar)
                    rows = length(Yrowsizes{i}) - 1;
                    cols = length(Ycolsizes{i}) - 1;
                    if rows < cols
                        m = cell(rows, 1);
                        for n = 1:rows
                            m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                        end
                        Yvar{i} = cat(1, m{:});
                    else
                        m = cell(cols, 1);
                        for n = 1:cols
                            m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                        end
                        Yvar{i} = cat(2, m{:});
                    end
                    offset = offset + rows * cols;
                end
                Y.var = Yvar;
                return
            end
            
            if b == 1
                offset = 0;
                for i = 1:length(Yvar)
                    rows = length(Yrowsizes{i}) - 1;
                    cols = length(Ycolsizes{i}) - 1;
                    if rows < cols
                        m = cell(rows, 1);
                        for n = 1:rows
                            m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                        end
                        Yvar{i} = Yvar{i} + cat(1, m{:});
                    else
                        m = cell(cols, 1);
                        for n = 1:cols
                            m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                        end
                        Yvar{i} = Yvar{i} + cat(2, m{:});
                    end
                    offset = offset + rows * cols;
                end
                return
            end
            
            offset = 0;
            for i = 1:length(Yvar)
                rows = length(Yrowsizes{i}) - 1;
                    cols = length(Ycolsizes{i}) - 1;
                if rows < cols
                    m = cell(rows, 1);
                    for n = 1:rows
                        m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                    end
                    Yvar{i} = b .* Yvar{i} + cat(1, m{:});
                else
                    m = cell(cols, 1);
                    for n = 1:cols
                        m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                    end
                    Yvar{i} = b .* Yvar{i} + cat(2, m{:});
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
            
            Acharge = A.charge;
            Bcharge = B.charge;
            Ccharge = C.charge;
            
            if isequal(Acharge, Ccharge)
                indA = true(size(Acharge));
                locA = 1:length(Acharge);
            else
                [indA, locA] = ismember_sorted(Ccharge, Acharge);
            end
            if isequal(Bcharge, Ccharge)
                indB = true(size(Bcharge));
                locB = 1:length(Bcharge);
            else
                [indB, locB] = ismember_sorted(Ccharge, Bcharge);
            end
            
            Avar = A.var;
            Bvar = B.var;
            Cvar = C.var;
            
            if nargin == 3 || (a == 1 && b == 1)
                for i = find(indA & indB)
                    Cvar{i} = Avar{locA(i)} * Bvar{locB(i)};
                end
                C.var = Cvar;
                return
            end
            
            if a == 1 % b ~= 1
                for i = find(indA & indB)
                    Cvar{i} = (Avar{locA(i)} .* a) * Bvar{locB(i)};
                end
                C.var = Cvar;
                return
            end
            
            if b == 1 % a ~= 1
                for i = find(indA & indB)
                    Cvar{i} = Avar{locA(i)} * (Bvar{locB(i)} .* b);
                end
                C.var = Cvar;
                return
            end
            
            % a ~= 1 && b ~= 1
            for i = find(indA & indB)
                Cvar{i} = (Avar{locA(i)} .* a) * (Bvar{locB(i)} .* b);
            end
            C.var = Cvar;
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
            for i = 1:length(b.var)
                nblocks = nblocks + size(b.tdims{i}, 1);
            end
            tblocks = cell(nblocks, 1);
            
            offset = 0;
            p = rankrange(b.rank);
            for i = 1:length(b.var)
                rowsz = b.rowsizes{i};
                colsz = b.colsizes{i};
                for k = 1:length(colsz) - 1
                    for j = 1:length(rowsz) - 1
                        offset = offset + 1;
                        tblocks{offset} = permute(reshape(...
                            b.var{i}(rowsz(j)+1:rowsz(j+1), colsz(k)+1:colsz(k+1)), ...
                            b.tdims{i}(j + (k-1) * (length(rowsz)-1), p)), ...
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
            
            mblocks = b.var;
            if nargout > 1
                mcharges = b.charge;
            end
        end
        
        function b = fill_matrix_data(b, vars, charges)
            if nargin < 3 || isempty(charges)
                assert(length(vars) == length(b.var), ...
                    'Invalid number of blocks');
                for i = 1:length(b.var)
                    b.var{i} = vars{i};
                end
                return
            end
            [lia, locb] = ismember(charges, b.charge);
            assert(all(lia));
            b.var(locb) = vars;
        end
        
        function b = fill_matrix_fun(b, fun, charges)
            if nargin < 3 || isempty(charges)
                for i = 1:length(b.var)
                    b.var{i} = fun(size(b.var{i}));
                end
            else
                [lia, locb] = ismember(charges, b.charge);
                for i = locb(lia)
                    b.var{i} = fun(size(b.var{i}), b.charge(i));
                end
            end
        end
        
        function b = fill_tensor_data(b, data)
            ctr = 0;
            p = rankrange(b.rank);
            
            for i = 1:length(b.var)
                rowsz = b.rowsizes{i};
                colsz = b.colsizes{i};
                for k = 1:length(colsz) - 1
                    for j = 1:length(rowsz) - 1
                        ctr = ctr + 1;
                        b.var{i}(rowsz(j)+1:rowsz(j+1), colsz(k)+1:colsz(k+1)) = ...
                            reshape(permute(data{ctr}, p), ...
                            rowsz(j+1)-rowsz(j), colsz(k+1)-colsz(k));
                    end
                end
            end
        end
        
        function b = fill_tensor_fun(b, fun, trees)
            
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
            
            for i = 1:length(Y.var)
                Y.var{i} = X.var{i} - Y.var{i};
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
            
            for i = 1:length(Y.var)
                Y.var{i} = X.var{i} + Y.var{i};
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
            
            for i = 1:length(Y.var)
                Y.var{i} = Y.var{i} .* a;
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
            
            for i = 1:length(Y.var)
                Y.var{i} = Y.var{i} ./ a;
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
            
            for i = 1:length(Y.var)
                Y.var{i} = -Y.var{i};
            end
        end
        
        function t = underlyingType(X)
            t = underlyingType(X.var{1});
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
            
            X.rank = fliplr(X.rank);
            for i = 1:length(X.var)
                X.var{i} = X.var{i}';
                ncols = length(X.colsizes{i}) - 1;
                nrows = length(X.rowsizes{i}) - 1;
                X.tdims{i} = reshape(permute(...
                    reshape(fliplr(X.tdims{i}), nrows, ncols, []), ...
                    [2 1 3]), ...
                    ncols * nrows, []);
            end
            
            [X.rowsizes, X.colsizes] = swapvars(X.rowsizes, X.colsizes);
        end
        
        function v = vectorize(X, type)
            if numel(X) > 1
                vs = cell(size(X));
                for i = 1:numel(vs)
                    vs{i} = vectorize(X(i), type);
                end
                v = vertcat(vs{:});
                return
            end
            
            qdims = sqrt(qdim(X.charge));
            switch type
                case 'complex'
                    blocks = cell(size(X.var));
                    for i = 1:length(X.var)
                        blocks{i} = reshape(X.var{i} .* qdims(i), [], 1);
                    end
                    v = vertcat(blocks{:});
                    
                case 'real'
                    blocks = cell(size(X.var));
                    for i = 1:length(X.var)
                        tmp = X.var{i} .* qdims(i);
                        blocks{i} = [reshape(real(tmp), [], 1)
                            reshape(imag(tmp), [], 1)];
                    end
                    v = vertcat(blocks{:});
            end
        end
        
        function X = devectorize(v, X, type)
            ctr = 0;
            for i = 1:numel(X)
                qdims = sqrt(qdim(X(i).charge));
                switch type
                    case 'complex'
                        for j = 1:length(X(i).var)
                            n = numel(X(i).var{j});
                            X(i).var{j} = reshape(v(ctr + (1:n)), size(X(i).var{j})) ./ qdims(j);
                            ctr = ctr + n;
                        end

                    case 'real'
                        for j = 1:length(X(i).var)
                            n = numel(X(i).var{j});
                            sz = size(X(i).var{j});
                            X(i).var{j} = complex(reshape(v(ctr + (1:n)), sz), ...
                                reshape(v(ctr + (n + 1:2 * n)), sz)) ./ qdims(j);
                            ctr = ctr + 2 * n;
                        end
                end
            end
        end
    end
    
    methods
        function assertBlocksizes(X)
            for i = 1:numel(X.var)
                assert(isequal(size(X.var{i}), [X.rowsizes{i}(end) X.colsizes{i}(end)]), ...
                    'kernel:dimerror', 'Wrong size of block');
                rows = length(X.rowsizes{i}) - 1;
                cols = length(X.colsizes{i}) - 1;
                matdims = [prod(X.tdims{i}(:, 1:X.rank(1)), 2) ...
                    prod(X(i).tdims(:, X.rank(1)+1:end), 2)];
                for k = 1:cols
                    for j = 1:rows
                        assert(matdims(j + (k-1) * rows, 1) == X.rowsizes{i}(j+1) - X.rowsizes{i}(j));
                        assert(matdims(j + (k-1) * rows, 2) == X.colsizes{i}(k+1) - X.colsizes{i}(k));
                    end
                end
            end
        end
    end
end
