classdef TrivialBlock < AbstractBlock
    % Data structure for tensors without symmetry.
    %
    % .. todo::
    %   Document properties and behavior.
    
    properties
        var
        tdims
        rank
    end
    
    methods
        function b = TrivialBlock(codomain, domain)
            if nargin == 0, return; end
            
            b.rank = [length(codomain) length(domain)];
            b.tdims = degeneracies([codomain, domain(length(domain):-1:1)]);
            b.var = uninit(b.tdims);
        end
        
        function Y = axpby(a, X, b, Y, p, ~)
            
            %%% Special case 1: scalar multiplication
            if a == 0
                if b == 1
                    return;
                end
                
                if b == -1
                    Y = -Y;
                    return
                end
                
                Y = Y .* b;
                return
            end
            
            
            %% Special case 2: addition without permutation
            if nargin == 4 || isempty(p) || all(p == 1:length(p)) && ...
                    all(Y.rank == X.rank)
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
                    Y.var = Y.var .* b + X.var;
                    return
                end
                
                if b == 1   % a ~= 1
                    Y.var = Y.var + X.var .* a;
                    return
                end
                
                % a ~= [0 1], b ~= [0 1]
                Y.var = Y.var .* b + X.var .* a;
                return
            end
            
            
            %% General case: addition with permutation
            if a == 1
                if b == 1
                    Y.var = Y.var + permute(X.var, p);
                    return
                end
                
                if b == 0
                    Y.var = permute(X.var, p);
                    return
                end
                
                Y.var = Y.var .* b + permute(X.var, p);
                return
            end
            
            if b == 1
                Y.var = Y.var + permute(X.var, p) .* a;
                return
            end
            if nargin == 0, return; end
            
            Y.var = Y.var .* b + permute(X.var, p) .* a;
        end
        
        function C = mul(C, A, B, a, b)
            if nargin > 3 && ~isempty(a)
                A = A .* a;
            end
            
            if nargin > 4 && ~isempty(b)
                B = B .* b;
            end
            
            C.var = tensorprod(A.var, B.var, ...
                A.rank(1) + (1:A.rank(2)), B.rank(1):-1:1, ...
                'NumDimensionsA', sum(A.rank));
        end
        
        function A = ctranspose(A)
            A.rank = flip(A.rank);
            A.tdims = flip(A.tdims);
            A.var = conj(permute(A.var, length(A.tdims):-1:1));
        end
        
        function A = times(A, a)
            if a == 1, return; end
            if a == -1, A = -A; return; end
            A.var = A.var .* a;
        end
        
        function Y = rdivide(Y, a)
            
            if a == 1, return; end
            if a == -1, Y = -Y; return; end
            
            Y.var = Y.var ./ a;
        end
        
        function tblocks = tensorblocks(b)
            tblocks = {b.var};
        end
        
        function [mblocks, mcharges] = matrixblocks(b)
            mblocks = {reshape(permute(b.var, rankrange(b.rank)), ...
                prod(b.tdims(1:b.rank(1)), 2), ...
                prod(b.tdims((1:b.rank(2)) + b.rank(1)), 2))};
            
            if nargout > 1
                mcharges = Z1;
            end
        end
        
        function Y = minus(X, Y)
            Y.var = X.var - Y.var;
        end
        
        function Y = plus(X, Y)
            Y.var = X.var + Y.var;
        end
        
        function A = uplus(A)
        end
        
        function Y = uminus(Y)
            Y.var = -Y.var;
        end
        
        function Y = fill_matrix_data(Y, dat, ~)
            rry = rankrange(Y.rank);
            Y.var = permute(reshape(dat{1}, Y.tdims(rry)), rry);
        end
        
        function Y = fill_matrix_fun(Y, fun, ~)
            rry = rankrange(Y.rank);
            mdims = [prod(Y.tdims(1:Y.rank(1)), 2), ...
                prod(Y.tdims(Y.rank(1) + (1:Y.rank(2))), 2)];
            Y.var = permute(reshape(fun(mdims), Y.tdims(rry)), rry);
        end
        
        function Y = fill_tensor_data(Y, dat, ~)
            Y.var = dat{1};
        end
        
        function Y = fill_tensor_fun(Y, fun, ~)
            Y.var = fun(Y.tdims);
        end
        
        function v = vectorize(X, type)
            if numel(X) > 1
                vs = cell(size(X));
                for i = 1:numel(X)
                    vs{i} = vectorize(X(i), type);
                end
                v = vertcat(vs{:});
                return
            end
            
            switch type
                case 'complex'
                    v = X.var(:);
            
                case 'real'
                    v = [real(X.var(:)); imag(X.var(:))];
            end
        end
        
        function X = devectorize(v, X, type)
            ctr = 0;
            for i = 1:numel(X)
                switch type
                    case 'complex'
                        X(i).var = reshape(v(ctr + (1:numel(X(i).var))), size(X(i).var));
                        ctr = ctr + numel(X(i).var);
                    case 'real'
                        v_ = v(ctr + 1:2*numel(X(i).var));
                        m = size(v_, 1) / 2;
                        sz = size(X(i).var);
                        X(i).var = complex(reshape(v_(1:m), sz), reshape(v_(m + 1:2 * m), sz));
                        ctr = ctr + numel(X(i).var) * 2;
                end
            end
        end
        
        function type = underlyingType(X)
            type = underlyingType(X.var);
        end
    end
end
