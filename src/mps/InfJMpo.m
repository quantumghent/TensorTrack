classdef InfJMpo < InfMpo
    % Infinite Mpo with a Jordan block structure
    
    methods
        function mpo = InfJMpo(varargin)
            mpo@InfMpo(varargin{:});
            if nargin > 0
                assert(istriu(mpo.O{1}));
                assert(iseye(mpo.O{1}(1, 1, 1, 1)) && iseye(mpo.O{1}(end, 1, end, 1)));
                assert(isconnected(mpo));
            end
        end
        
        function bool = isconnected(mpo)
            bool = true;
        end
        
        function [GL, lambda] = leftenvironment(mpo, mps1, mps2, GL, linopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GL = cell(1, period(mps1))
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            linkwargs = namedargs2cell(linopts);
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'LL');
            
            if isempty(GL) || isempty(GL{1})
                GL = cell(1, period(mps1));
                GL{1} = SparseTensor.zeros(domain(T), []);
                pSpace = space(T(1).O{1}(:,:,:,1), 4);
                fp1 = insert_onespace(fixedpoint(mps1, 'l_LL'), ...
                    2, ~isdual(pSpace(1)));
                GL{1}(1) = repartition(fp1, [nspaces(fp1) 0]);
            end
            
            for i = 2:size(GL{1}, 2)
                rhs = apply(slice(T, i, 1:i-1), GL{1}(1, 1:i-1, 1));
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GL{1}(i) = rhs;
                elseif iseye(T, i)
                    fp_left  = repartition(insert_onespace(fixedpoint(mps1, 'l_LL'), ...
                        2, isdual(space(rhs, 2))), rank(rhs));
                    fp_right = insert_onespace(fixedpoint(mps1, 'r_LL'), ...
                        2, ~isdual(space(rhs, 2)));
                    lambda = contract(rhs, 1:3, fp_right, 3:-1:1);
                    
                    rhs = rhs - lambda * fp_left;
                    [GL{1}(i), ~] = linsolve(@(x) x - apply(Tdiag, x), rhs, GL{1}(i), ...
                        linkwargs{:});
                    GL{1}(i) = GL{1}(i) - ...
                        contract(GL{1}(i), 1:3, fp_right, 3:-1:1) * fp_left;
                else
                    [GL{1}(i), ~] = linsolve(@(x) x - apply(Tdiag, x), rhs, GL{1}(i), ...
                        linkwargs{:});
                end
            end
            
            if nnz(GL{1}) == numel(GL{1})
                GL{1} = full(GL{1});
            end
            
            for w = 1:period(mps1)-1
                T = transfermatrix(mpo, mps1, mps2, w, 'Type', 'LL');
                GL{next(w, period(mps1))} = apply(T, GL{w});
            end
        end
        
        function [GR, lambda] = rightenvironment(mpo, mps1, mps2, GR, linopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GR = cell(1, period(mps1))
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            linkwargs = namedargs2cell(linopts);
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'RR').';
            N = size(T(1).O{1}, 2);
            
            if isempty(GR) || isempty(GR{1})
                GR = cell(1, period(mps1));
                GR{1} = SparseTensor.zeros(domain(T), []);
                pSpace = space(T(1).O{1}(:, end, :, :), 2);
                fp1 = insert_onespace(fixedpoint(mps1, 'r_RR'), ...
                    2, isdual(pSpace(end)));
                GR{1}(1, N, 1) = repartition(fp1, [nspaces(fp1) 0]);
            end
            
            for i = N-1:-1:1
                rhs = apply(slice(T, i, i+1:N), GR{1}(1, i+1:N, 1));
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GR{1}(i) = rhs;
                elseif iseye(T, i)
                    fp_left  = insert_onespace(fixedpoint(mps1, 'l_RR'), ...
                        2, ~isdual(space(rhs, 2)));
                    fp_right = repartition(insert_onespace(fixedpoint(mps1, 'r_RR'), ...
                        2, isdual(space(rhs, 2))), rank(rhs));
                    lambda = contract(rhs, 1:3, fp_left, 3:-1:1);
                    
                    rhs = rhs - lambda * fp_right;
                    [GR{1}(i), ~] = ...
                        linsolve(@(x) x - apply(Tdiag, x), rhs, GR{1}(i), linkwargs{:});
                    
                    GR{1}(i) = GR{1}(i) - ...
                        contract(GR{1}(i), 1:3, fp_left, 3:-1:1) * fp_right;
                else
                    [GR{1}(i), ~] = linsolve(@(x) x - apply(Tdiag, x), rhs, GR{1}(i), ...
                        linkwargs{:});
                end
            end
            
            if nnz(GR{1}) == numel(GR{1})
                GR{1} = full(GR{1});
            end
            
            for w = period(mps1):-1:2
                T = transfermatrix(mpo, mps1, mps2, w, 'Type', 'RR').';
                GR{w} = apply(T, GR{next(w, period(mps1))});
            end
        end
        
        function [GL, GR, lambda] = environments(mpo, mps1, mps2, GL, GR, linopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GL = cell(1, period(mps1))
                GR = cell(1, period(mps1))
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            kwargs = namedargs2cell(linopts);
            [GL, lambdaL] = leftenvironment(mpo, mps1, mps2, GL, kwargs{:});
            [GR, lambdaR] = rightenvironment(mpo, mps1, mps2, GR, kwargs{:});
            lambda = (lambdaL + lambdaR) / 2;
            if abs(lambdaL - lambdaR)/abs(lambda) > eps(lambda)^(1/3)
                warning('lambdas disagree');
            end
        end
        
        function mpo = horzcat(varargin)
            Os = cellfun(@(x) x.O, varargin, 'UniformOutput', false);
            mpo = InfJMpo([Os{:}]);
        end
        
        function mpo = plus(a, b)
            if isa(a, 'InfJMpo') && isnumeric(b)
                if period(a) > 1 && isscalar(b)
                    b = repmat(b, 1, period(a));
                end
                
                for i = 1:period(a)
                    a.O{i}(1, 1, end, 1) = a.O{i}(1, 1, end, 1) + b(i);
                end
                
            elseif isnumeric(a) && isa(b, 'InfJMpo')
                mpo = b + a;
            end
        end
        
        function mpo = mtimes(mpo, b)
            if isnumeric(mpo) || isnumeric(b)
                mpo = mpo .* b;
                return
            end
            
            mpo = [mpo; b];
        end
        
        function finitempo = open_boundary_conditions(mpo, L)
            Os = repmat(mpo.O, 1, L);
            
            Os{1}   = Os{1}(1, :, :, :);
            Os{end} = Os{end}(:, :, end, :);
            
            rspace = subspaces(rightvspace(mpo, period(mpo)), size(Os{end}, 3));
            rightedge = MpsTensor(Tensor.eye([one(rspace) rspace'], one(rspace)));
            
            lspace = subspaces(leftvspace(mpo, 1), size(Os{1}, 1));
            leftedge = MpsTensor(Tensor.eye([one(lspace) lspace], one(lspace))');
            
            finitempo = FiniteMpo(leftedge, Os, rightedge);
        end
    end
end
