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
                GL{1} = SparseTensor.zeros(1, size(T(1).O{1}, 2), 1);
                pSpace = space(T(1).O{1}(:,:,:,1), 4);
                GL{1}(1) = insert_onespace(fixedpoint(mps1, 'l_LL'), ...
                    2, ~isdual(pSpace(1)));
            end
            
            for i = 2:size(GL{1}, 2)
                rhs = apply(slice(T, i, 1:i-1), GL{1}(1, 1:i-1, 1));
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GL{1}(i) = rhs;
                elseif iseye(T, i)
                    fp_left  = insert_onespace(fixedpoint(mps1, 'l_LL'), ...
                        2, isdual(space(rhs, 2)));
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
                GR{1} = SparseTensor.zeros(1, N, 1);
                pSpace = space(T(1).O{1}(:, end, :, :), 2);
                GR{1}(1, N, 1) = insert_onespace(fixedpoint(mps1, 'r_RR'), ...
                    2, isdual(pSpace(end)));
            end
            
            for i = N-1:-1:1
                rhs = apply(slice(T, i, i+1:N), GR{1}(1, i+1:N, 1));
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GR{1}(i) = rhs;
                elseif iseye(T, i)
                    fp_left  = insert_onespace(fixedpoint(mps1, 'l_RR'), ...
                        2, ~isdual(space(rhs, 2)));
                    fp_right = insert_onespace(fixedpoint(mps1, 'r_RR'), ...
                        2, isdual(space(rhs, 2)));
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
    end
    
    methods (Static)
        function mpo = Ising(J, h, kwargs)
            arguments
                J = 1
                h = 1
                kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'Z2'})} = 'Z1'
            end
            
            sigma_x = [0 1; 1 0];
            sigma_z = [1 0; 0 -1];
            
            if strcmp(kwargs.Symmetry, 'Z1')
                pSpace = CartesianSpace.new(2);
                S = Tensor([one(pSpace) pSpace], [pSpace one(pSpace)]);
                Sx = fill_matrix(S, sigma_x);
                Sz = fill_matrix(S, sigma_z);
                
                O = MpoTensor.zeros(3, 1, 3, 1);
                O(1, 1, 1, 1) = 1;
                O(3, 1, 3, 1) = 1;
                O(1, 1, 2, 1) = -J * Sx;
                O(2, 1, 3, 1) = Sx;
                O(1, 1, 3, 1) = (-J * h) * Sz;
                
            else
                pSpace = GradedSpace.new(Z2(0, 1), [1 1], false);
                vSpace = GradedSpace.new(Z2(1), 1, false);
                trivSpace = one(pSpace);
                
                Sx_l = fill_matrix(Tensor([trivSpace pSpace], [pSpace vSpace]), {1 1});
                Sx_r = fill_matrix(Tensor([vSpace pSpace], [pSpace trivSpace]), {1 1});
                Sz = fill_matrix(Tensor([trivSpace pSpace], [pSpace trivSpace]), {1 -1});
                
                O = MpoTensor.zeros(3, 1, 3, 1);
                O(1, 1, 1, 1) = 1;
                O(3, 1, 3, 1) = 1;
                O(1, 1, 2, 1) = -J * Sx_l;
                O(2, 1, 3, 1) = Sx_r;
                O(1, 1, 3, 1) = (-J * h) * Sz;
            end
            
            mpo = InfJMpo(O);
        end
        
        function mpo = Heisenberg(J, h, kwargs)
            arguments
                J = 1
                h = 0
                kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'U1', 'SU2'})} = 'Z1'
                kwargs.Spin = SU2(3)
            end
            
            switch kwargs.Symmetry
                case 'SU2'
                    assert(isscalar(J) || all(J == J(1)), ...
                        'Different spin couplings not invariant under SU2');
                    assert(h == 0, 'Magnetic field not invariant under SU2');
                    
                    pSpace = GradedSpace.new(kwargs.Spin, 1, false);
                    aSpace = GradedSpace.new(SU2(3), 1, false);
                    tSpace = one(aSpace);
                    
                    s = spin(kwargs.Spin);
                    L = Tensor.ones([tSpace pSpace], [pSpace aSpace]);
                    L = L * (-J(1) * (s^2 + s));
                    R = Tensor.ones([aSpace pSpace], [pSpace tSpace]);
                    
                    
                    O = MpoTensor.zeros(3, 1, 3, 1);
                    O(1, 1, 1, 1) = 1;
                    O(3, 1, 3, 1) = 1;
                    O(1, 1, 2, 1) = L;
                    O(2, 1, 3, 1) = R;
                    
                otherwise
                    error('TBA');
            end
            mpo = InfJMpo(O);
        end
    end
end
