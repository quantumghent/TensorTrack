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
        
        function [GL, lambda] = leftenvironment(mpo, mps1, mps2, GL, eigopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GL = []
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'LL');
            
            if isempty(GL)
                GL = SparseTensor.zeros(1, size(T.O{1}, 2), 1);
                pSpace = pspace(mpo.O{1});
                GL(1) = insert_onespace(fixedpoint(mps1, 'l_LL'), ...
                    2, ~isdual(pSpace(1)));
            end
            
            for i = 2:size(mpo.O{1}, 1)
                rhs = apply(slice(T, i, 1:i-1), GL(1, 1:i-1, 1));
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GL(i) = rhs;
                elseif iseye(Tdiag)
                    fp_left  = insert_onespace(fixedpoint(mps1, 'l_LL'), ...
                        2, isdual(space(rhs, 2)));
                    fp_right = insert_onespace(fixedpoint(mps1, 'r_LL'), ...
                        2, ~isdual(space(rhs, 2)));
                    lambda = contract(rhs, 1:3, fp_right, 3:-1:1);
                    
                    rhs = rhs - lambda * fp_left;
                    GL(i) = linsolve(@(x) x - apply(Tdiag, x), rhs, ...
                        GL(i));
                    GL(i) = GL(i) - contract(GL(i), 1:3, fp_right, 3:-1:1) * fp_left;
                else
                    GL(i) = linsolve(@(x) x - apply(Tdiag, x), rhs, ...
                        GL(i));
                end
            end
        end
        
        function [GR, lambda] = rightenvironment(mpo, mps1, mps2, GR, eigopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GR = []
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'RR').';
            N = size(T.O{1}, 2);
            if isempty(GR)
                GR = SparseTensor.zeros(1, N, 1);
                pSpace = pspace(mpo.O{1});
                GR(1, N, 1) = insert_onespace(fixedpoint(mps1, 'r_RR'), ...
                    2, isdual(pSpace(end)));
            end
            
            for i = N-1:-1:1
                rhs = apply(slice(T, i, i+1:N), GR(1, i+1:N, 1));
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GR(i) = rhs;
                elseif iseye(Tdiag)
                    fp_left  = insert_onespace(fixedpoint(mps1, 'l_RR'), ...
                        2, ~isdual(space(rhs, 2)));
                    fp_right = insert_onespace(fixedpoint(mps1, 'r_RR'), ...
                        2, isdual(space(rhs, 2)));
                    lambda = contract(rhs, 1:3, fp_left, 3:-1:1);
                    
                    rhs = rhs - lambda * fp_right;
                    GR(i) = linsolve(@(x) x - apply(Tdiag, x), rhs, ...
                        GR(i));
                    GR(i) = GR(i) - contract(GR(i), 1:3, fp_left, 3:-1:1) * fp_right;
                else
                    GR(i) = linsolve(@(x) x - apply(Tdiag, x), rhs, ...
                        GR(i));
                end
            end
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
    end
end
