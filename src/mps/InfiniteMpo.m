classdef InfiniteMpo
    % Infinite translation invariant matrix product operator.
    
    properties
        O
    end
    
    methods
        function mpo = InfiniteMpo(varargin)
            if nargin == 0, return; end
            
            if nargin == 1
                O = varargin{1};
                
                if isa(O, 'InfiniteMpo')
                    for i = numel(O):1:1
                        mpo(i, 1).O = O(i, 1).O;
                    end
                    mpo = reshape(mpo, size(O));
                    
                elseif isa(O, 'MpoTensor')
                    mpo.O = {O};
                
                elseif isa(O, 'AbstractTensor')
                    for i = height(O):-1:1
                        mpo(i, 1).O = arrayfun(@MpoTensor, O(i, :), 'UniformOutput', false);
                    end
                    
                elseif iscell(O)
                    mpo.O = O;
                end
                return
            end
        end
    end
    
    methods
        function p = period(mpo)
            p = length(mpo(1).O);
        end
        
        function d = depth(mpo)
            d = length(mpo);
        end
        
        function mpo = block(mpo)
            if depth(mpo) == 1, return; end
            O_ = mpo(1, 1).O;
            for d = 2:depth(mpo)
                
                for w = period(mpo):-1:1
                    vspaces = [rightvspace(O_{w})' rightvspace(mpo(d, 1).O{w})'];
                    fuser(w) = Tensor.eye(vspaces, prod(vspaces));
                end
                
                for w = 1:period(mpo)
                    O_{w} = MpoTensor(contract(O_{w}, [1 2 5 -4], mpo(d, 1).O{w}, [3 -2 4 2], ...
                        fuser(prev(w, period(mpo)))', [-1 3 1], fuser(w), [5 4 -3], ...
                        'Rank', [2 2]));
                end
            end
            mpo = mpo(1, 1);
            mpo.O = O_;
        end

        function [GL, lambda] = leftenvironment(mpo, mps1, mps2, GL, eigopts)
            arguments
                mpo
                mps1
                mps2
                GL = []
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'LL');
            [GL, lambda] = eigsolve(T, GL, 1, 'largestabs', ...
                'KrylovDim', eigopts.KrylovDim, 'Tol', eigopts.Tol, ...
                'ReOrth', eigopts.ReOrth, 'MaxIter', eigopts.MaxIter);
        end
        
        function [GR, lambda] = rightenvironment(mpo, mps1, mps2, GR, eigopts)
            arguments
                mpo
                mps1
                mps2
                GR = []
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'RR').';
            
            [GR, lambda] = eigsolve(T, GR, 1, 'largestabs', ...
                'KrylovDim', eigopts.KrylovDim, 'Tol', eigopts.Tol, ...
                'ReOrth', eigopts.ReOrth, 'MaxIter', eigopts.MaxIter);
        end
        
        function [GL, GR, lambda] = environments(mpo, mps1, mps2, GL, GR, eigopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GL = []
                GR = []
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            kwargs = [fieldnames(eigopts).'; struct2cell(eigopts).'];
            [GL, lambdaL] = leftenvironment(mpo, mps1, mps2, GL, kwargs{:});
            [GR, lambdaR] = rightenvironment(mpo, mps1, mps2, GR, kwargs{:});
            lambda = (lambdaL + lambdaR) / 2;
            if abs(lambdaL - lambdaR)/abs(lambda) > eps(lambda)^(1/3)
                warning('lambdas disagree');
            end
            
            overlap = sqrt(contract(GL, [1 3 2], mps1.C, [2 4], mps2.C', [5 1], GR, [4 3 5]));
            GL = GL / overlap;
            GR = GR / overlap;
        end
    end
    %% Derived operators
    methods
        function T = transfermatrix(mpo, mps1, mps2, sites, kwargs)
            arguments
                mpo
                mps1
                mps2 = mps1
                sites = 1:period(mpo)
                kwargs.Type {mustBeMember(kwargs.Type, {'LL' 'LR' 'RL' 'RR'})} = 'RR'
            end
            
            assert(all(diff(sites) == 1), 'sites must be neighbouring and increasing.');
            
            if kwargs.Type(1) == 'L'
                A1 = mps1.AL(sites);
            else
                A1 = mps1.AR(sites);
            end
            if kwargs.Type(2) == 'L'
                A2 = mps2.AL(sites);
            else
                A2 = mps2.AR(sites);
            end
            
            A2 = twist(A2, 1 + find(isdual(space(A1, 2:nspaces(A1)-1))));
            T = FiniteMpo(A2', {rot90(mpo.O{1})}, A1);
            
%             T = transfermatrix(mps1, mps2, sites, 'Type', kwargs.Type);
%             T.O = {rot90(mpo.O{1})};
        end
        
        function H = AC_hamiltonian(mpo, mps, GL, GR)
            arguments
                mpo
                mps
                GL = fixedpoint(transfermatrix(mpo, mps, 'Type', 'LL'))
                GR = fixedpoint(transfermatrix(mpo, mps, 'Type', 'RR').')
            end
            GR = twist(GR, find(isdual(space(GR, nspaces(GR)))) + nspaces(GR)-1);
            GL = twist(GL, find(isdual(space(GL, 1))));
            H = FiniteMpo(GL, mpo.O, GR);
        end
        
        function H = C_hamiltonian(mpo, mps, GL, GR)
            arguments
                mpo
                mps
                GL = fixedpoint(transfermatrix(mpo, mps, 'Type', 'LL'))
                GR = fixedpoint(transfermatrix(mpo, mps, 'Type', 'RR').')
            end
            GR = twist(GR, find(isdual(space(GR, nspaces(GR)))) + nspaces(GR)-1);
            GL = twist(GL, find(isdual(space(GL, 1))));
            H = FiniteMpo(GL, {}, GR);
        end
    end
    
    methods (Static)
        function mpo = Ising(beta, kwargs)
            arguments
                beta = log(1 + sqrt(2)) / 2;
                kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'Z2'})} = 'Z1'
            end
            
            if strcmp(kwargs.Symmetry, 'Z1')
                t = [exp(beta) exp(-beta); exp(-beta) exp(beta)];
                [v, d] = eig(t);
                t = v * sqrt(d) * v;
                
                o = zeros(2, 2, 2, 2);
                o(1, 1, 1, 1) = 1;
                o(2, 2, 2, 2) = 1;
                
                o = contract(o, 1:4, t, [-1 1], t, [-2 2], t, [-3 3], t, [-4 4]);
                
                O = fill_tensor(Tensor.zeros([2 2 2 2]), o);
                
            else
                s = GradedSpace.new(Z2(0, 1), [1 1], false);
                O = fill_tensor(Tensor([s s], [s s]), ...
                    @(~, f) 2 * sqrt(prod(logical(f.uncoupled) .* sinh(beta) + ...
                    ~logical(f.uncoupled) .* cosh(beta))));
            end
            
            mpo = InfiniteMpo(O);
        end
        
        function mpo = fDimer()
            pspace = GradedSpace.new(fZ2(0, 1), [1 1], false);
            O = Tensor([pspace pspace], [pspace pspace]);
            O1 = fill_tensor(O, @(~, f) ~any(f.uncoupled) || ...
                (f.uncoupled(2) && sum(f.uncoupled) == 2));
            O2 = fill_tensor(O, @(~, f) ~any(f.uncoupled) || ...
                (f.uncoupled(4) && sum(f.uncoupled) == 2));
            
            mpo = InfiniteMpo([O1; O2]);
        end
    end
end
