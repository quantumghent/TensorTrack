classdef InfMpo
    % Infinite translation invariant matrix product operator.
    
    properties
        O
    end
    
    methods
        function mpo = InfMpo(varargin)
            if nargin == 0, return; end
            
            if nargin == 1
                O = varargin{1};
                
                if isa(O, 'InfMpo')
                    mpo.O = O.O;
                    
                elseif isa(O, 'MpoTensor')
                    mpo.O = {O};
                
                elseif isa(O, 'AbstractTensor')
                    mpo.O = arrayfun(@MpoTensor, O, 'UniformOutput', false);
                    
                elseif iscell(O)
                    mpo.O = O;
                end
                return
            end
        end
    end
    
    methods
        function p = period(mpo)
            p = size(mpo.O, 2);
        end
        
        function d = depth(mpo)
            d = size(mpo.O, 1);
        end
        
        function mpo = block(mpo)
            if depth(mpo) == 1, return; end
            O_ = mpo.O(1, :);
            for d = 2:depth(mpo)
                for w = period(mpo):-1:1
                    vspaces = [rightvspace(O_{w})' rightvspace(mpo.O{d, w})'];
                    fuser(w) = Tensor.eye(vspaces, prod(vspaces));
                end
                
                for w = 1:period(mpo)
                    O_{w} = MpoTensor(contract(O_{w}, [1 2 5 -4], mpo.O{d, w}, [3 -2 4 2], ...
                        fuser(prev(w, period(mpo)))', [-1 3 1], fuser(w), [5 4 -3], ...
                        'Rank', [2 2]));
                end
            end
            mpo.O = O_;
        end

        function s = pspace(mpo, x)
            s = pspace(mpo.O{x});
        end
        
        function mpo = horzcat(mpo, varargin)
            Os = cellfun(@(x) x.O, varargin, 'UniformOutput', false);
            mpo.O = horzcat(mpo.O, Os{:});
        end
        
        function mpo = vertcat(mpo, varargin)
            Os = cellfun(@(x) x.O, varargin, 'UniformOutput', false);
            mpo.O = vertcat(mpo.O, Os{:});
        end
        
        function mpo = slice(mpo, d)
            mpo.O = mpo.O(d, :);
        end
        
        function mps = initialize_mps(mpo, vspaces)
            arguments
                mpo
            end
            arguments (Repeating)
                vspaces
            end
            
            pspaces = arrayfun(@(x) subspaces(pspace(mpo, x)), 1:period(mpo), 'UniformOutput', false);
            args = [pspaces; vspaces];
            mps = UniformMps.randnc(args{:});
        end
    end
    
    
    %% Environments
    methods
        function [GL, lambda] = leftenvironment(mpo, mps1, mps2, GL, eigopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GL = cell(1, period(mps1))
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
                eigopts.Verbosity = Verbosity.warn
            end
            N = period(mps1);
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'LL');
            [GL{1}, lambda] = eigsolve(T, GL{1}, 1, 'largestabs', ...
                'KrylovDim', eigopts.KrylovDim, 'Tol', eigopts.Tol, ...
                'ReOrth', eigopts.ReOrth, 'MaxIter', eigopts.MaxIter);
            
            for i = 2:N
                T = transfermatrix(mpo, mps1, mps2, i-1, 'Type', 'LL');
                GL{i} = apply(T, GL{i-1}) ./ lambda^(1/N);
            end
        end
        
        function [GR, lambda] = rightenvironment(mpo, mps1, mps2, GR, eigopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GR = cell(1, period(mps1))
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
                eigopts.Verbosity = Verbosity.warn
            end
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'RR').';
            [GR{1}, lambda] = eigsolve(T, GR{1}, 1, 'largestabs', ...
                'KrylovDim', eigopts.KrylovDim, 'Tol', eigopts.Tol, ...
                'ReOrth', eigopts.ReOrth, 'MaxIter', eigopts.MaxIter);
            N = period(mps1);
            for i = N:-1:2
                T = transfermatrix(mpo, mps1, mps2, i, 'Type', 'RR').';
                GR{i} = apply(T, GR{next(i, N)}) ./ lambda^(1/N);
            end
        end
        
        function [GL, GR, lambda] = environments(mpo, mps1, mps2, GL, GR, eigopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GL = cell(1, period(mps1))
                GR = cell(1, period(mps1))
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
                eigopts.Verbosity = Verbosity.warn
            end
            
            kwargs = namedargs2cell(eigopts);
            [GL, lambdaL] = leftenvironment(mpo, mps1, mps2, GL, kwargs{:});
            [GR, lambdaR] = rightenvironment(mpo, mps1, mps2, GR, kwargs{:});
            lambda = (lambdaL + lambdaR) / 2;
            if abs(lambdaL - lambdaR)/abs(lambda) > eps(lambda)^(1/3)
                warning('lambdas disagree');
            end
            
            for w = 1:period(mps1)
                overlap = sqrt(contract(GL{w}, [1 3 2], ...
                    mps1.C(prev(w, period(mps1))), [2 4], ...
                    mps2.C(prev(w, period(mps1)))', [5 1], ...
                    GR{w}, [4 3 5]));
                GL{w} = GL{w} ./ overlap;
                GR{w} = GR{w} ./ overlap;
            end
        end
    end
    
    
    %% Derived operators
    methods
        function T = transfermatrix(mpo, mps1, mps2, sites, kwargs)
            arguments
                mpo
                mps1
                mps2 = mps1
                sites = 1:period(mps1)
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
            O = mpo.O(sites);                                                               %#ok<PROPLC>
            T = FiniteMpo.mps_channel_operator(A1, O, A2);                                  %#ok<PROPLC>
        end
        
        function H = AC_hamiltonian(mpo, mps, GL, GR, sites)
            arguments
                mpo
                mps
                GL = fixedpoint(transfermatrix(mpo, mps, 'Type', 'LL'))
                GR = fixedpoint(transfermatrix(mpo, mps, 'Type', 'RR').')
                sites = 1:period(mps)
            end
            
            H = cell(1, length(sites));
            for i = 1:length(sites)
                for d = depth(mpo):-1:1
                    gl = twistdual(GL{d, sites(i)}, 1);
                    gr = GR{d, next(sites(i), period(mps))};
                    gr = twistdual(gr, nspaces(gr));
                    H{i}(d, 1) = FiniteMpo(gl, mpo.O(d, sites(i)), gr);
                end
            end
        end
        
        function H = AC2_hamiltonian(mpo, mps, GL, GR, sites)
            arguments
                mpo
                mps
                GL = fixedpoint(transfermatrix(mpo, mps, 'Type', 'LL'))
                GR = fixedpoint(transfermatrix(mpo, mps, 'Type', 'RR').')
                sites = 1:period(mps)
            end
            
            H = cell(1, length(sites));
            for i = 1:length(sites)
                for d = depth(mpo):-1:1
                    gl = GL{d, sites(i)};
                    gl = twistdual(gl, 1);
                    gr = GR{d, mod1(sites(i) + 2, period(mps))};
                    gr = twistdual(gr, nspaces(gr));
                    H{i}(d, 1) = FiniteMpo(gl, mpo.O(d, mod1(sites(i) + [0 1], period(mps))), gr);
                end
            end
        end
        
        function H = C_hamiltonian(mpo, mps, GL, GR, sites)
            arguments
                mpo
                mps
                GL = fixedpoint(transfermatrix(mpo, mps, 'Type', 'LL'))
                GR = fixedpoint(transfermatrix(mpo, mps, 'Type', 'RR').')
                sites = 1:period(mps)
            end
            
            H = cell(1, length(sites));
            for i = 1:length(sites)
                for d = depth(mpo):-1:1
                    gl = GL{d, next(sites(i), period(mps))};
                    for j = 1:numel(gl)
                        if nnz(gl(j)) ~= 0
                            gl(j) = twist(gl(j), find(isdual(space(gl(j), 1))));
                        end
                    end
                    gr = GR{d, next(sites(i), period(mps))};
                    for j = 1:numel(gr)
                        if nnz(gr(j)) ~= 0
                            gr(j) = twist(gr(j), find(isdual(space(gr(j), nspaces(gr(j))))) + ...
                                nspaces(gr(j)) - 1);
                        end
                    end
                    H{i}(d, 1) = FiniteMpo(gl, {}, gr);
                end
            end
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
            
            mpo = InfMpo(O);
        end
        
        function mpo = fDimer()
            pspace = GradedSpace.new(fZ2(0, 1), [1 1], false);
            O = Tensor([pspace pspace], [pspace pspace]);
            O1 = fill_tensor(O, @(~, f) ~any(f.uncoupled) || ...
                (f.uncoupled(2) && sum(f.uncoupled) == 2));
            O2 = fill_tensor(O, @(~, f) ~any(f.uncoupled) || ...
                (f.uncoupled(4) && sum(f.uncoupled) == 2));
            
            mpo = InfMpo([O1; O2]);
        end
    end
end
