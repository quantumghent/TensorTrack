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

        function s = pspace(mpo, w)
            s = pspace(mpo.O{w});
        end
        
        function s = leftvspace(mpo, w)
            s = leftvspace(mpo.O{w});
        end
        
        function s = rightvspace(mpo, w)
            s = rightvspace(mpo.O{w});
        end
        
        function mpo = repmat(mpo, varargin)
            mpo.O = repmat(mpo.O, varargin{:});
        end
        
        function mpo = horzcat(mpo, varargin)
            Os = cellfun(@(x) x.O, varargin, 'UniformOutput', false);
            mpo.O = horzcat(mpo.O, Os{:});
        end
        
        function mpo = vertcat(mpo, varargin)
            Os = cellfun(@(x) x.O, varargin, 'UniformOutput', false);
            mpo.O = vertcat(mpo.O, Os{:});
        end
        
        function mpo = mrdivide(mpo, scalar)
            assert(depth(mpo) == 1);
            if isscalar(scalar)
                scalar = scalar .^ (ones(size(mpo.O)) ./ length(mpo));
            end
            for i = 1:length(mpo.O)
                mpo.O{i} = mpo.O{i} / scalar(i);
            end
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
            if length(vspaces) ~= length(pspaces)
                assert(length(vspaces) == 1, 'ArgError', ...
                    'Invalid length of input spaces')
                vspaces = repmat(vspaces, size(pspaces));
            end
            args = [pspaces; vspaces];
            mps = UniformMps.randnc(args{:});
        end
    end
    
    
    %% Environments
    methods
        function fp = fixedpoint(operator, state, type, w)
            arguments
                operator
                state
                type
                w = 1
            end
            
            fp = fixedpoint(state, type(1:4), w);
            
            % add leg to fit operator
            switch type(1)
                case 'l'
                    fp = insert_onespace(fp, 2, ~isdual(leftvspace(operator, w)));
                case 'r'
                    fp = insert_onespace(fp, 2, ~isdual(rightvspace(operator, w)));
                otherwise
                    error('invalid fixedpoint type (%s)', type);
            end
            
            % add leg to fit quasiparticle auxiliary leg
            if isa(state, 'InfQP')
                switch type(6)
                    case '0'
                        dual = isdual(auxspace(state, w));
                    case '1'
                        dual = ~isdual(auxspace(state, w));
                    otherwise
                        error('invalid type (%s)', type);
                end
                fp = MpsTensor(insert_onespace(fp, nspaces(fp) + 1, dual), 1);
            end
        end
        
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
            GL = cellfun(@MpsTensor, GL, 'UniformOutput', false);
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
            GR = cellfun(@MpsTensor, GR, 'UniformOutput', false);
        end
        
        function GBL = leftquasienvironment(mpo, qp, GL, GR, GBL, linopts)
            arguments
                mpo
                qp
                GL
                GR
                GBL = cell(1, period(mpo))
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(qp))^(3/4)
            end
            
            assert(period(qp) == period(mpo), 'quasiparticles have different period');
            
            linkwargs = namedargs2cell(linopts);
            expP = exp(-1i * qp.p);
            L = period(mpo);
            
            needsRegularization = istrivial(qp);
            if needsRegularization || true
                fp_left  = fixedpoint(mpo, qp, 'l_RL_0', 1);
                fp_right = fixedpoint(mpo, qp, 'r_RL_1', L);
            end
            
            T = transfermatrix(mpo, qp, qp, 'Type', 'RL');
            TB = transfermatrix(mpo, qp, qp, 'Type', 'BL');
            
            % initialize and precompute GL * TB
            GBL = cell(size(GL));
            GBL{1} = MpsTensor(SparseTensor.zeros(domain(T), auxspace(qp, 1)), 1);
            for w = 1:L
                TB_w = transfermatrix(mpo, qp, qp, w, 'Type', 'BL');
                T_w = transfermatrix(mpo, qp, qp, w, 'Type', 'RL');
                GBL{next(w, L)} = ...
                    (apply(TB_w, GL{w}) + apply(T_w, GBL{w})) * (expP^(1/L));
            end
            
%             rho = GL{1};
%             for pos = 1:L
%                 T_B = transfermatrix(mpo, qp, qp, pos, 'Type', 'BL');
%                 if pos == 1
%                     rho = apply(T_B, GL{pos});
%                 else
%                     T_R = transfermatrix(mpo, qp, qp, pos, 'Type', 'RL');
%                     rho = apply(T_B, GL{pos}) + apply(T_R, rho);
%                 end
%             end
            T = transfermatrix(mpo, qp, qp, 'Type', 'RL');
            TB = transfermatrix(mpo, qp, qp, 'Type', 'BL');
            
            
            rhs = GBL{1};
            if needsRegularization
                lambda = overlap(rhs, fp_right);
                rhs = rhs - lambda * fp_left;
                GBL{1} = linsolve(@(x) x - expP * apply_regularized(T, fp_left, fp_right, x), ...
                    rhs, [], linkwargs{:});
            else
                GBL{1} = expP * linsolve(@(x) x - expP * apply(T, x), ...
                    rhs, [], linkwargs{:});
            end
            
            if nnz(GBL{1}) == numel(GBL{1}.var)
                GBL{1}.var = full(GBL{1}.var);
            end
            for w = 1:L-1
                GBL{next(w, L)} = expP^(1 / L) * ...
                    (apply(TB(w), GL{w}) + apply(T(w), GBL{w}));
            end
        end
        
        function GBR = rightquasienvironment(mpo, qp, GL, GR, GBR, linopts)
            arguments
                mpo
                qp
                GL
                GR
                GBR = cell(1, period(mpo))
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(qp))^(3/4)
            end
            
            assert(period(qp) == period(mpo), 'quasiparticles have different period');
            
            linkwargs = namedargs2cell(linopts);
            
            expP = exp(1i * qp.p);
            
            rho = GR{1};
            for pos = period(mpo):-1:1
                T_B = transfermatrix(mpo, qp, qp, pos, 'Type', 'BR').';
                if pos == period(mpo)
                    rho = apply(T_B, GR{pos});
                else
                    T_L = transfermatrix(mpo, qp, qp, pos, 'Type', 'LR').';
                    rho = apply(T_B, GR{pos}) + apply(T_L, rho);
                end
            end
            
            T_L = transfermatrix(mpo, qp, qp, 1:period(mpo), 'Type', 'LR').';
            if istrivial(qp)
                C = qp.mpsright.C(1);
                FL = insert_onespace(multiplyleft(MpsTensor(GL{1}), C'), ...
                    1, ~isdual(auxspace(qp, 1)));
                FR = insert_onespace(multiplyleft(MpsTensor(GR{1}), C), ...
                    nspaces(GR{1}) + 1, isdual(auxspace(qp, 1)));
                
                rho = rho - overlap(rho, FL) * FR;
                GBR{1} = linsolve(@(x) x - expP * apply_regularized(T_L, FR, FL, x), ...
                    expP * rho, [], linkwargs{:});
            else
                GBR{1} = linsolve(@(x) x - expP * apply(T_L, x), ...
                    expP * rho, [], linkwargs{:});
            end
            
            N = period(mpo);
            GBR{1} = MpsTensor(GBR{1}, 1);
            for i = N:-1:2
                T_L = transfermatrix(mpo, qp, qp, i, 'Type', 'LR').';
                T_B = transfermatrix(mpo, qp, qp, i, 'Type', 'BR').';
                GBR{i} = apply(T_L, GBR{next(i, N)}) + apply(T_B, GR{next(i, N)});
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
                overlap = sqrt(contract(GL{w}, [1, 3:nspaces(GL{w}), 2], ...
                    mps1.C{prev(w, period(mps1))}, [2, nspaces(GL{w})+1], ...
                    mps2.C{prev(w, period(mps1))}', [nspaces(GL{w})+2, 1], ...
                    GR{w}, [nspaces(GL{w})+1, flip(3:nspaces(GL{w})), nspaces(GL{w})+2]));
                GL{w} = GL{w} ./ overlap;
                GR{w} = GR{w} ./ overlap;
            end
        end
    end
    
    
    %% Derived operators
    methods
        function T = transfermatrix(operator, state1, state2, sites, kwargs)
            arguments
                operator
                state1
                state2 = state1
                sites = 1:period(state1)
                kwargs.Type (1,2) char = 'RR'
            end
            
            assert(all(diff(sites) == 1), 'sites must be neighbouring and increasing.');
            assert(length(kwargs.Type) == 2);
            
            % unpack tensors top state
            switch kwargs.Type(1)
                case 'L'
                    top = state1.AL(sites);
                case 'R'
                    top = state1.AR(sites);
                case 'B'
                    top = state1.B(sites);
                otherwise
                    error('invalid transfermatrix top type (%s)', kwargs.Type(1));
            end
            
            % unpack tensors bottom state
            switch kwargs.Type(2)
                case 'L'
                    bot = state2.AL(sites);
                case 'R'
                    bot = state2.AR(sites);
                case 'B'
                    bot = state2.B(sites);
                otherwise
                    error('invalid transfermatrix bot type (%s)', kwargs.Type(2));
            end
            
            % unpack tensors operator
            mid = operator.O(sites);
            
            T = FiniteMpo.mps_channel_operator(top, mid, bot);
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
                    p = 1:nspaces(gl);
                    p(2:(nspaces(gl)-1-gl.alegs)) = flip(p(2:(nspaces(gl)-1-gl.alegs)));
                    gl = tpermute(gl, p, rank(gl));
                    
                    gr = GR{d, next(sites(i), period(mpo))};
                    gr = twistdual(gr, nspaces(gr)-gr.alegs);
                    p = 1:nspaces(gr);
                    p(2:(nspaces(gr)-1-gr.alegs)) = flip(p(2:(nspaces(gr)-1-gr.alegs)));
                    gr = tpermute(gr, p, rank(gr));
                    
                    H{i}(d, 1) = FiniteMpo(gl, mpo.O(d, sites(i)), gr);
                end
            end
        end
        
        function H = AC2_hamiltonian(mpo, mps, GL, GR, sites)
            arguments
                mpo
                mps
                GL = fixedpoint(transfermatrix(mpo, mps, 'Type', 'LL')) % BROKEN
                GR = fixedpoint(transfermatrix(mpo, mps, 'Type', 'RR').') % BROKEN
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
                    gl = twistdual(gl, 1);
                    gr = GR{d, next(sites(i), period(mps))};
                    gr = twistdual(gr, nspaces(gr));
                    H{i}(d, 1) = FiniteMpo(gl, {}, gr);
                end
            end
        end
        
        function H = B_hamiltonian(mpo, qp, GL, GR, sites, envopts, kwargs)
            arguments
                mpo
                qp
                GL
                GR
                sites = 1:period(mpo)
                envopts = {}
                kwargs.Type
            end
            
            switch kwargs.Type
                case {'l', 'left'}
                    GBL = leftquasienvironment(mpo, qp, GL, GR, envopts{:});
                    H = AC_hamiltonian(mpo, qp, GBL, GR, sites);
                case {'c', 'center'}
                    H = AC_hamiltonian(mpo, qp, GL, GR, sites);
                case {'r', 'right'}
                    GBR = rightquasienvironment(mpo, qp, GL, GR, envopts{:});
                    H = AC_hamiltonian(mpo, qp, GL, GBR, sites);
                otherwise
                    error('unknown type %s', kwargs.Type)
            end
        end
    end
    
    methods (Static)
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
