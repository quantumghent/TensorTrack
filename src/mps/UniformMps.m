classdef UniformMps
    % Implementation of infinite translation invariant MPS
    %
    % The center gauge is defined to have:
    %
    % .. math::
    %	AL_w \cdot C_w = AC_w = C_{w-1} \cdot AR_w
    %
    % Properties
    % ----------
    % AL : :class:`.MpsTensor`
    %   left-gauged mps tensors.
    %
    % AR : :class:`.MpsTensor`
    %   right-gauged mps tensors.
    %
    % C : :class:`.Tensor`
    %   center gauge transform.
    %
    % AC : :class:`.MpsTensor`
    %   center-gauged mps tensors.
    %
    % Todo
    % ----
    % Document all methods.
    
    properties
        AL (1,:) cell
        AR (1,:) cell
        C  (1,:) cell
        AC (1,:) cell
    end
    
    
    %% Constructors
    methods
        function mps = UniformMps(varargin)
            % Usage
            % -----
            % :code:`mps = UniformMps(A)`
            %
            % :code:`mps = UniformMps(AL, AR, C [, AC])`
            %
            % Arguments
            % ---------
            % A : :class:`cell` of :class:`.MpsTensor`
            %   set of tensors per site that define an MPS to be gauged.
            %
            % AL, AR, AC : :class:`cell` of :class:`.MpsTensor`
            %   set of gauged MpsTensors.
            %
            % C : :class:`cell` of :class:`.Tensor`
            %   gauge tensor.
            %
            % Returns
            % -------
            % mps : :class:`.UniformMps`
            %   gauged uniform MPS.
            
            narginchk(0, 4);
            
            if nargin == 0, return; end % default empty constructor
            
            if nargin == 1
                if isa(varargin{1}, 'UniformMps') % copy constructor
                    for i = numel(varargin{1}):-1:1
                        mps(i).AL = varargin{1}(i).AL;
                        mps(i).AR = varargin{1}(i).AR;
                        mps(i).C  = varargin{1}(i).C;
                        mps(i).AC = varargin{1}(i).AC;
                    end
                    mps = reshape(mps, size(mps));
                    
                elseif isa(varargin{1}, 'Tensor')
                    for i = height(varargin{1}):-1:1
                        mps.AR = MpsTensor(varargin{1});
                        mps.AL = MpsTensor(varargin{1});
                    end
                    mps = canonicalize(mps);
                    
                elseif isa(varargin{1}, 'MpsTensor')
                    % by default we take the input as AR, as canonicalize right
                    % orthogonalizes before left orthogonalization
                    for i = height(varargin{1}):-1:1
                        mps.AR = varargin{1};
                        mps.AL = varargin{1};
                    end
                    mps = canonicalize(mps);
                    
                elseif iscell(varargin{1})
                    for j = height(varargin{1}):-1:1
                        mps(j).AR = varargin{1}(j, :);
                        mps(j).AL = varargin{1}(j, :);
                    end
                    mps = canonicalize(mps);
                else
                    error('Invalid constructor for UniformMps.')
                end
                
            else
                if iscell(varargin{1}), mps.AL = varargin{1}; else, mps.AL = varargin(1); end
                if iscell(varargin{2}), mps.AR = varargin{2}; else, mps.AR = varargin(2); end
                if iscell(varargin{3}), mps.C = varargin{3};  else, mps.C = varargin(3); end
                if nargin == 4
                    if iscell(varargin{4})
                        mps.AC = varargin{4};
                    else
                        mps.AC = varargin(4);
                    end
                end
            end
            
            
        end
    end
    
    methods (Static)
        function mps = new(fun, pspaces, vspaces)
            % Create a uniform matrix product state with data using a function handle.
            %
            % Usage
            % -----
            % :code:`UniformMps.new(fun, pspaces, vspaces)`
            %
            % Arguments
            % ---------
            % fun : :class:`function_handle`
            %   function to initialize the tensor.
            %
            % Repeating Arguments
            % -------------------
            % pspaces : :class:`.AbstractSpace`
            %   physical spaces for each site.
            %
            % vspaces : :class:`.AbstractSpace`
            %   virtual spaces between each site. (entry `i` corresponds to left of site
            %   `i`.)
            
            arguments
                fun = []
            end
            arguments (Repeating)
                pspaces
                vspaces
            end
            
            if isempty(fun), fun = @randnc; end
            L = length(pspaces);
            
            for w = length(pspaces):-1:1
                rankdeficient = vspaces{w} * prod(pspaces{w}) < vspaces{next(w, L)} || ...
                        vspaces{w} > prod(pspaces{w})' * vspaces{next(w, L)};
                if rankdeficient
                    error('mps:rank', ...
                        'Cannot create a full rank mps with given spaces.');
                end
                
                A{w} = MpsTensor(Tensor.new(fun, [vspaces{w} pspaces{w}], ...
                    vspaces{next(w, L)}), 0);
            end
            
            mps = UniformMps(A);
        end
        
        function mps = randnc(pspaces, vspaces)
            % Create a uniform matrix product state with random entries.
            %
            % See Also
            % --------
            % :meth:`.UniformMps.new`
            
            arguments (Repeating)
                pspaces
                vspaces
            end
            args = [pspaces; vspaces];
            mps = UniformMps.new(@randnc, args{:});
        end
    end
    
    
    %% Properties
    methods
        function p = period(mps)
            % Period over which the mps is translation invariant.
            p = length(mps(1).AL);
        end
        
        function d = depth(mps)
            % Number of lines in a multi-line mps.
            d = size(mps, 1);
        end
        
        function mps = horzcat(mps, mps2, varargin)
            mps.AL = horzcat(mps.AL, mps2.AL);
            mps.AR = horzcat(mps.AR, mps2.AR);
            mps.C  = horzcat(mps.C, mps2.C);
            mps.AC = horzcat(mps.AC, mps2.AC);
            if nargin > 2
                mps = horzcat(mps, varargin{:});
            end
        end
        
        function mps = vertcat(varargin)
            for i = 2:length(varargin)
                assert(period(varargin{1}) == period(varargin{i}), ...
                    'Can only stack UniformMps with matching periods.')
            end
            mps = builtin('vertcat', varargin{:});
        end
        
        function s = leftvspace(mps, w)
            % Return the virtual space to the left of site `w`.
            if nargin == 1 || isempty(w), w = 1:period(mps); end
            s = arrayfun(@leftvspace, mps.AL{w});
        end
        
        function s = pspace(mps, w)
            % Return the physical space at site `w`.
            if nargin == 1 || isempty(w), w = 1:period(mps); end
            s = pspace(mps.AL{w});
        end
        
        function s = rightvspace(mps, w)
            % Return the virtual space to the right of site `w`.
            if nargin == 1 || isempty(w), w = 1:period(mps); end
            s = arrayfun(@rightvspace, mps.AL{w});
        end
        
        function type = underlyingType(mps)
            type = underlyingType(mps(1).AR{1});
        end
    end
    
    
    %% Methods
    methods
        function [mps, lambda] = canonicalize(mps, kwargs)
            % Compute the center-gauged form of an mps.
            %
            % Usage
            % -----
            % :code:`[mps, lambda] = canonicalize(mps, kwargs)`
            %
            % Arguments
            % ---------
            % mps : :class:`.UniformMps`
            %   input mps, from which AL or AR is used as the state, and optionally C as an
            %   initial guess for the gauge.
            %
            % Keyword Arguments
            % -----------------
            % Tol : :class:`numeric`
            %   tolerance for the algorithm.
            %
            % MaxIter : :class:`integer`
            %   maximum number of iterations.
            %
            % Method : :class:`char`
            %   algorithm used for decomposition. Must be 'polar', 'qr' or 'qrpos'.
            %
            % Verbosity : :class:`.Verbosity`
            %   level of output.
            %
            % DiagC : :class:`logical`
            %   flag to indicate if `C` needs to be diagonalized.
            %
            % ComputeAC : :class:`logical`
            %   flag to indicate if `AC` needs to be computed.
            %
            % Order : :class:`char`, 'lr' or 'rl'
            %   order of gauge fixing:
            %
            %   - 'lr' uses AL as input tensors, first :code:`leftorth`, then :class:`rightorth`.
            %   - 'rl' uses AR as input tensors, first :class:`rightorth`, then :code:`leftorth`.
            
            arguments
                mps
                kwargs.Tol = eps(underlyingType(mps))^(3/4)
                kwargs.MaxIter = 1e3
                kwargs.Method {mustBeMember(kwargs.Method, {'polar', 'qr', 'qrpos'})} ...
                    = 'polar'
                kwargs.Verbosity = Verbosity.warn
                kwargs.DiagC = false
                kwargs.ComputeAC = true
                kwargs.Order {mustBeMember(kwargs.Order, {'lr', 'rl'})} = 'lr'
            end
            
            for i = 1:depth(mps)
                if strcmp(kwargs.Order, 'rl')
                    [mps(i).AR, ~, ~, eta1]             = uniform_rightorth(...
                        mps(i), [], ...
                        'Tol', kwargs.Tol, 'MaxIter', kwargs.MaxIter, ...
                        'Method', kwargs.Method, 'Verbosity', kwargs.Verbosity);
                    mps(i).AL = mps(i).AR;
                    [mps(i).AL, mps(i).C, lambda, eta2] = uniform_leftorth(...
                        mps(i), mps(i).C, ...
                        'Tol', kwargs.Tol, 'MaxIter', kwargs.MaxIter, ...
                        'Method', kwargs.Method, 'Verbosity', kwargs.Verbosity);
                else
                    [mps(i).AL, ~, ~, eta1]             = uniform_leftorth(...
                        mps(i), {}, ...
                        'Tol', kwargs.Tol, 'MaxIter', kwargs.MaxIter, ...
                        'Method', kwargs.Method, 'Verbosity', kwargs.Verbosity);
                    mps(i).AR = mps(i).AL;
                    [mps(i).AR, mps(i).C, lambda, eta2] = uniform_rightorth(...
                        mps(i), mps(i).C, ...
                        'Tol', kwargs.Tol, 'MaxIter', kwargs.MaxIter, ...
                        'Method', kwargs.Method, 'Verbosity', kwargs.Verbosity);
                end
            end
            
            if kwargs.DiagC
                mps = diagonalizeC(mps);
            end
            
            if kwargs.ComputeAC
                for d = 1:depth(mps)
                    for w = 1:period(mps)
                        mps(d).AC{w} = computeAC(mps, d, w);
                    end
                end
            end
        end
        
        function AC = computeAC(mps, d, w, type)
            arguments
                mps
                d
                w
                type = 'L'
            end
            if strcmp(type, 'L')
                AC = multiplyright(mps(d).AL{w}, mps(d).C{w});
            elseif strcmp(type, 'R')
                AC = multiplyleft(mps(d).AR{w}, mps(d).C{prev(w, period(mps))});
            else
                error('mps:argerror', 'invalid type %s', type)
            end 
        end
        
        function AC2 = computeAC2(mps, d, w, type)
            arguments
                mps
                d
                w
                type = 'R'
            end
            
            ww = next(w, period(mps));

            if strcmp(type, 'R')
                assert(mps(d).AC{w}.alegs == 0, 'tba');
                AC2 = MpsTensor(contract(...
                        mps(d).AC{w}, ...
                        [-(1:mps(d).AC{w}.plegs+1), 1], ...
                        mps(d).AR{ww}, ...
                        [1, -(1:mps(d).AR{ww}.plegs+1) - 1 - mps(d).AR{ww}.plegs], ...
                        'Rank', [1 + mps(d).AC{w}.plegs, 1 + mps(d).AR{ww}.plegs]));
            elseif strcmp(type, 'L')
                assert(mps(d).AC{ww}.alegs == 0, 'tba');
                AC2 = MpsTensor(contract(...
                        mps(d).AL{w}, ...
                        [-(1:mps(d).AL{w}.plegs+1), 1], ...
                        mps(d).AC{ww}, ...
                        [1, -(1:mps(d).AR{ww}.plegs+1) - 1 - mps(d).AC{ww}.plegs], ...
                        'Rank', [1 + mps(d).AL{w}.plegs, 1 + mps(d).AR{ww}.plegs]));
            else
                error('mps:argerror', 'invalid type %s', type)
            end
        end
        
        function mps = diagonalizeC(mps)
            % Gauge transform an mps such that C is diagonal.
            
            for i = 1:depth(mps)
                for w = 1:period(mps(i))
                    C_iw = mps(i).C{w};
                    [U, S, V] = tsvd(C_iw, 1, 2);
                    
                    if isdual(C_iw.codomain) ~= isdual(S.codomain)
                        U.domain = conj(U.domain);
                        S.codomain = conj(S.codomain);
                        S = twist(S, 1);
                    end
                    
                    if isdual(C_iw.domain) ~= isdual(S.domain)
                        V.codomain = conj(V.codomain);
                        S.domain = conj(S.domain);
                        S = twist(S, 2);
                    end
                    
                    ww = next(w, period(mps(i)));
                    mps(i).C{w} = S;
                    
                    mps(i).AL{ww} = multiplyleft(mps(i).AL{ww}, U');
                    mps(i).AL{w}  = multiplyright(mps(i).AL{w}, U);
                    
                    mps(i).AR{ww} = multiplyleft(mps(i).AR{ww}, V);
                    mps(i).AR{w}  = multiplyright(mps(i).AR{w}, V');
                end
            end
        end
        
        function mps = normalize(mps)
            % Normalize an mps state.
            
            mps.C = arrayfun(@normalize, mps.C);
            mps.AC = arrayfun(@normalize, mps.AC);
        end
        
        function mps = circshift(mps, k)
            if length(k) > 1
                error('tba');
            end
            mps.AR = circshift(mps.AR, k);
            mps.AL = circshift(mps.AL, k);
            mps.C = circshift(mps.C, k);
            mps.AC = circshift(mps.AC, k);
        end
        
        function T = transfermatrix(mps1, mps2, sites, kwargs)
            % A finite matrix product operator that represents the transfer matrix of an
            % mps.
            %
            % Usage
            % -----
            % :code:`T = transfermatrix(mps1, mps2, sites, kwargs)`
            %
            % Arguments
            % ---------
            % mps1 : :class:`.UniformMps`
            %   input mps for top layer.
            %
            % mps2 : :class:`.UniformMps`
            %   input mps for bottom layer, by default equal to the top.
            %
            % sites : :class:`int`
            %   optionally slice the unit cell of the mps and only define the transfer
            %   matrix for this slice.
            %
            % Keyword Arguments
            % -----------------
            % Type : :class:`char`
            %   'LL', 'LR', 'RL', 'RR' to determine if the top or bottom respectively are AL
            %   or AR.
            %
            % Returns
            % -------
            % T : :class:`.FiniteMpo`
            %   transfer matrix of an mps, acting to the left.
            
            arguments
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
            
            T = cellfun(@transfermatrix, A1, A2);
        end
        
        function rho = fixedpoint(mps, type, w)
            % compute the fixed point of the transfer matrix of an mps.
            %
            % Usage
            % -----
            % :code:`rho = fixedpoint(mps, type, w)`
            %
            % Arguments
            % ---------
            % mps : :class:`.UniformMps`
            %   input state.
            %
            % type : :class:`char`
            %   specification of the type of transfer matrix: 
            %   general format: :code:`sprintf(%c_%c%c, side, top, bot)` where side is 'l' or 'r' to
            %   determine which fixedpoint, and top and bot are 'L' or 'R' to specify
            %   whether to use AL or AR in the transfer matrix.
            %
            % w : :class:`int`
            %   position within the mps unitcell of the fixed point.
            
            arguments
                mps
                type {mustBeMember(type, ...
                    {'l_LL' 'l_LR' 'l_RL' 'l_RR' 'r_LL' 'r_LR' 'r_RL' 'r_RR'})}
                w = strcmp(type(1), 'l') * 1 + strcmp(type(1), 'r') * period(mps)
            end
            
            ww = prev(w, period(mps));
            switch type
                case 'l_RR'
                    rho = contract(mps.C{ww}', [-1 1], mps.C{ww}, [1 -2], 'Rank', [1 1]);
%                     if isdual(space(rho, 1)), rho = twist(rho, 1); end
                case 'l_RL'
                    rho = mps.C{ww};
%                     if isdual(space(rho, 1)), rho = twist(rho, 1); end
                case 'l_LR'
                    rho = mps.C{ww}';
%                     if isdual(space(rho, 1)), rho = twist(rho, 1); end
                case 'l_LL'
                    rho = mps.C{w}.eye(leftvspace(mps, w), leftvspace(mps, w));
                    if isdual(space(rho, 1)), rho = twist(rho, 1); end
                    
                case 'r_RR'
                    rho = mps.C{w}.eye(rightvspace(mps, w)', rightvspace(mps, w)');
                    if isdual(space(rho, 2)), rho = twist(rho, 2); end
                case 'r_RL'
                    rho = twist(mps.C{w}', 2);
%                     if isdual(space(rho, 1)), rho = twist(rho, 2); end
                case 'r_LR'
                    rho = twist(mps.C{w}, 2);
%                     if ~isdual(space(rho, 2)), rho = twist(rho, 2); end
                case 'r_LL'
                    rho = contract(mps.C{w}, [-1 1], mps.C{w}', [1 -2], 'Rank', [1 1]);
                    rho = twist(rho, 2);
            end
        end
        
        function [V, D] = transfereigs(mps1, mps2, howmany, which, eigopts, kwargs)
            % Compute the eigenvalues of the transfer matrix of an mps.
            %
            % Usage
            % -----
            % :code:`[V, D] = transfereigs(mps1, mps2, howmany, which, eigopts, kwargs)`
            %
            % Arguments
            % ---------
            % mps1 : :class:`.UniformMps`
            %   input mps for top layer.
            %
            % mps2 : :class:`.UniformMps`
            %   input mps for bottom layer. Default value equal to `mps1`.
            %
            % howmany : :class:`int`
            %   number of eigenvectors and eigenvalues to compute.
            %
            % which : :class:`char`
            %   type of eigenvectors to target.
            %
            % Keyword Arguments
            % -----------------
            % eigopts
            %   see keyword arguments for :func:`.eigsolve`.
            %
            % Verbosity : :class:`int`
            %   detail level for output.
            %
            % Type : :class:`char`
            %   type of transfer matrix to construct, see
            %   :meth:`.UniformMps.transfermatrix`.
            %
            % Charge : :class:`.AbstractCharge`
            %   charge of eigenvectors to target.
            
            arguments
                mps1
                mps2 = mps1
                howmany = min(20, dim(leftvspace(mps1, 1) * leftvspace(mps2, 1)'))
                which = 'largestabs'
                eigopts.Algorithm = 'Arnoldi'
                eigopts.KrylovDim = 100
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
                kwargs.Verbosity = 0
                kwargs.Type {mustBeMember(kwargs.Type, ...
                    {'l_LL' 'l_LR' 'l_RL' 'l_RR' 'r_LL' 'r_LR' 'r_RL' 'r_RR'})} = 'r_RR'
                kwargs.Charge = []
            end
            
            T = transfermatrix(mps1, mps2, 'Type', kwargs.Type(3:4));
            if kwargs.Type(1) == 'r', T = T'; end
            
            if ~isempty(kwargs.Charge)
                Tdomain = domain(T);
                auxspace = Tdomain.new(...
                    struct('charges', kwargs.Charge, 'degeneracies', 1), ...
                    false);
                v0 = Tensor.randnc([Tdomain auxspace], []); 
            else
                v0 = [];
            end
            
            eigopts.Verbosity = kwargs.Verbosity;
            eigkwargs = namedargs2cell(eigopts);
            [V, D] = eigsolve(T, v0, howmany, which, eigkwargs{:});

            if kwargs.Type(1) == 'r', V = V'; end
            if nargout < 2, V = diag(D); end
        end
        
        function f = fidelity(mps1, mps2, kwargs)
            % Compute the fidelity between two uniform MPSs.
            
            arguments
                mps1
                mps2
                kwargs.KrylovDim = 30
                kwargs.MaxIter = 500
                kwargs.ReOrth = 2
                kwargs.Tol = eps(underlyingType(mps1))^(3/4)
                kwargs.Verbosity = 0
            end
            
            eigkwargs = [fieldnames(kwargs).'; struct2cell(kwargs).'];
            f = transfereigs(mps1, mps2, 1, 'largestabs', eigkwargs{:});
        end
        
        function E = expectation_value(mps1, O, mps2)
            % Compute the expectation value of an operator.
            
            arguments
                mps1
                O
                mps2 = mps1
            end
            
            if isa(O, 'InfJMpo')
                [GL, GR] = environments(O, mps1, mps2);
                Hs = AC_hamiltonian(O, mps1, GL, GR, period(mps1));
                H = Hs{1};
                N = size(H.R.var, 2);
                H.R.var = H.R.var(1, N, 1);
                H.O{1} = H.O{1}(:, :, N, :);
                AC_ = apply(H, mps1.AC{end});
                E = dot(AC_, mps2.AC{end});
                
            elseif isa(O, 'InfMpo')
                [GL, GR] = environments(O, mps1, mps2); % should be normalized
                Hs = AC_hamiltonian(O, mps1, GL, GR);
                E = zeros(size(Hs));
                for i = 1:length(Hs)
                    AC_ = apply(Hs{i}, mps1.AC{i});
                    E(i) = dot(AC_, mps2.AC{i});
                end
                E = prod(E);
                
            elseif isa(O, 'AbstractTensor')
                E = local_expectation_value(mps1, O);
            else
                error('Unknown operator type (%s)', class(O));
            end
        end
        
        function E = local_expectation_value(mps, O, offset)
            % Compute the expectation value of a local operator.
            
            arguments
                mps
                O
                offset = 0   % site offset
            end
            
            local_ops = MpoTensor.decompose_local_operator(O);
            N = length(local_ops);
            
            A = [mps.AC(1 + offset) mps.AR(2:(N - 1) + offset)];
            
            T = FiniteMpo.mps_channel_operator(A, local_ops, A);
            rhoL = insert_onespace(fixedpoint(mps, 'l_LL', offset+ 1), 2, ...
                ~isdual(space(T(1).O{1}, 4)));
            rhoR = insert_onespace(fixedpoint(mps, 'r_RR', offset + N), 2, ...
                ~isdual(space(T(1).O{1}, 2)));
            E1 = apply(T, rhoL);
            E = contract(E1, 1:3, rhoR, flip(1:3));
        end
        
        function [svals, charges] = schmidt_values(mps, w)
            % Compute the Schmidt values and corresponding charges for an entanglement cut
            % to the right of site :code:`w`.
            %
            % Usage
            % -----
            % :code:`[svals, charges] = schmidt_values(mps, w)`
            
            arguments
                mps
                w = 1
            end
            
            [svals, charges] = matrixblocks(tsvd(mps.C{w}));
            svals = cellfun(@diag, svals, 'UniformOutput', false);
        end
        
        function plot_entanglementspectrum(mps, d, w, ax, kwargs)
            % Plot the entanglement spectrum of a uniform MPS.
            
            arguments
                mps
                d = 1:depth(mps)
                w = 1:period(mps)
                ax = []
                kwargs.SymmetrySort = true
                kwargs.ExpandQdim = false
            end
            if isempty(ax)
                figure;
                ax = gobjects(depth(mps), width(mps));
                for dd = 1:length(d)
                    for ww = 1:length(w)
                        ax(dd, ww) = subplot(length(d), length(w), ww + (dd-1)*length(w));
                    end
                end
            end
            for dd = 1:length(d)
                for ww = 1:length(w)
                    hold(ax(dd,ww), 'off');
                    [svals, charges] = schmidt_values(mps(dd), w(ww));
                    if kwargs.ExpandQdim
                        for i = 1:length(svals)
                            svals{i} = reshape(repmat(svals{i}, 1, qdim(charges(i))), [], 1);
                        end
                    end
                    ctr = 0;
                    labels = arrayfun(@string, charges, 'UniformOutput', false);
                    lengths = cellfun(@length, svals);
                    ticks = cumsum(lengths);
                    if kwargs.SymmetrySort
                        for i = 1:length(svals)
                            semilogy(ax(dd, ww), ctr+(1:lengths(i)), svals{i}, '.', 'MarkerSize', 10, 'Color', colors(i));
                            if i == 1, hold(ax(dd,ww), 'on'); end
                            ctr = ctr + lengths(i);
                        end
                        set(ax(dd, ww), 'Xtick', ticks, 'fontsize', 10, ...
                            'XtickLabelRotation', 60, 'Xgrid', 'on');
                    else
                        [~, p] = sort(vertcat(svals{:}), 'descend');
                        p = invperm(p);
                        for i = 1:length(svals)
                            semilogy(ax(dd, ww), p(ctr+(1:lengths(i))), svals{i}, '.', 'MarkerSize', 10, 'Color', colors(i));
                            if i == 1, hold(ax(dd,ww), 'on'); end
                            ctr = ctr + lengths(i);
                        end
                        
                    end
                    legend(ax(dd, ww), labels)
                    set(ax(dd, ww), 'TickLabelInterpreter', 'latex');
                    xlim(ax(dd, ww), [1 - 1e-8 ticks(end) + 1e-8]);
                end
            end
            hold off
            linkaxes(ax, 'y');
        end
        
        function mps = desymmetrize(mps)
            if numel(mps) > 1
                mps = arrayfun(@desymmetrize, mps);
                return
            end
            mps.AL = desymmetrize(mps.AL);
            mps.AR = desymmetrize(mps.AR);
            mps.C = desymmetrize(mps.C);
            mps.AC = desymmetrize(mps.AC);
        end
        
        function [xi, theta] = correlation_length(mps, charge)
            % Compute the correlation length of an MPS in a given charge sector.
            %
            % Usage
            % -----
            % :code:`[xi, theta] = correlation_length(mps, charge)`
            %
            % Arguments
            % ---------
            % mps : :class:`.UniformMps`
            %   input mps.
            %
            % charge : :class:`.AbstractCharge`
            %   charge sector for correlation length to target.
            %
            % Returns
            % -------
            % xi : :class:`numeric`
            %   correlation length in the given charge sector.
            %
            % theta : :class:`numeric`
            %   angle of the corresponding oscillation period.
            
            arguments
                mps
                charge = []
            end
            
            if isempty(charge) || charge == one(charge)
                f = transfereigs(mps, mps, 5, 'KrylovDim', 30);
                if abs(f(1) - 1) > 1e-12
                    warning('mps:noninjective', ...
                        'mps might be non-injective:\n%s', num2str(f));
                end
                epsilon = -log(abs(f(2)));
                theta = angle(f(2));
            else
                f = transfereigs(mps, mps, 1, 'KrylovDim', 20, 'Charge', charge);
                epsilon = -log(abs(f));
                theta = angle(f);
            end
            
            xi = 1 / epsilon;
        end
        
        function [epsilon, delta, spectrum] = marek_gap(mps, charge, kwargs)
            % Compute the Marek gap of an MPS in a given charge sector.
            %
            % Usage
            % -----
            % :code:`[epsilon, delta, spectrum] = marek_gap(mps, charge, kwargs)`
            %
            % Arguments
            % ---------
            % mps : :class:`.UniformMps`
            %   input mps.
            %
            % charge : :class:`.AbstractCharge`
            %   charge sector for correlation length to target.
            %
            % Keyword Arguments
            % -----------------
            % HowMany : :class:`int`
            %   amount of transfer matrix eigenvalues to compute.
            %
            % Angle : :class:`numeric`
            %   angle in radians around which the gap should be computed.
            %
            % AngleTol : :class:`numeric`
            %   tolerance in radians for angles to be considered equal.
            %
            % Returns
            % -------
            % epsilon : :class:`numeric`
            %   inverse correlation length in the given charge sector.
            %
            % delta : :class:`numeric`
            %   refinement parameter.
            %
            % spectrum : :class:`numeric`
            %   computed partial transfer matrix spectrum.
            arguments
                mps
                charge = []
                kwargs.Angle
                kwargs.AngleTol = 1e-1
                kwargs.HowMany = 5
            end
            
            spectrum = transfereigs(mps, mps, kwargs.HowMany, 'largestabs', 'Charge', charge);
            [d, p] = sort(abs(spectrum), 'descend');
            inds = d > 1 - 1e-12;
            
            if ((isempty(charge) || charge == one(charge)) && sum(inds) > 1) || ...
                    sum(inds) > 0
                warning('mps:noninjective', ...
                    'mps might be non-injective:\n%s', num2str(spectrum));
            end
            
            d(inds) = [];
            p(inds) = [];
            
            if abs(diff(unwrap(angle(spectrum(p(1:2)))))) > 1e-1
                warning('comparing values with different angles:\n%s', num2str(spectrum));
            end
            if isfield('Angle', kwargs)
                error('tba');
            end
            
            epsilon = -log(d(1));
            delta = log(d(1) / d(2));
        end

        function S = entanglement_entropy(mps, w)
            % Compute the entanglement entropy of a uniform MPS for a cut to the right of
            % site :code:`w`
            arguments
                mps
                w = 1
            end
            
            [svals, charges] = schmidt_values(mps, w);
            
            S = 0;
            for i = 1:length(svals)
                S = S - qdim(charges(i)) * sum(svals{i}.^2 .* log(svals{i}.^2));
            end
        end
        
        function S = renyi_entropy(mps, n, w)
            % Compute the :code:`n`-th Renyi entropy of a uniform MPS for a cut to the right
            % of site :code:`w`
           arguments
                mps
                n
                w = 1
            end
            
            [svals, charges] = schmidt_values(mps, w);
            
            S = 0;
            for i = 1:length(svals)
                S = S + qdim(charges(i)) * sum(svals{i}.^(2 * n));
            end
            S = 1 / (1 - n) * log(S);
        end

        function out = truncate(mps, trunc)
            % Truncate a uniform MPS according to the options specified in :code:`trunc`
            % (see :meth:`.Tensor.tsvd` for details on the truncation options).
            arguments
                mps
                trunc.TruncDim
                trunc.TruncTotalDim
                trunc.TruncBelow
                trunc.TruncSpace
            end

            trunc = [fieldnames(trunc),struct2cell(trunc)]';       
            out = repmat(UniformMps, size(mps, 1), 1);
            for d = 1:depth(mps)
                for w = 1:period(mps(d))
                    [~, ~, V] = tsvd(mps(d).C(w), 1, 2, trunc{:});
                    ww = next(w, period(mps(d)));
                    mps(d).AR(ww) = multiplyleft(mps(d).AR(ww), V);
                    mps(d).AR(w)  = multiplyright(mps(d).AR(w), V');
                end
                % bring truncated mps to canonical form
                out(d) = UniformMps(mps(d).AR);
            end
        end
        
%         S = EntanglementEntropy(mps, loc);
%         S = RenyiEntropy(mps,n, loc);
%         E = ExpectationValue(mps, W, GL, GR)
%         rho = LeftFixedPoint(mps1, mps2, w, choice)
%         rho = RightFixedPoint(mps1, mps2, w, choice)
%         sf=StaticStructureFactor(mps,S,k)
%         
%         out = Block(mps, opts)
%         out = Split(mps, varargin)
%         %[out, lambda] = Truncate(mps, control, opts)
%         
%         [f, rho] = Fidelity(mps1, mps2, tol)
%         
%         mps = Conj(mps)
%         
%         mps = mtimes(mps, lambda)
%         
%         mps = ShiftUnitCell(mps,dd,dw)
%         out = Rotate180(mps)
%         out = Transpose(mps)
%         
%         mps = SendToGpu(mps)
%         mps = GetFromGpu(mps)
%         
%         [mps, xi] = Retract(mps, eta, alpha)
%         n = Inner(x, eta, xi)
        
        
    end
    
    %% Subroutines
    methods (Access = protected)
        function [AL, CL, lambda, eta] = uniform_leftorth(mps, CL, kwargs)
            % Bring a uniform MPS into left canonical form.
            %
            % Usage
            % -----
            % :code:`[AL, CL, lambda, eta] = uniform_leftorth(mps, CL, kwargs)`
            arguments
                mps
                CL = {}
                kwargs.Tol = eps(underlyingType(mps))^(3/4)
                kwargs.MaxIter = 500
                kwargs.Method = 'polar'
                kwargs.Verbosity = Verbosity.off
                kwargs.Normalize = true
                kwargs.EigsInit = 3
                kwargs.EigsFrequence = 2
            end
            
            A = mps.AL;
%             
            % constants
            EIG_TOLFACTOR = 1/50;
            EIG_MAXTOL = 1e-4;
            MINKRYLOVDIM = 8;
            MAXKRYLOVDIM = 30;
            
            
            
            % initialization
            N = period(mps);
            if isempty(CL)
                CL = initializeC(A{:});
            end
            
            if kwargs.Normalize, CL{1} = normalize(CL{1}); end
            
            for i = 1:numel(A)
                A{i} = repartition(A{i}, [nspaces(A{i}) - 1, 1]);
            end
            AL = A;
            
            eta_best = Inf;
            ctr_best = 0;
            AL_best = AL;
            C_best = CL;
            lambda_best = 0;
            
            for ctr = 1:kwargs.MaxIter
                if ctr > kwargs.EigsInit && mod(ctr, kwargs.EigsFrequence) == 0
                    C_ = repartition(CL{end}, [2 0]);
                    T = cellfun(@transfermatrix, A, AL);
                    
                    eigalg = Arnoldi('Tol', min(eta * EIG_TOLFACTOR, EIG_MAXTOL), ...
                        'KrylovDim', between(MINKRYLOVDIM, ...
                            MAXKRYLOVDIM - ctr / 2 + 4, ...
                            MAXKRYLOVDIM), ...
                        'NoBuild', 4, ...
                        'Verbosity', kwargs.Verbosity - 1);
                    
                    [C_, ~] = eigsolve(eigalg, @(x) T.apply(x), C_, 1, 'largestabs');
                    [~, CL{end}] = leftorth(C_, 1, 2, kwargs.Method);
                    if isdual(space(CL{end}, 2)) == isdual(space(CL{end}, 1))
                        CL{end}.codomain = conj(CL{end}.codomain);
                        CL{end} = twist(CL{end}, 1);
                    end
                end
                
                C_ = CL{end};
                lambdas = ones(1, N);
                for w = 1:N
                    ww = prev(w, N);
                    CA = multiplyleft(A{w}, CL{ww});
                    [AL{w}, CL{w}] = leftorth(CA, kwargs.Method);
                    lambdas(w) = norm(CL{w});
                    if kwargs.Normalize, CL{w} = CL{w} ./ lambdas(w); end
                end
                lambda = prod(lambdas);
                eta = norm(C_ - CL{end}, Inf);
                if eta < kwargs.Tol
                    if kwargs.Verbosity >= Verbosity.conv
                        fprintf('Conv %2d:\terror = %0.4e\n', ctr, eta);
                    end
                    break;
                elseif eta < eta_best
                    eta_best = eta;
                    ctr_best = ctr;
                    AL_best = AL;
                    C_best = CL;
                    lambda_best = lambda;
                elseif ctr > 40 && ctr - ctr_best > 5
                    warning('uniform_orthright:stagnate', 'Algorithm stagnated');
                    eta = eta_best;
                    AL = AL_best;
                    CL = C_best;
                    lambda = lambda_best;
                    break;
                end
                
                if kwargs.Verbosity >= Verbosity.iter
                    fprintf('Iter %2d:\terror = %0.4e\n', ctr, eta);
                end
            end
            
            if kwargs.Verbosity >= Verbosity.warn && eta > kwargs.Tol
                fprintf('Not converged %2d:\terror = %0.4e\n', ctr, eta);
            end
        end
        
        function [AR, CR, lambda, eta] = uniform_rightorth(mps, CR, kwargs)
            % Bring a uniform MPS into right canonical form.
            %
            % Usage
            % -----
            % :code:`[AR, CR, lambda, eta] = uniform_rightorth(mps, CR, kwargs)`
            arguments
                mps
                CR = mps.C
                kwargs.Tol = eps(underlyingType(mps))^(3/4)
                kwargs.MaxIter = 500
                kwargs.Method = 'polar'
                kwargs.Verbosity = 0
                kwargs.Normalize = true
                kwargs.EigsInit = 3
                kwargs.EigsFrequence = 2
            end
            
            opts = namedargs2cell(kwargs);
            mps_d = mps;
            mps_d.AL = flip(cellfun(@ctranspose, mps.AR, 'UniformOutput', false));
            
            if isempty(CR)
                Cd = {};
            else
                Cd = circshift(flip(cellfun(@ctranspose, CR, 'UniformOutput', false)), -1);
            end
            
            [ARd, CRd, lambda, eta] = uniform_leftorth(mps_d, Cd, opts{:});
            
            AR = flip(cellfun(@ctranspose, ARd, 'UniformOutput', false));
            CR  = circshift(flip(cellfun(@ctranspose, CRd, 'UniformOutput', false)), -1);
            lambda = conj(lambda);
        end
    end
end

