classdef UniformMps
    %UNIFORMMPS Implementation of infinite translation invariant MPS
    %   MPS is stored in center gauge, where
    %       AL(w) * C(w) = AC(w) = C(w-1) * AR(w)
    
    properties
        AL (1,:) MpsTensor  = MpsTensor.empty(1, 0)
        AR (1,:) MpsTensor  = MpsTensor.empty(1, 0)
        C  (1,:) Tensor     = Tensor.empty(1, 0)
        AC (1,:) MpsTensor  = MpsTensor.empty(1, 0)
    end
    
    
    %% Constructors
    methods
        function mps = UniformMps(varargin)
            
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
                    % by default we take the input as AR, as canonicalize right
                    % orthogonalizes before left orthogonalization
                    for i = height(varargin{1}):-1:1
                        mps(i).AR = varargin{1}(i, :);
                        mps(i).AL = varargin{1}(i, :);
                    end
                    mps = canonicalize(mps);
                    
                else
                    error('Invalid constructor for UniformMps.')
                end
                
            elseif nargin == 4
                mps.AL = varargin{1};
                mps.AR = varargin{2};
                mps.C = varargin{3};
                mps.AC = varargin{4};
                
            else
                error('Invalid constructor for UniformMps.')
            end
        end
    end
    
    methods (Static)
        function mps = new(fun, pspaces, vspaces)
            if isempty(fun), fun = @randnc; end
            
            if isscalar(pspaces) && ~isscalar(vspaces)
                pspaces = repmat(pspaces, size(vspaces));
            elseif isscalar(vspaces) && ~isscalar(pspaces)
                vspaces = repmat(vspaces, size(pspaces));
            else
                assert(isequal(size(vspaces), size(pspaces)), 'mps:dimagree', ...
                    'Invalid sizes of input spaces.');
            end
            
            for w = length(pspaces):-1:1
                A(w) = Tensor.new(fun, [vspaces(w) pspaces(w)], ...
                    vspaces(next(w, length(pspaces))));
            end
            
            mps = UniformMps(A);
        end
        
        function mps = randnc(pspaces, vspaces)
            mps = UniformMps.new(@randnc, pspaces, vspaces);
        end
    end
    
    
    %% Properties
    methods
        function p = period(mps)
            for i = numel(mps):-1:1
                p(i) = length(mps(i).AL);
            end
        end
        
        function d = depth(mps)
            d = size(mps, 1);
        end
        
        function s = leftvspace(mps, w)
            if nargin == 1 || isempty(w), w = 1:period(mps); end
            s = arrayfun(@leftvspace, mps.AL(w));
        end
        
        function s = rightvspace(mps, w)
            if nargin == 1 || isempty(w), w = 1:period(mps); end
            s = arrayfun(@rightvspace, mps.AL(w));
        end
        
        function type = underlyingType(mps)
            type = underlyingType(mps(1).AR(1));
        end
    end
    
    
    %% Methods
    methods
        function [mps, lambda] = canonicalize(mps, kwargs)
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
            
            for i = 1:length(mps)
                if strcmp(kwargs.Order, 'lr')
                    [mps(i).AR, ~, ~, eta1]             = uniform_rightorth(...
                        mps(i).AR, [], ...
                        'Tol', kwargs.Tol, 'MaxIter', kwargs.MaxIter, ...
                        'Method', kwargs.Method, 'Verbosity', kwargs.Verbosity);
                    [mps(i).AL, mps(i).C, lambda, eta2] = uniform_leftorth(...
                        mps(i).AR, mps(i).C, ...
                        'Tol', kwargs.Tol, 'MaxIter', kwargs.MaxIter, ...
                        'Method', kwargs.Method, 'Verbosity', kwargs.Verbosity);
                else
                    [mps(i).AL, ~, ~, eta1]             = uniform_leftorth(...
                        mps(i).AL, [], ...
                        'Tol', kwargs.Tol, 'MaxIter', kwargs.MaxIter, ...
                        'Method', kwargs.Method, 'Verbosity', kwargs.Verbosity);
                    [mps(i).AR, mps(i).C, lambda, eta2] = uniform_rightorth(...
                        mps(i).AL, mps(i).C, ...
                        'Tol', kwargs.Tol, 'MaxIter', kwargs.MaxIter, ...
                        'Method', kwargs.Method, 'Verbosity', kwargs.Verbosity);
                end
            end
            
            if kwargs.DiagC
                mps = diagonalizeC(mps);
            end
            if kwargs.ComputeAC
                for i = 1:height(mps)
                    for w = period(mps(i)):-1:1
                        mps(i).AC(w) = multiplyleft(mps(i).AR(w), ...
                            mps(i).C(prev(w, period(mps(i)))));
                    end
                end
            end
        end
        
        function mps = diagonalizeC(mps)
            for i = 1:height(mps)
                for w = 1:period(mps(i))
                    C_iw = mps(i).C(w);
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
                    mps(i).C(w) = S;
                    
                    mps(i).AL(ww) = multiplyleft(mps(i).AL(ww), U');
                    mps(i).AL(w)  = multiplyright(mps(i).AL(w), U);
                    
                    mps(i).AR(ww) = multiplyleft(mps(i).AR(ww), V);
                    mps(i).AR(w)  = multiplyright(mps(i).AR(w), V');
                end
            end
        end
        
        function [V, D] = transfereigs(mps1, mps2, howmany, which, eigopts, kwargs)
            arguments
                mps1
                mps2 = mps1
                howmany = min(20, dim(leftvspace(mps1, 1) * leftvspace(mps2, 1)'))
                which = 'largestabs'
                eigopts.KrylovDim = 100
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
                kwargs.Verbosity = 0
                kwargs.Type {mustBeMember(kwargs.Type, ...
                    {'l_LL' 'l_LR' 'l_RL' 'l_RR' 'r_LL' 'r_LR' 'r_RL' 'r_RR'})} = 'l_LL'
                kwargs.Charge = []
            end
            
            % extract tensors
            if strcmp(kwargs.Type(3), 'L'), T = mps1.AL; else, T = mps1.AR; end
            if strcmp(kwargs.Type(4), 'L'), B = mps2.AL; else, B = mps2.AR; end
            
            % generate initial guess
            if strcmp(kwargs.Type(1), 'l')
                vspace1 = leftvspace(mps2, 1);
                vspace2 = leftvspace(mps1, 1);
            else
                vspace1 = rightvspace(mps1, period(mps));
                vspace2 = rightvspace(mps2, period(mps));
            end
            if isempty(kwargs.Charge) || isequal(kwargs.Charge, one(kwargs.Charge))
                v0 = Tensor.randnc([vspace1 vspace2'], []);
            else
                dims = struct('charges', kwargs.Charge, ...
                    'degeneracies', ones(size(kwargs.Charge)));
                auxspace = vspace1.new(dims, false);
                v0 = Tensor.randnc([vspace1 vspace2' auxspace], []);
            end
            
            % find eigenvectors
            eigkwargs = [fieldnames(eigopts).'; struct2cell(eigopts).'];
            if strcmp(kwargs.Type(1), 'l')
                [V, D] = eigsolve(@(x) applyleft(T, conj(B), x), v0, howmany, which, ...
                    eigkwargs{:}, 'Verbosity', kwargs.Verbosity);
            else
                [V, D] = eigsolve(@(x) applyright(T, conj(B), x), v0, howmany, which, ...
                    eigkwargs{:}, 'Verbosity', kwargs.Verbosity);
            end
            
            if nargout < 2, V = D; end
        end
        
        function f = fidelity(mps1, mps2, kwargs)
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
        
        function rho = fixedpoint(mps, type, w)
            arguments
                mps
                type {mustBeMember(type, ...
                    {'l_LL' 'l_LR' 'l_RL' 'l_RR' 'r_LL' 'r_LR' 'r_RL' 'r_RR'})}
                w = strcmp(type(1), 'l') * 1 + strcmp(type(1), 'r') * period(mps)
            end
            
            ww = prev(w, period(mps));
            switch type
                case 'l_RR'
                    rho = mps.C(ww)' * mps.C(ww);
                case 'l_RL'
                    rho = mps.C(ww);
                case 'l_LR'
                    rho = mps.C(ww)';
                case 'l_LL'
                    rho = mps.C.eye(leftvspace(mps, w), leftvspace(mps, w));
                case 'r_RR'
                    rho = mps.C.eye(rightvspace(mps, w)', rightvspace(mps, w)');
                case 'r_RL'
                    rho = mps.C(w)';
                case 'r_LR'
                    rho = mps.C(w);
                case 'r_LL'
                    rho = mps.C(w) * mps.C(w)';
            end
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
        
        % essential
        [mps, lambda] = Canonical(mps, options);
        [f, rho] = TransferEigs(mps1, mps2, x0, num, charge, choice, options);
        
        % utility
        PlotTransferEigs(mps1, mps2, x0, num, charges);
        [schmidt, charges] = SchmidtValues(mps, loc)
        PlotEntSpectrum(mps, svalue, opts)
        xi = CorrelationLength(mps, charge);
        [epsilon, delta, spectrum] = MarekGap(mps, charge, angle, num)
        S = EntanglementEntropy(mps, loc);
        S = RenyiEntropy(mps,n, loc);
        E = ExpectationValue(mps, W, GL, GR)
        rho = LeftFixedPoint(mps1, mps2, w, choice)
        rho = RightFixedPoint(mps1, mps2, w, choice)
        sf=StaticStructureFactor(mps,S,k)
        
        out = Block(mps, opts)
        out = Split(mps, varargin)
        [out, lambda] = Truncate(mps, control, opts)
        
        [f, rho] = Fidelity(mps1, mps2, tol)
        
        mps = Conj(mps)
        
        mps = mtimes(mps, lambda)
        
        mps = ShiftUnitCell(mps,dd,dw)
        out = Rotate180(mps)
        out = Transpose(mps)
        
        mps = SendToGpu(mps)
        mps = GetFromGpu(mps)
        
        [mps, xi] = Retract(mps, eta, alpha)
        n = Inner(x, eta, xi)
        
        mps = TensorNone(mps)
        
        % to be done later
        
        Add; % is this useful?
        MultiplyLeft; % no clue
        MultiplyRight; % no clue
        
    end
    
    methods (Static)
        
        [AR, C, lambda, error] = OrthRight(A, C, options)
        [AL, C, lambda, error] = OrthLeft (A, C, options)
        out = Random(virtualLeg, varargin);
        
        % conversion from legacy datastructure
        mps = FromCell(A);
        
    end
    
    
    methods (Access = protected)
        
        % anything?
        
    end
    
end

