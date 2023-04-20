classdef (InferiorClasses = {?Tensor, ?SparseTensor}) MpsTensor < AbstractTensor
    % Generic mps tensor objects that have a notion of virtual, physical and auxiliary legs.
    
    properties
        var
        plegs
        alegs = 0
    end
    
    
    %% Constructors
    methods
        function A = MpsTensor(tensor, alegs)
            arguments
                tensor = []
                alegs = 0
            end
            
            if iscell(tensor)
                for i = length(tensor):-1:1
                    A(i) = MpsTensor(tensor{i}, alegs);
                end
                return
            end
            
            if ~isempty(tensor)
                A.var = tensor;
                A.plegs = nspaces(tensor) - alegs - 2;
                A.alegs = alegs;
            end
        end
    end
    
    
    %% Properties
    methods
        function s = space(A, varargin)
            s = space(A.var, varargin{:});
        end
        
        function n = nspaces(A)
            n = nspaces(A.var);
        end
        
        function s = pspace(A)
            s = space(A, 1 + (1:A.plegs));
        end
        
        function s = leftvspace(A)
            s = space(A.var(1), 1);
        end
        
        function s = rightvspace(A)
            s = space(A.var(1), nspaces(A.var(1)) - A.alegs);
        end
        
        function cod = codomain(A)
            cod = A.var.codomain;
        end
        
        function dom = domain(A)
            dom = A.var.domain;
        end
        
        function r = rank(A)
            r = rank([A.var]);
        end
        
        function tdst = insert_onespace(tsrc, varargin)
            % insert a trivial space at position i.
            tdst = MpsTensor(insert_onespace(tsrc.var, varargin{:}), tsrc.alegs + 1);
        end
    end
    
    
    %% Linear Algebra
    methods
        function [A, L] = leftorth(A, alg)
            arguments
                A
                alg = 'polar'
            end
            
            if numel(A) > 1
                for i = numel(A):-1:1
                    [A(i), L(i)] = leftorth(A(i), alg);
                end
                return
            end
            
            if A.alegs == 0
                [A.var, L] = leftorth(A.var, 1:nspaces(A)-1, nspaces(A), alg);
                if isdual(space(L, 1)) == isdual(space(L, 2))
                    L.codomain = conj(L.codomain);
                    L = twist(L, 1);
                    A.var.domain = conj(A.var.domain);
                end
            else
                [A.var, L] = leftorth(A.var, [1:A.plegs+1 A.plegs+3], A.plegs+2, alg);
                A.var = permute(A.var, [1:A.plegs+1 A.plegs+3 A.plegs+2], rank(A));
            end
        end
        
        function [R, A] = rightorth(A, alg)
            arguments
                A
                alg = 'rqpos'
            end
            
            if numel(A) > 1
                for i = numel(A):-1:1
                    [R(i), A(i)] = rightorth(A, alg);
                end
                return
            end
            
            [R, A.var] = rightorth(A.var, 1, 2:nspaces(A), alg);
            if isdual(space(R, 1)) == isdual(space(R, 2))
                R.domain = conj(R.domain);
                R = twist(R, 2);
                A.var.codomain = conj(A.var.codomain);
            end
        end
        
        function A = leftnull(A, alg)
            arguments
                A
                alg = 'svd'
            end
            
            if numel(A) > 1
                for i = numel(A):-1:1
                    A(i) = leftnull(A(i), alg);
                end
                return
            end
            
            p = 1:nspaces(A);
            p1 = p(1:(end - A.alegs - 1));
            p2 = p((end - A.alegs):end);
            A.var = leftnull(A.var, p1, p2, alg);
        end
        
        function A = rightnull(A, alg)
            arguments
                A
                alg = 'svd'
            end
            
            if numel(A) > 1
                for i = numel(A):-1:1
                    A(i) = rightnull(A(i), alg);
                end
                return
            end
            
            A.var = rightnull(A.var, 1, 2:nspaces(A), alg);
        end
        
        function d = dot(A, B)
            if isa(A, 'MpsTensor')
                A = A.var;
            end
            if isa(B, 'MpsTensor')
                B = B.var;
            end
            d = dot(A, B);
        end
        
        function d = overlap(A, B)
            arguments
                A
                B MpsTensor
            end
            
            inds = 1:nspaces(A);
            
            d = contract(A, inds, B, [fliplr(inds(1:end-B.alegs)) inds(end-B.alegs+1:end)]);
        end
        
        function C = mtimes(A, B)
            if isnumeric(A)
                C = B;
                C.var = A * C.var;
                return
            end
            
            if isnumeric(B)
                C = A;
                C.var = C.var * B;
                return
            end
            
            error('not implemented when both inputs not numeric.');
        end
        
        function A = repartition(A, varargin)
            for i = 1:numel(A)
                A(i).var = repartition(A(i).var, varargin{:});
            end
        end
        
        function A = plus(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    alegs = varargin{i}.alegs;
                    varargin{i} = varargin{i}.var;
                end
            end
            A = MpsTensor(plus(varargin{:}), alegs);
        end
        
        function A = minus(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    alegs = varargin{i}.alegs;
                    varargin{i} = varargin{i}.var;
                end
            end
            A = MpsTensor(minus(varargin{:}), alegs);
        end
        
        function A = rdivide(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    alegs = varargin{i}.alegs;
                    varargin{i} = varargin{i}.var;
                end
            end
            A = MpsTensor(rdivide(varargin{:}), alegs);
        end
        
        function A = ldivide(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    alegs = varargin{i}.alegs;
                    varargin{i} = varargin{i}.var;
                end
            end
            A = MpsTensor(ldivide(varargin{:}), alegs);
        end
        
        function n = norm(A)
            n = norm([A.var]);
        end
        
        function A = conj(A)
            [A.var] = conj([A.var]);
        end

        function A = twist(A, varargin)
            [A.var] = twist([A.var], varargin{:});
        end
        
        function t = ctranspose(t)
            % Compute the adjoint of a tensor. This is defined as swapping the codomain and
            % domain, while computing the adjoint of the matrix blocks.
            %
            % Usage
            % -----
            % :code:`t = ctranspose(t)`
            % :code:`t = t'`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   adjoint tensor.
            
            for i = 1:numel(t)
                t(i).var = t(i).var';
            end
            
            t = permute(t, ndims(t):-1:1);
        end
        
        function A = tpermute(A, varargin)
            A.var = tpermute(A.var, varargin{:});
        end
        
        function C = tensorprod(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    varargin{i} = varargin{i}.var;
                end
            end
            C = tensorprod(varargin{:});
        end
        
        function [AL, CL, lambda, eta] = uniform_leftorth(A, CL, kwargs)
            arguments
                A
                CL = []
                kwargs.Tol = eps(underlyingType(A))^(3/4)
                kwargs.MaxIter = 500
                kwargs.Method = 'polar'
                kwargs.Verbosity = 0
                kwargs.Normalize = true
                kwargs.EigsInit = 3
                kwargs.EigsFrequence = 2
            end
%             
            % constants
            EIG_TOLFACTOR = 1/50;
            EIG_MAXTOL = 1e-4;
            MINKRYLOVDIM = 8;
            MAXKRYLOVDIM = 30;
            
            % initialization
            N = size(A, 2);
            if isempty(CL), CL = initializeC(A, A); end
            if kwargs.Normalize, CL(1) = normalize(CL(1)); end
            for i = 1:numel(A)
                A(i).var = repartition(A(i).var, [nspaces(A(i).var) - 1, 1]);
            end
            AL = A;
            
            eta_best = Inf;
            ctr_best = 0;
            AL_best = AL;
            C_best = CL;
            lambda_best = 0;
            
            for ctr = 1:kwargs.MaxIter
                if ctr > kwargs.EigsInit && mod(ctr, kwargs.EigsFrequence) == 0
                    C_ = repartition(CL(end), [2 0]);
                    T = transfermatrix(A, AL);
                    [C_, ~] = eigsolve(T, C_, ...
                        1, 'largestabs', ...
                        'Tol', min(eta * EIG_TOLFACTOR, EIG_MAXTOL), ...
                        'KrylovDim', between(MINKRYLOVDIM, ...
                            MAXKRYLOVDIM - ctr / 2 + 4, ...
                            MAXKRYLOVDIM), ...
                        'NoBuild', 4, ...
                        'Verbosity', kwargs.Verbosity - 1);
                    [~, CL(end)] = leftorth(C_, 1, 2, kwargs.Method);
                    if isdual(space(CL(end), 2)) == isdual(space(CL(end), 1))
                        CL(end).codomain = conj(CL(end).codomain);
                        CL(end) = twist(CL(end), 1);
                    end
                end
                
                C_ = CL(end);
                lambdas = ones(1, N);
                for w = 1:N
                    ww = prev(w, N);
                    CA = multiplyleft(A(w), CL(ww));
                    [AL(w), CL(w)] = leftorth(CA, kwargs.Method);
                    lambdas(w) = norm(CL(w));
                    if kwargs.Normalize, CL(w) = CL(w) ./ lambdas(w); end
                end
                lambda = prod(lambdas);
                eta = norm(C_ - CL(end), Inf);
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
        
        function [AR, CR, lambda, eta] = uniform_rightorth(A, CR, kwargs)
            arguments
                A
                CR = []
                kwargs.Tol = eps(underlyingType(A))^(3/4)
                kwargs.MaxIter = 500
                kwargs.Method = 'polar'
                kwargs.Verbosity = 0
                kwargs.Normalize = true
                kwargs.EigsInit = 3
                kwargs.EigsFrequence = 2
            end
            
            opts = namedargs2cell(kwargs);
            
            Ad = flip(arrayfun(@ctranspose, A));
            if isempty(CR)
                Cd = [];
            else
                Cd = circshift(flip(arrayfun(@ctranspose, CR)), -1);
            end
            
            [ARd, CRd, lambda, eta] = uniform_leftorth(Ad, Cd, opts{:});
            
            AR = flip(arrayfun(@ctranspose, ARd));
            CR  = circshift(flip(arrayfun(@ctranspose, CRd)), -1);
            lambda = conj(lambda);
        end
        
        function T = transfermatrix(A, B)
            arguments
                A MpsTensor
                B MpsTensor = A
            end
            
            for i = length(A):-1:1
                B(i).var = twist(B(i).var, [isdual(space(B(i), 1:2)) ~isdual(space(B(i), 3))]);
                T(i, 1) = FiniteMpo(B(i).var', {}, A(i));
            end
        end
        
        function v = applytransfer(L, R, v)
            arguments
                L MpsTensor
                R MpsTensor
                v = []
            end
            
            if isempty(v)
                v = tracetransfer(L, R);
                return
            end
            
            auxlegs_v = nspaces(v) - 2;
            auxlegs_l = L.alegs;
            auxlegs_r = R.alegs;
            newrank = rank(v); newrank(2) = newrank(2) + auxlegs_l + auxlegs_r;
            
            v = contract( ...
                v, [1, R.plegs+2, (-(1:auxlegs_v) - 2 - auxlegs_l)], ...
                L, [-1, (1:L.plegs)+1, 1 (-(1:auxlegs_l) - (1+L.plegs))], ...
                R, [R.plegs+2, (R.plegs:-1:1)+1, -2, (-(1:auxlegs_r) - (2+R.plegs) - auxlegs_l - auxlegs_v)], ...
                'Rank', newrank);
        end
        
        function v = contracttransfer(L, R)
            arguments
                L MpsTensor
                R MpsTensor
            end
            
            auxlegs_l = L.alegs;
            auxlegs_r = R.alegs;
            assert(R.plegs == L.plegs);
            plegs = L.plegs; %#ok<PROPLC>
            
            v = contract(...
                L, [-1, 1:plegs, -4, (-(1:auxlegs_l) - 4)], ...
                R, [-3, flip(1:plegs), -2 (-(1:auxlegs_r) - 4 - auxlegs_l)], ...
                'Rank', [2, 2] + [0, auxlegs_l + auxlegs_r]); %#ok<PROPLC>
        end
        
        function v = tracetransfer(L, R)
            arguments
                L MpsTensor
                R MpsTensor
            end
            
            auxlegs_l = L.alegs;
            auxlegs_r = R.alegs;
            assert(R.plegs == L.plegs);
            plegs = L.plegs; %#ok<PROPLC>
            newrank = [1 1 + auxlegs_l + auxlegs_r];
            
            v = contract(...
                L, [-1 1:(plegs + 1) (-(1:auxlegs_l) - 2)], ...
                R, [flip(1:(plegs + 1)) -2 (-(1:auxlegs_r) - 2 - auxlegs_l)], ...
                'Rank', newrank); %#ok<PROPLC>
        end
        
        function rho = applyleft(T, B, rho)
            if nargin == 2
                rho = [];
                return
            end
            assert(isequal(size(T), size(B)), 'mpstensor:dimagree', ...
                'dimensions should agree.');
            if length(T) > 1
                for i = 1:length(T)
                    rho = applyleft(T(i), B(i), rho);
                end
                return
            end
            
            if isempty(rho)
                switch num2str([T1.plegs T1.alegs T2.alegs])
                    case '1  0  0'
                        indices = {[1 2 -2] [1 2 -1]};
                        
                    case '2  0  0'
                        indices = {[1 2 3 -2] [1 2 3 -1]};
                    case '1  1  0'
                        indices = {[1 2 -2 -3] [1 2 -1]};
                    case '1  0  1'
                        indices = {[1 2 -2] [1 2 -1 -3]};
                    case '1  1  1'
                        indices = {[1 2 -2 3] [1 2 -1 3]};
                    otherwise
                        error('mps:ArgError', 'leftapply ill-defined.');
                end
                rho = contract(T1, indices{1}, T2, indices{2});
                return
            end
            
            switch num2str([nspaces(rho) T.plegs T.alegs B.alegs])
                case '2  1  0  0'
                    tmp = repartition(rho, [1 1]) * repartition(T, [1 2]);
                    tmp = tpermute(B, [3 2 1], [1 2]) * repartition(tmp, [2 1]);
                    rho = repartition(tmp, rank(rho));
%                     indices = {[1 2] [2 3 -2] [1 3 -1]};
%                     r = rank(rho);
%                     T = twist(T, [isdual(space(T, [1 2])) false]);
                case '3  1  0  0'
                    tmp = tpermute(rho, [1 3 2], [2 1]) * repartition(T, [1 2]);
                    tmp = tpermute(B, [3 2 1], [1 2]) * tpermute(tmp, [1 3 4 2], [2 2]);
                    rho = repartition(tmp, rank(rho));
                otherwise
                    error('mps:ArgError', 'applyleft ill-defined');
            end
%             rho = contract(rho, indices{1}, T, indices{2}, B, indices{3}, ...
%                 'Rank', r);
        end
        
        function rho = applyright(T, B, rho)
            if nargin == 2
                rho = [];
                return
            end
            assert(isequal(size(T), size(B)), 'mpstensor:dimagree', ...
                'dimensions should agree.');
            if length(T) > 1
                for i = length(T):-1:1
                    rho = applyright(T(i), B(i), rho);
                end
                return
            end
            if isempty(rho)
                switch num2str([T.plegs T.alegs B.alegs])
                    case '1  0  0'
                        rho = contract(T, [-1 2 1], B, [-2 2 1]);
%                         rhoR=Contract({T1,T2},{[-1,2,1],[-2,2,1]});
                    case '2  0  0'
                        rhoR=Contract({T1,T2},{[-1,2,3,1],[-2,2,3,1]});
                    case '1  1  0'
                        rhoR=Contract({T1,T2},{[-1,2,1,-3],[-2,2,1]});
                    case '1  0  1'
                        rhoR=Contract({T1,T2},{[-1,2,1],[-2,2,1,-3]});
                    case '1  1  1'
                        rhoR=Contract({T1,T2},{[-1,2,1,3],[-2,2,1,3]});
                    otherwise
                        error('mps:ArgError', 'applyright ill-defined.');
                end
                return
            end
            
            switch num2str([nspaces(rho) T.plegs T.alegs B.alegs])
                case '2  1  0  0'
                    tmp = repartition(T, [2 1]) * repartition(rho, [1 1]);
                    rho = repartition(...
                        repartition(tmp, [1 2]) * tpermute(B, [3 2 1], [2 1]), rank(rho));
                    
%                     T = twist(T, [false ~isdual(space(T, 2:3))]);
%                     rho = contract(rho, [1 2], T, [-1 3 1], B, [-2 3 2], ...
%                         'Rank', rank(rho));
                otherwise
                    error('mps:ArgError', 'applyright ill-defined.');
            end
        end
        
        function A = multiplyleft(A, C)
            [A.var] = contract(C, [-1 1], ...
                [A.var], [1 -((1:A.plegs)+1) -((1:A.alegs+1)+A.plegs+1)], 'Rank', rank(A));
        end
        
        function A = multiplyright(A, C)
            [A.var] = contract([A.var], [-(1:A.plegs+1) 1 -((1:A.alegs)+A.plegs+2)], ...
                C, [1 -(2+A.plegs)], 'Rank', rank(A));
        end
        
        function C = initializeC(AL, AR)
            for i = length(AL):-1:1
                C(i) = AL(i).var.eye(rightvspace(AL(i))', leftvspace(AR(next(i, length(AR)))));
            end
        end
        
        function A = expand(A, addspace, noisefactor)
            arguments
                A
                addspace
                noisefactor = 1e-3
            end
            
            for i = length(A):-1:1
                spaces = space(A(i));
                spaces(nspaces(A(i)) - A(i).alegs) = addspace(i);
                spaces(1) = conj(addspace(prev(i, length(A))));
                r = rank(A(i));
                A(i).var = embed(A(i).var, ...
                    noisefactor * ...
                    normalize(A(i).var.randnc(spaces(1:r(1)), spaces(r(1)+1:end)')));
            end
        end
        
        function type = underlyingType(A)
            type = underlyingType(A.var);
        end
        
        function disp(t)
            builtin('disp', t);
        end
        
        function n = nnz(t)
            n = nnz(t.var);
        end
    end
    
    
    %% Converters
    methods
        function t = Tensor(A)
            for i = numel(A):-1:1
                t(i) = full(A(i).var);
            end
        end
        
        function t = SparseTensor(A)
            t = reshape([A.var], size(A));
            t = sparse(t);
        end
    end
    
    
    %% Solvers
    methods
        function v = vectorize(t, type)
            % Collect all parameters in a vector, weighted to reproduce the correct
            % inproduct.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % type : 'real' or 'complex'
            %   optionally specify if complex entries should be seen as 1 or 2 parameters.
            %   Defaults to 'complex', with complex parameters.
            %
            % Returns
            % -------
            % v : numeric
            %   real or complex vector containing the parameters of the tensor.
            
            arguments
                t
                type = 'complex'
            end
            
            v = vectorize(t.var, type);
        end
        
        function t = devectorize(v, t, type)
            % Collect all parameters from a vector, and insert into a tensor.
            %
            % Arguments
            % ---------
            % v : numeric
            %   real or complex vector containing the parameters of the tensor.
            %
            % t : :class:`Tensor`
            %   input tensor.
            %
            % type : 'real' or 'complex'
            %   optionally specify if complex entries should be seen as 1 or 2 parameters.
            %   Defaults to 'complex', with complex parameters.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor, filled with the parameters.
            
            arguments
                v
                t
                type = 'complex'
            end
            
            t.var = devectorize(v, t.var, type);
        end
    end
    
    
    %%
    methods (Static)
        function local_tensors = decompose_local_state(psi, kwargs)
            % convert a tensor into a product of local operators.
            %
            % Usage
            % -----
            % :code:`local_operators = MpoTensor.decompose_local_operator(H, kwargs)`.
            %
            % Arguments
            % ---------
            % H : :class:`AbstractTensor`
            %   tensor representing a local operator on N sites.
            %
            % Keyword Arguments
            % -----------------
            % 'Trunc' : cell
            %   optional truncation method for the decomposition. See also
            %   :method:`Tensor.tsvd`
            arguments
                psi
                kwargs.Trunc = {'TruncBelow', 1e-14}
            end
            
            L = nspaces(psi);
            assert(L >= 3, 'argerror', ...
                sprintf('state must have at least 3 legs. (%d)', L));
            
            local_tensors = cell(1, L - 2);
            for i = 1:length(local_tensors)-1
                [u, s, psi] = tsvd(psi, [1 2], [3:nspaces(psi)], kwargs.Trunc{:});
                local_tensors{i} = multiplyright(MpsTensor(u), s);
            end
            local_tensors{end} = MpsTensor(repartition(psi, [2 1]));
        end
    end
end

