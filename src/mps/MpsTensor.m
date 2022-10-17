classdef MpsTensor < Tensor
    %MPSTENSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        plegs = 1
        alegs = 0
    end
    
    
    %% Constructors
    methods
        function A = MpsTensor(tensor, alegs)
            arguments
                tensor = []
                alegs = zeros(size(tensor))
            end
            
            if isscalar(alegs) && ~isscalar(tensor)
                alegs = repmat(alegs, size(tensor));
            else
                assert(isequal(size(alegs, 1:ndims(tensor)), size(tensor)), 'mps:ArgError', ...
                    'Input sizes incompatible.');
            end
            
            if isempty(tensor)
                args = {};
            else
                args = {tensor};
            end
            A@Tensor(args{:});
            if ~isempty(tensor)
                for i = numel(A):-1:1
                    A(i).plegs = nspaces(tensor(i)) - alegs(i) - 2;
                    A(i).alegs = alegs(i);
                end
            end
        end
    end
    
    
    %% Properties
    methods
        function s = pspace(A)
            s = space(A, 1 + (1:A.plegs));
        end
        
        function s = leftvspace(A)
            s = space(A, 1);
        end
        
        function s = rightvspace(A)
            s = space(A, nspaces(A) - A.alegs);
        end
    end
    
    
    %% Linear Algebra
    methods
        function [AL, L] = leftorth(A, alg)
            arguments
                A
                alg = 'qr'
            end
            
            if A.alegs == 0
                [AL, L] = leftorth@Tensor(A, 1:nspaces(A)-1, nspaces(A), alg);
            else
                [AL, L] = leftorth@Tensor(A, [1:A.plegs+1 A.plegs+3], A.plegs+2, alg);
                AL = permute(AL, [1:A.plegs+1 A.plegs+3 A.plegs+2], rank(A));
            end
            AL = MpsTensor(AL, A.alegs);
        end
        
        function [AL, C, lambda, eta] = uniform_leftorth(A, C, kwargs)
            arguments
                A
                C = initializeC(A, circshift(A, 1))
                kwargs.Tol = eps(underlyingType(A))^(3/4)
                kwargs.MaxIter = 500
                kwargs.Method = 'polar'
                kwargs.Verbosity = 0
                kwargs.Normalize = true
                kwargs.EigsInit = 3
                kwargs.EigsFrequence = 2
            end
            
            % constants
            EIG_TOLFACTOR = 1/50;
            EIG_MAXTOL = 1e-4;
            MINKRYLOVDIM = 8;
            MAXKRYLOVDIM = 30;
            
            % initialization
            N = width(A);
            if isempty(C), C = initializeC(A, circshift(A, 1)); end
            if kwargs.Normalize, C(1) = normalize(C(1)); end
            A = arrayfun(@(a) MpsTensor(repartition(a, [nspaces(a)-1 1])), A);
            AL = A;
            
            for ctr = 1:kwargs.MaxIter
                if ctr > kwargs.EigsInit && mod(ctr, kwargs.EigsFrequence) == 0
                    C_ = repartition(C(end), [2 0]);
                    [C_, ~] = eigsolve(@(x) applyleft(A, conj(AL), x), C_, ...
                        1, 'largestabs', ...
                        'Tol', min(eta * EIG_TOLFACTOR, EIG_MAXTOL), ...
                        'KrylovDim', between(MINKRYLOVDIM, ...
                            MAXKRYLOVDIM - ctr / 2 + 4, ...
                            MAXKRYLOVDIM), ...
                        'NoBuild', 4, ...
                        'Verbosity', kwargs.Verbosity - 1);
                    [~, C(end)] = leftorth(C_, 1, 2, kwargs.Method);
                end
                
                C_ = C(end);
                lambdas = ones(1, N);
                for w = 1:N
                    ww = prev(w, N);
                    CA = multiplyleft(A(w), C(ww));
                    [AL(w), C(w)] = leftorth(CA, kwargs.Method);
                    lambdas(w) = norm(C(w));
                    if kwargs.Normalize, C(w) = C(w) ./ lambdas(w); end
                end
                eta = norm(C_ - C(end), Inf);
                
                if eta < kwargs.Tol
                    if kwargs.Verbosity >= Verbosity.conv
                        fprintf('Conv %2d:\terror = %0.4e\n', ctr, eta);
                    end
                    break;
                else
                    if kwargs.Verbosity >= Verbosity.iter
                        fprintf('Iter %2d:\terror = %0.4e\n', ctr, eta);
                    end
                end
            end
            
            if kwargs.Verbosity >= Verbosity.warn && eta > kwargs.Tol
                fprintf('Not converged %2d:\terror = %0.4e\n', ctr, eta);
            end
            
            if nargout > 2, lambda = prod(lambdas); end
        end
        
        function [R, AR] = rightorth(A, alg)
            arguments
                A
                alg = 'lq'
            end
            
            [R, AR] = rightorth@Tensor(A, 1, 2:nspaces(A), alg);
            AR = MpsTensor(repartition(AR, rank(A)), A.alegs);
        end
        
        function [AR, C, lambda, eta] = uniform_rightorth(A, C, kwargs)
            arguments
                A
                C = initializeC(A, circshift(A, 1))
                kwargs.Tol = eps(underlyingType(A))^(3/4)
                kwargs.MaxIter = 500
                kwargs.Method = 'polar'
                kwargs.Verbosity = 0
                kwargs.Normalize = true
                kwargs.EigsInit = 3
                kwargs.EigsFrequence = 2
            end
            
            % constants
            EIG_TOLFACTOR = 1/50;
            EIG_MAXTOL = 1e-4;
            MINKRYLOVDIM = 8;
            MAXKRYLOVDIM = 30;
            
            % initialization
            N = width(A);
            if isempty(C), C = initializeC(A, circshift(A, 1)); end
            if kwargs.Normalize, C(1) = normalize(C(1)); end
            A = arrayfun(@(a) MpsTensor(repartition(a, [1 nspaces(a)-1])), A);
            AR = A;
            
            for ctr = 1:kwargs.MaxIter
                if ctr > kwargs.EigsInit && mod(ctr, kwargs.EigsFrequence) == 0
                    C_ = repartition(C(end), [nspaces(C(end)) 0]);
                    [C_, ~] = eigsolve(@(x) applyright(A, conj(AR), x), C_, ...
                        1, 'largestabs', ...
                        'Tol', min(eta * EIG_TOLFACTOR, EIG_MAXTOL), ...
                        'KrylovDim', between(MINKRYLOVDIM, ...
                            MAXKRYLOVDIM - ctr / 2 + 4, ...
                            MAXKRYLOVDIM), ...
                        'NoBuild', 4, ...
                        'Verbosity', kwargs.Verbosity - 1);
                    [C(end), ~] = rightorth(C_, 1, 2, kwargs.Method);
                end
                
                C_ = C(end);
                lambdas = ones(1, N);
                for w = N:-1:1
                    ww = prev(w, N);
                    AC = multiplyright(A(w), C(w));
                    [C(ww), AR(w)] = rightorth(AC, kwargs.Method);
                    lambdas(w) = norm(C(ww));
                    if kwargs.Normalize, C(ww) = C(ww) ./ lambdas(w); end
                end
                eta = norm(C_ - C(end), Inf);
                
                if eta < kwargs.Tol
                    if kwargs.Verbosity >= Verbosity.conv
                        fprintf('Conv %2d:\terror = %0.4e\n', ctr, eta);
                    end
                    break;
                else
                    if kwargs.Verbosity >= Verbosity.iter
                        fprintf('Iter %2d:\terror = %0.4e\n', ctr, eta);
                    end
                end
            end
            
            if kwargs.Verbosity >= Verbosity.warn && eta > kwargs.Tol
                fprintf('Not converged %2d:\terror = %0.4e\n', ctr, eta);
            end
            
            if nargout > 2, lambda = prod(lambdas); end
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
            A = MpsTensor(repartition(...
                repartition(C, [1 1]) * repartition(A, [1 nspaces(A)-1]), ...
                rank(A)), A.alegs);
%             if ~isdual(space(C, 2))
%                 C = twist(C, 2);
%             end
%             A = MpsTensor(contract(C, [-1 1], A, [1 -(2:nspaces(A))], 'Rank', rank(A)));
        end
        
        function A = multiplyright(A, C)
            if A.alegs == 0
                A = MpsTensor(repartition(...
                    repartition(A, [nspaces(A)-1 1]) * repartition(C, [1 1]), ...
                    rank(A)), 0);
                return
            end
            if isdual(space(C, 1)), C = twist(C, 1); end
            Alegs = nspaces(A);
            if A.alegs == 0
                A = contract(A, [-(1:Alegs-1) 1], C, [1 -Alegs], 'Rank', rank(A));
            else
                A = contract(A, [-(1:Alegs-1-A.alegs) 1 -(Alegs-A.legs+1:Alegs)], ...
                    C, [1 Alegs - A.alegs], 'Rank', rank(A));
            end
            A = MpsTensor(A);
        end
        
        function C = initializeC(AL, AR)
            for i = length(AL):-1:1
                C(i) = AL.eye(rightvspace(AL(i))', leftvspace(AR(i)));
            end
        end
        
        function A = expand(A, addspace)
            NOISE_FACTOR = 1e-3;
            for i = length(A):-1:1
                spaces = space(A(i));
                spaces(1) = addspace(i);
                spaces(nspaces(A(i)) - A(i).alegs) = conj(addspace(next(i, length(A))));
                r = rank(A(i));
                A(i) = embed(A(i), ...
                    NOISE_FACTOR * ...
                    normalize(A.randnc(spaces(1:r(1)), spaces(r(1)+1:end)')));
            end
        end
    end
end

