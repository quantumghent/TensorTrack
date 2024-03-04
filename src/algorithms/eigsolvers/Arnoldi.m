classdef Arnoldi
    % Arnoldi Krylov algorithm for linear algebra problems.
    %
    % Properties
    % ----------
	% tol : :class:`double`
    %   convergence tolerance, defaults to :code:`1e-10`.
    %
	% maxiter : :class:`int`
    %   maximum number of iterations, defaults to :code:`100`.
    %
	% krylovdim : :class:`int`
    %   Krylov subspace dimension, defaults to :code:`20`.
    %
	% deflatedim : :class:`int`
    %   number of Krylov vectors to keep when deflating.
    %
	% reorth : :class:`int`
    %   reorthogonalize basis if larger than this number, defaults to :code:`20`.
    %
	% nobuild : :class:`int`
    %   frequency of convergence check when building, defaults to :code:`3`.
    %
	% verbosity : :class:`.Verbosity`
    %   display information, defaults to :code:`Verbosity.warn`.
    
    properties
        tol         = 1e-10
        maxiter     = 100
        krylovdim   = 20
        deflatedim
        reorth      = 20
        nobuild     = 3
        verbosity   = Verbosity.warn
    end
    
    methods
        function alg = Arnoldi(kwargs)
            arguments
                kwargs.?Arnoldi
            end
            fields = fieldnames(kwargs);
            if ~isempty(fields)
                for field = fields.'
                    alg.(field{1}) = kwargs.(field{1});
                end
            end
        end
        
        function varargout = eigsolve(alg, A, v0, howmany, sigma)
            % Find a few eigenvalues and eigenvectors of an operator using an Arnoldi
            % routine.
            %
            % Usage
            % -----
            % :code:`[V, D, flag] = eigsolve(A, v, howmany, sigma)`
            %
            % :code:`D = eigsolve(A, v, ...)`
            %
            % Arguments
            % ---------
            % A : :class:`matrix` or :class:`function_handle`
            %   A square matrix.
            %   A function handle which implements one of the following, depending on sigma:
            %
            %   - :code:`A \ x`, if `sigma` is 0 or 'smallestabs'
            %   - :code:`(A - sigma * I) \ x`, if sigma is a nonzero scalar
            %   - :code:`A * x`, for all other cases
            %
            % v : :class:`vector`
            %   initial guess for the eigenvector.
            %
            % howmany : :class:`int`
            %   amount of eigenvalues and eigenvectors that should be computed. By default
            %   this is 1, and this should not be larger than the total dimension of A.
            %
            % sigma : :class:`char` or numeric
            %   selector for the eigenvalues, should be either one of the following:
            %
            %   - 'largestabs', 'lm': default, eigenvalues of largest magnitude
            %   - 'largestreal', 'lr': eigenvalues with largest real part
            %   - 'largestimag', 'li': eigenvalues with largest imaginary part.
            %   - 'smallestabs', 'sm': default, eigenvalues of smallest magnitude
            %   - 'smallestreal', 'sr': eigenvalues with smallest real part
            %   - 'smallestimag', 'si': eigenvalues with smallest imaginary part.
            %   - numeric : eigenvalues closest to sigma.
            %
            % Returns
            % -------
            % V : (1, howmany) :class:`vector`
            %   vector of eigenvectors.
            %
            % D : :class:`numeric`
            %   vector of eigenvalues if only a single output argument is asked, diagonal
            %   matrix of eigenvalues otherwise.
            %
            % flag : :class:`int`
            %   convergence info flag:
            %
            %   - flag = 0: all eigenvalues are converged.
            %   - flag = 1: invariant subspace was found and the algorithm was aborted.
            %   - flag = 2: algorithm did not converge after maximum number of iterations.

            arguments
                alg
                A
                v0
                howmany = 1
                sigma = 'lm'
            end
            
            t_total = tic;
            
            if isnumeric(v0)
                v0_vec = v0;
                if isa(A, 'function_handle')
                    A_fun = @A;
                else
                    A_fun = @(v) A * v;
                end
            else
                v0_vec = vectorize(v0);
                if isa(A, 'function_handle')
                    A_fun = @(v) vectorize(A(devectorize(v, v0)));
                else
                    A_fun = @(v) vectorize(A * devectorize(v, v0));
                end
            end
            
            if norm(v0_vec) < eps(underlyingType(v0_vec))^(3/4)
                error('eigsolve:inputnorm', 'starting vector should not have zero norm.');
            end

            sz = size(v0_vec);
            
            if sz(1) < howmany
                howmany = sz(1);
                if alg.verbosity >= Verbosity.warn
                    warning('eigsolve:size', 'requested %d out of %d eigenvalues.', ...
                        howmany, sz(1));
                end
            end
            
            if sz(1) < alg.krylovdim
                alg.krylovdim = sz(1);
                if alg.verbosity >= Verbosity.warn
                    warning('eigsolve:size', ...
                        'Krylov subspace dimension is larger than total number of eigenvalues, reducing Krylov dimension to %d.', ...
                        sz(1));
                end
            end
            
            % some input validation
            if isempty(alg.nobuild), alg.nobuild = ceil(alg.krylovdim / 10); end
            if isempty(alg.deflatedim), alg.deflatedim = max(round(3/5 * alg.krylovdim), howmany); end
            alg.deflatedim = max(alg.deflatedim, howmany);
            assert(alg.deflatedim < alg.krylovdim, 'eigsolve:argerror', ...
                'Deflate size should be smaller than krylov dimension.')
            
            v = v0_vec / norm(v0_vec, 'fro');
            
            % preallocation
            V = zeros(length(v), alg.krylovdim, 'like', v); % Krylov subspace basis
            H = zeros(alg.krylovdim, alg.krylovdim, underlyingType(v)); % Hessenberg matrix
            
            ctr_outer = 0;
            ctr_inner = 0;
            flag = 0;
            
            while ctr_outer < alg.maxiter
                t_outer = tic;
                ctr_outer = ctr_outer + 1;
                
                flag_inner = 0;
                while ctr_inner < alg.krylovdim  % build Krylov subspace
                    t_inner = tic;
                    ctr_inner = ctr_inner + 1;
                    
                    V(:, ctr_inner) = v;
                    v = A_fun(v);
                    for i = 1:ctr_inner
                        H(i, ctr_inner) = dot(V(:, i), v);
                    end
                    v = v - V * H(:, ctr_inner);
                    
                    % reorthogonalize new vector
                    if ctr_inner >= alg.reorth
                        if isnumeric(V)
                            c = V' * v;
                        else
                            c = zeros(size(V, 2), 1);
                            for i = 1:ctr_inner
                                c(i) = dot(V(:, i), v);
                            end
                        end
                        H(:, ctr_inner) = H(:, ctr_inner) + c;
                        v = v - V * c;
                    end
                    
                    % normalize
                    beta = norm(v, 'fro');
                    v = v ./ beta;
                    
                    if ctr_inner == alg.krylovdim, break; end
                    
                    if ctr_inner >= howmany
                        invariantsubspace = beta < eps(underlyingType(beta))^(3/4);
                        if invariantsubspace || ctr_inner == alg.krylovdim
                            break;
                        end
                        
                        % check for convergence during subspace build
                        if ~mod(ctr_inner, alg.nobuild)
                            [U, lambda] = eig(H(1:ctr_inner, 1:ctr_inner), 'vector');
                            select = selecteigvals(lambda, howmany, sigma);
                            conv = beta * max(abs(U(ctr_inner, select)));
                            
                            if conv < alg.tol
                                V = V(:, 1:ctr_inner) * U(:, select);
                                D = diag(lambda(select));
                                if alg.verbosity >= Verbosity.conv
                                    fprintf('Conv %2d (%2d/%2d):\tlambda = %.5e + %.5ei;\terror = %.5e;\ttime = %s.\n', ...
                                        ctr_outer, ctr_inner, alg.krylovdim, ...
                                         real(lambda(1)), imag(lambda(1)), conv, ...
                                         time2str(toc(t_total)));
                                end
                                flag_inner = 1;
                                break
                            end
                            
                            if alg.verbosity >= Verbosity.detail
                                fprintf('Iter %2d (%2d/%2d):\tlambda = %.5e + %.5ei;\terror = %.5e;\ttime = %s.\n', ...
                                    ctr_outer, ctr_inner, alg.krylovdim, ...
                                    real(lambda(1)), imag(lambda(1)), conv, ...
                                    time2str(toc(t_inner)));
                            end
                        end
                    end
                    
                    H(ctr_inner + 1, ctr_inner) = beta;
                end
                if flag_inner
                    break
                end
                
                % stopping criterium reached - irrespective of convergence
                if ctr_outer == alg.maxiter || ctr_inner ~= alg.krylovdim
                    [U, lambda] = eig(H(1:ctr_inner, 1:ctr_inner), 'vector');
                    select = selecteigvals(lambda, howmany, sigma);
                    conv = max(abs(beta * U(end, select)));
                    V = V(:, 1:ctr_inner) * U(:, select);
                    D = diag(lambda(select));
                    
                    if conv > alg.tol
                        if invariantsubspace
                            flag = 1;
                        else
                            flag = 2;
                        end
                    end
                    break
                end
                
                % deflate Krylov subspace
                [U1, T] = schur(H, 'real');
                E = ordeig(T);
                select1 = false(size(E));
                select1(selecteigvals(E, alg.deflatedim, sigma)) = true;
                [U1, T] = ordschur(U1, T, select1);
                
                V = V * U1;
                [U, lambda] = eig(T(1:alg.deflatedim, 1:alg.deflatedim), 'vector');
                select = selecteigvals(lambda, howmany, sigma);

                conv = max(abs(beta * U1(alg.krylovdim, 1:alg.deflatedim) * U(:, select)));
                
                % check for convergence
                if conv < alg.tol
                    V = V(:, 1:alg.deflatedim) * U(:, select);
                    D = diag(lambda(select));
                    if alg.verbosity >= Verbosity.conv
                        fprintf('Conv %2d:\tlambda = %.5e + %.5ei;\terror = %.5e;\ttime = %s.\n', ...
                            ctr_outer, real(lambda(1)), imag(lambda(1)), conv, time2str(toc(t_outer)));
                    end
                    break
                end
                
                if alg.verbosity >= Verbosity.iter
                    fprintf('Iter %2d:\tlambda = %.5e + %.5ei;\terror = %.5e;\ttime = %s.\n', ...
                        ctr_outer, real(lambda(1)), imag(lambda(1)), conv, time2str(toc(t_total)));
                end
                
                % deflate Krylov subspace
                H = zeros(alg.krylovdim, alg.krylovdim, underlyingType(v));
                H(1:alg.deflatedim, 1:alg.deflatedim) = ...
                    T(1:alg.deflatedim, 1:alg.deflatedim);
                H(alg.deflatedim + 1, 1:alg.deflatedim) = ...
                    beta * U1(alg.krylovdim, 1:alg.deflatedim);
                V(:, alg.deflatedim + 1:end) = 0 * V(:, alg.deflatedim + 1:end);
                ctr_inner = alg.deflatedim;
            end
            
            % process results
            if nargout <= 1
                varargout = {diag(D)};
            else
                if ~isnumeric(v0)
                    for i = howmany:-1:1
                        varargout{1}(:, i) = devectorize(V(:, i), v0);
                    end
                else
                    for i = howmany:-1:1
                        varargout{1}(:, i) = V(:, i);
                    end
                end
                varargout{2} = D;
                if nargout == 3
                    varargout{3} = flag;
                end
            end
            
            % display
            if nargout < 3
                if flag == 1
                    warning('Found invariant subspace (error = %.5e).\n', conv);
                elseif flag == 2
                     warning('Reached maxiter without convergence.\n');
                end
            end
            if ~flag && alg.verbosity > Verbosity.warn
                fprintf('eigsolve converged.\n');
            end
        end
        
        function [w0, flag] = expsolve(alg, A, t, v0)
            if isempty(alg.nobuild), alg.nobuild = ceil(alg.krylovdim / 10); end
            if isempty(alg.deflatedim), alg.deflatedim = max(round(3/5 * alg.krylovdim), howmany); end
            
            beta0 = norm(v0, 'fro');
            w0 = v0;
            
            w = w0;
            v = w0;
            beta = norm(w, 'fro');
            v = w / beta;
            
            % preallocation
            V = zeros(1, alg.krylovdim, 'like', v0);                 % Krylov subspace basis
            H = zeros(alg.krylovdim, alg.krylovdim, underlyingType(v0)); % Hessenberg matrix
            
            % time step parameters
            eta = alg.tol;
            totalerr = 0;
            sgn = sign(t);
            tau = abs(t);
            tau0 = 0;
            dtau = tau - tau0;
            
            delta = 1.2;
            gamma = 0.8;
            
            ctr_outer = 0;
            ctr_inner = 0;
            flag = 0;
            
            while ctr_outer < alg.maxiter
                ctr_outer = ctr_outer + 1;
                
                while ctr_inner < alg.krylovdim  % build Krylov subspace
                    ctr_inner = ctr_inner + 1;
                    
                    V(1:length(w0), ctr_inner) = w;
                    w = A(w);
                    for i = 1:ctr_inner
                        H(i, ctr_inner) = dot(V(:, i), w);
                    end
                    w = w - V * H(:, ctr_inner);
                    
                    % reorthogonalize new vector
                    if ctr_inner >= alg.reorth
                        c = zeros(size(V, 2), 1);
                        for i = 1:ctr_inner
                            c(i) = dot(V(:, i), w);
                        end
                        H(:, ctr_inner) = H(:, ctr_inner) + c;
                        w = w - V * c;
                    end
                    
                    % normalize
                    beta = norm(w, 'fro');
                    w = w / beta;
                    
                    if ctr_inner == alg.krylovdim
                        dtau = min(dtau, tau - tau0);
                        
                        HH = zeros(ctr_inner + 1);
                        HH(1:ctr_inner, 1:ctr_inner) = H(1:ctr_inner, 1:ctr_inner) * (sgn * dtau);
                        HH(1, ctr_inner + 1) = 1;
                        expH = expm(HH);
                        
                        vareps = abs(beta * H(ctr_inner, ctr_inner) * expH(ctr_inner, ctr_inner + 1));
                        omega = vareps / (dtau * eta);
                        q = ctr_inner / 2;
                        
                        while omega > 1
                            vareps_prev = vareps;
                            dtau_prev = dtau;
                            dtau = dtau * (gamma / omega)^(1 / (q + 1));
                            
                            HH = zeros(ctr_inner + 1);
                            HH(1:ctr_inner, 1:ctr_inner) = H(1:ctr_inner, 1:ctr_inner) * (sgn * dtau);
                            HH(1, ctr_inner+1) = 1;
                            expH = expm(HH);
                            
                            vareps = abs(beta * H(ctr_inner, ctr_inner) * expH(ctr_inner, ctr_inner + 1));
                            omega = vareps / (dtau * eta);
                            q = max(0, log(vareps / vareps_prev) / log(dtau / dtau_prev) - 1);
                        end
                        
                        % take time step
                        totalerr = totalerr + vareps;
                        w = V * expH(1:ctr_inner, ctr_inner);
                        w = w + expH(ctr_inner, end) * w0;
                        w0 = w0 + beta * w;
                        tau0 = tau0 + dtau;
                        
                        % increase time step for next iteration
                        if omega < gamma, dtau = dtau * (gamma / omega) ^(1 / (q + 1)); end
                        
                        if alg.verbosity > Verbosity.iter
                            fprintf('Iter %d: t = %.2e,\terror = %.5e.\n', ...
                                ctr_inner, tau0, totalerr);
                        end
                        
                    elseif H(ctr_inner, ctr_inner) <= (tau - tau0) * eta || alg.nobuild
                        HH = zeros(ctr_inner + 1);
                        HH(1:ctr_inner, 1:ctr_inner) = H(1:ctr_inner, 1:ctr_inner) * (sgn * (tau - tau0));
                        HH(1, ctr_inner+1) = 1;
                        expH = expm(HH);
                        
                        vareps = abs(beta * H(ctr_inner, ctr_inner) * expH(ctr_inner, ctr_inner + 1));
                        omega = vareps / ((tau - tau0) * eta);
                        
                        if omega < 1 % take time step
                            totalerr = totalerr + vareps;
                            w = V(:, 1:ctr_inner) * expH(1:ctr_inner, ctr_inner);
                            w = w + expH(ctr_inner, end) * v0;
                            w0 = w0 + beta * w;
                            tau0 = tau;
                        end
                    end
                    
                    if tau0 >= tau
                        if alg.verbosity >= Verbosity.conv
                            fprintf('Conv %d: error = %.5e.\n', ctr_outer, totalerr);
                        end
                        return
                    end
                end
                
                if ctr_outer == alg.maxiter
                    if alg.verbosity >= Verbosity.conv
                        fprintf('Maxiter %d: error = %.5e\n', ctr_outer, totalerr);
                    end
                    return
                end
                
                % reinitialize
                beta = norm(w, 'fro');
                if beta < alg.tol
                    warning('converged to fixed point.');
                    return
                end
                v = w / beta;
                V = zeros(1, alg.krylovdim, 'like', v0);                 % Krylov subspace basis
                H = zeros(alg.krylovdim, alg.krylovdim, underlyingType(v0)); % Hessenberg matrix
                
            end
        end
    end
end

