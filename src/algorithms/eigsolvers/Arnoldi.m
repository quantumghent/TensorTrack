classdef Arnoldi
    % Arnoldi Krylov algorithm for linear algebra problems.
    
    properties
        tol         = 1e-10             % convergence tolerance
        maxiter     = 100               % maximum iterations
        krylovdim   = 20                % Krylov subspace dimension
        deflatedim  = 3                 % number of Krylov vectors to keep when deflating
        reorth      = 20                % reorthogonalize basis if larger than this number
        nobuild     = 1                 % frequency of convergence check when building
        verbosity   = Verbosity.warn    % display information
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
        
        function [V, D, flag] = eigsolve(alg, A, v, howmany, sigma)
            arguments
                alg
                A
                v
                howmany = 1
                sigma = 'lm'
            end
            
            if isempty(alg.nobuild), alg.nobuild = ceil(alg.krylovdim / 10); end
            if isempty(alg.deflatedim), alg.deflatedim = max(round(3/5 * alg.krylovdim), howmany); end
            
            v = v / norm(v, 'fro');
            
            % preallocation
            if isnumeric(v) && isvector(v)
                V = zeros(length(v), alg.krylovdim, 'like', v);
            else
                V = zeros(1, alg.krylovdim, 'like', v);                 % Krylov subspace basis
            end
            H = zeros(alg.krylovdim, alg.krylovdim, underlyingType(v)); % Hessenberg matrix
            
            ctr_outer = 0;
            ctr_inner = 0;
            flag = 0;
            
            while ctr_outer < alg.maxiter
                ctr_outer = ctr_outer + 1;
                
                while ctr_inner < alg.krylovdim  % build Krylov subspace
                    ctr_inner = ctr_inner + 1;
                    
                    V(:, ctr_inner) = v;
                    v = A(v);
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
                    v = v / beta;
                    
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
                                    fprintf('Conv %2d (%2d/%2d): error = %.5e.\n', ctr_outer, ...
                                        ctr_inner, alg.krylovdim, conv);
                                end
                                return
                            end
                            
                            if alg.verbosity >= Verbosity.detail
                                fprintf('Iter %2d (%2d/%2d):\tlambda = %.5e + %.5ei;\terror = %.5e\n', ...
                                    ctr_outer, ctr_inner, alg.krylovdim, ...
                                    real(lambda(1)), imag(lambda(1)), conv);
                            end
                        end
                    end
                    
                    H(ctr_inner + 1, ctr_inner) = beta;
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
                            fprintf('Found invariant subspace (error = %.5e).\n', conv);
                            flag = 1;
                        else
                            fprintf('Reached maxiter without convergence.\n');
                            flag = 2;
                        end
                    end
                    return
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
                        fprintf('Conv %2d: error = %.5e.\n', ctr_outer, conv);
                    end
                    return
                end
                
                if alg.verbosity >= Verbosity.iter
                    fprintf('Iter %2d:\tlambda = %.5e + %.5ei;\terror = %.5e\n', ...
                        ctr_outer, real(lambda(1)), imag(lambda(1)), conv);
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

