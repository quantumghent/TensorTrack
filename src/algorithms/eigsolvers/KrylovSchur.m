classdef KrylovSchur
    % KrylovSchur wrapper for Matlab implementation of eigs
    
    properties
        tol         = 1e-10             % convergence tolerance
        maxiter     = 100               % maximum iterations
        krylovdim   = 20                % Krylov subspace dimension
        verbosity   = Verbosity.warn    % display information
    end
    
    methods
        function alg = KrylovSchur(kwargs)
            arguments
                kwargs.?KrylovSchur
            end
            fields = fieldnames(kwargs);
            if ~isempty(fields)
                for field = fields.'
                    alg.(field{1}) = kwargs.(field{1});
                end
            end
        end
        
        function varargout = eigsolve(alg, A, v0, howmany, sigma, kwargs)
            arguments
                alg
                A
                v0
                howmany = 1
                sigma = 'lm'
                kwargs.IsSymmetric logical = false
            end
            
            assert(isnumeric(sigma) || ismember(sigma, {'largestabs', 'smallestabs', ...
                'largestreal', 'smallestreal', 'bothendsreal', ...
                'largestimag', 'smallestimag', 'bothendsimag'}), ...
                'tensors:ArgumentError', 'Invalid choice of eigenvalue selector.');
            nargoutchk(0, 3);
            
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
            
            sz = size(v0_vec);
            
            if sz(1) < howmany
                warning('eigsolve:size', 'requested %d out of %d eigenvalues.', ...
                    howmany, sz(1)); 
                howmany = sz(1);
            end
            
            if sz(1) < alg.krylovdim
                warning('eigsolve:size', 'requested %d out of %d eigenvalues.', ...
                    howmany, sz(1)); 
                alg.krylovdim = sz(1);
            end
            
            % call builtin eigs
            if howmany > sz(1) - 2
                % annoying bug with eigs and subspace dimension specification
                [V, D, flag] = eigs(A_fun, sz(1), howmany, sigma, ...
                    'Tolerance', alg.tol, 'MaxIterations', alg.maxiter, ...
                    'IsFunctionSymmetric', ...
                    kwargs.IsSymmetric, 'StartVector', v0_vec, ...
                    'Display', alg.verbosity >= Verbosity.iter);
            else
                [V, D, flag] = eigs(A_fun, sz(1), howmany, sigma, ...
                    'Tolerance', alg.tol, 'MaxIterations', alg.maxiter, ...
                    'SubspaceDimension', alg.krylovdim, 'IsFunctionSymmetric', ...
                    kwargs.IsSymmetric, 'StartVector', v0_vec, ...
                    'Display', alg.verbosity >= Verbosity.iter);
            end
            
            % process results
            if nargout <= 1
                varargout = {D};
            else
                if ~isnumeric(v0)
                    for i = howmany:-1:1
                        varargout{1}(:, i) = devectorize(V(:, i), v0);
                    end
                end
                varargout{2} = D;
                if nargout == 3
                    varargout{3} = flag;
                end
            end
            
            % display
            if flag
                warning('eigsolve did not converge.');
            elseif ~flag && alg.verbosity > Verbosity.warn
                fprintf('eigsolve converged.\n');
            end
        end
    end
end


