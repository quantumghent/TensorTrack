classdef IDmrg
    % Infinite Density Matrix Renormalization Group algorithm
    
    properties
        tol = 1e-10
        miniter = 5
        maxiter = 100
        verbosity = Verbosity.iter
        which = 'largestabs'
        
        dynamical_tols = true
        tol_min                 = 1e-12
        tol_max                 = 1e-6
        eigs_tolfactor          = 1e-4
    end
    
    properties (Access = private)
        alg_eigs = struct('MaxIter', 100, 'KrylovDim', 20)
    end
    
    methods
        function alg = IDmrg(kwargs)
            arguments
                kwargs.?IDmrg
            end
            
            fields = fieldnames(kwargs);
            if ~isempty(fields)
                for field = fields.'
                    alg.(field{1}) = kwargs.(field{1});
                end
            end
            
            if ~isfield('alg_eigs', kwargs)
                alg.alg_eigs.Tol = sqrt(alg.tol_min * alg.tol_max);
                alg.alg_eigs.Verbosity = alg.verbosity - 2;
            end
        end
        
        function [mps, lambda, GL, GR] = fixedpoint(alg, mpo, mps)
            if period(mpo) ~= period(mps)
                error('idmrg:argerror', ...
                    'periodicity of mpo (%d) should be equal to that of the mps (%d)', ...
                    period(mpo), period(mps));
            end
            
            t_total = tic;
            disp_init(alg);
            
            [GL, GR] = environments(mpo, mps, mps);
            
            for iter = 1:alg.maxiter
                t_iter = tic;
                
                C_ = mps.C(end);
                kwargs = {};
                lambdas = zeros(1, period(mps));
                for pos = 1:period(mps)
                    H = AC_hamiltonian(mpo, mps, GL, GR, pos);
                    [mps.AC(pos).var, lambdas(pos)] = ...
                        eigsolve(H{1}, mps.AC(pos).var, 1, alg.which, kwargs{:});
                    [mps.AL(pos), mps.C(pos)] = leftorth(mps.AC(pos));
                    
                    T = transfermatrix(mpo, mps, mps, pos, 'Type', 'LL');
                    GL{next(pos, period(mps))} = apply(T, GL{pos}) / lambdas(pos);
                end
                
                for pos = period(mps):-1:1
                    H = AC_hamiltonian(mpo, mps, GL, GR, pos);
                    [mps.AC(pos).var, lambdas(pos)] = ...
                        eigsolve(H{1}, mps.AC(pos).var, 1, alg.which, kwargs{:});
                    [mps.C(prev(pos, period(mps))), mps.AR(pos)] = rightorth(mps.AC(pos));
                    
                    T = transfermatrix(mpo, mps, mps, pos, 'Type', 'RR').';
                    GR{pos} = apply(T, GR{next(pos, period(mps))}) / lambdas(pos);
                end

                eta = distance(C_, mps.C(end));
                lambda = prod(lambdas);
                
                if iter > alg.miniter && eta < alg.tol
                    mps = canonicalize(mps);
                    disp_conv(alg, iter, lambda, eta, toc(t_total));
                    return
                end
                
                alg = updatetols(alg, iter, eta);
                disp_iter(alg, iter, lambda, eta, toc(t_iter));
            end
            
            mps = canonicalize(mps);
            disp_maxiter(alg, iter, lambda, eta, toc(t_total));
        end
    end
    
    
    methods
        function [GL, GR, lambda] = environments(alg, mpo, mps, GL, GR)
            arguments
                alg
                mpo
                mps
                GL = cell(1, period(mps))
                GR = cell(1, period(mps))
            end
            
            [GL, GR, lambda] = environments(mpo, mps, mps, GL, GR);
        end
        
        function alg = updatetols(alg, iter, eta)
            if alg.dynamical_tols
                alg.alg_eigs.Tol = between(alg.tol_min, eta * alg.eigs_tolfactor, ...
                    alg.tol_max / iter);
                
                if alg.verbosity > Verbosity.iter
                    fprintf('Updated subalgorithm tolerances: (%e,\t%e,\t%e)\n', ...
                        alg.alg_eigs.Tol, alg.alg_canonical.Tol, alg.alg_environments.Tol);
                end
            end
        end
    end
    
    %% Display
    methods (Access = private)
        function disp_init(alg)
            if alg.verbosity < Verbosity.conv, return; end
            fprintf('---- IDmrg ----\n');
        end
        
        function disp_iter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.iter, return; end
            
            s = settings;
            if abs(imag(lambda)) < eps(lambda)^(3/4) * abs(lambda)
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('IDmrg %2d:\tE = %-0.4f\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t));
                    otherwise
                        fprintf('IDmrg %4d:\tE = %-0.15f\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t, 's'));
                        
                end
            else
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('IDmrg %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                    otherwise
                        fprintf('IDmrg %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                        
                end
            end
        end
        
        function disp_conv(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.conv, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('IDmrg converged %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                otherwise
                    fprintf('IDmrg converged %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
        
        function disp_maxiter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.warn, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('IDmrg max iterations %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                otherwise
                    fprintf('IDmrg max iterations %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
    end
end

