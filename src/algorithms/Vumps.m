classdef Vumps < handle
    % Variational fixed point algorithm for uniform matrix product states.
    
    %% Options
    properties
        tol = 1e-10
        miniter = 5
        maxiter = 100
        verbosity = Verbosity.iter
        doplot = false
        which = 'largestabs'
        
        dynamical_tols = true
        tol_min                 = 1e-12
        tol_max                 = 1e-6
        eigs_tolfactor          = 1e-4
        canonical_tolfactor     = 1e-8
        environments_tolfactor  = 1e-4
        
        multiAC = 'parallel'
        dynamical_multiAC = false;
        tol_multiAC = Inf
    end
    
    properties (Access = private)
        alg_eigs = struct('MaxIter', 100, 'KrylovDim', 20)
        alg_canonical = struct('Method', 'polar')
        alg_environments = struct
        
        progressfig
    end
    
    
    %%
    methods
        function alg = Vumps(kwargs)
            arguments
                kwargs.?Vumps
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
            
            if ~isfield('alg_canonical', kwargs)
                alg.alg_canonical.Tol = sqrt(alg.tol_min * alg.tol_max);
                alg.alg_canonical.Verbosity = alg.verbosity - 2;
            end
            
            if ~isfield('alg_environments', kwargs)
                alg.alg_environments.Tol = sqrt(alg.tol_min * alg.tol_max);
                alg.alg_environments.Verbosity = alg.verbosity - 2;
            end
        end
        
        function [mps, lambda, GL, GR] = fixedpoint(alg, mpo, mps)
            
            if period(mpo) ~= period(mps)
                error('vumps:argerror', ...
                    'periodicity of mpo (%d) should be equal to that of the mps (%d)', ...
                    period(mpo), period(mps));
            end
            
            t_total = tic;
            disp_init(alg);
            
            mps = canonicalize(mps);
            [GL, GR] = environments(alg, mpo, mps);
            
            for iter = 1:alg.maxiter
                t_iter = tic;

                AC = updateAC(alg, iter, mpo, mps, GL, GR);
                C  = updateC (alg, iter, mpo, mps, GL, GR);
                mps = updatemps(alg, iter, mps, AC, C);
                
                [GL, GR, lambda] = environments(alg, mpo, mps, GL, GR);
                eta = convergence(alg, mpo, mps, GL, GR);
                
                if iter > alg.miniter && eta < alg.tol
                    disp_conv(alg, iter, lambda, eta, toc(t_total));
                    return
                end
                alg = updatetols(alg, iter, eta);
                plot(alg, iter, mps, eta);
                disp_iter(alg, iter, lambda, eta, toc(t_iter));
            end
            
            disp_maxiter(alg, iter, lambda, eta, toc(t_total));
        end
    end
    
    
    %% Subroutines
    methods
        function AC = updateAC(alg, iter, mpo, mps, GL, GR)
            kwargs = namedargs2cell(alg.alg_eigs);
            if strcmp(alg.multiAC, 'sequential')
                sites = mod1(iter, period(mps));
            else
                sites = 1:period(mps);
            end
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR, sites);
            for i = length(sites):-1:1
                [AC(i), ~] = eigsolve(H_AC{i}, mps.AC(sites(i)), 1, alg.which, ...
                    kwargs{:});
            end
        end
        
        function C = updateC(alg, iter, mpo, mps, GL, GR)
            kwargs = namedargs2cell(alg.alg_eigs);
            if strcmp(alg.multiAC, 'sequential')
                sites = mod1(iter, period(mps));
            else
                sites = 1:period(mps);
            end
            
            H_C = C_hamiltonian(mpo, mps, GL, GR, sites);
            for i = length(sites):-1:1
                [C(i), ~] = eigsolve(H_C{i}, mps.C(sites(i)), 1, alg.which, ...
                    kwargs{:});
            end
        end
        
        function mps = updatemps(alg, iter, mps, AC, C)
            if strcmp(alg.multiAC, 'sequential')
                sites = mod1(iter, period(mps));
            else
                sites = 1:period(mps);
            end
            
            for i = length(AC):-1:1
                [Q_AC, ~] = leftorth(AC(i));
                [Q_C, ~]  = leftorth(C(i), 1, 2);
                mps.AL(sites(i)) = multiplyright(Q_AC, Q_C');
            end
            
            kwargs = namedargs2cell(alg.alg_canonical);
            mps = canonicalize(mps, kwargs{:});
        end
        
        function [GL, GR, lambda] = environments(alg, mpo, mps, GL, GR)
            arguments
                alg
                mpo
                mps
                GL = cell(1, period(mps))
                GR = cell(1, period(mps))
            end
            
            kwargs = namedargs2cell(alg.alg_environments);
            [GL, GR, lambda] = environments(mpo, mps, mps, GL, GR, ...
                kwargs{:});
        end
        
        function eta = convergence(alg, mpo, mps, GL, GR)
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            H_C  = C_hamiltonian(mpo, mps, GL, GR);
            eta = zeros(1, period(mps));
            for w = 1:period(mps)
                AC_ = apply(H_AC{w}, mps.AC(w));
                lambda_AC = dot(AC_, mps.AC(w));
                AC_ = normalize(AC_ ./ lambda_AC);
                
                ww = prev(w, period(mps));
                C_ = apply(H_C{ww}, mps.C(ww));
                lambda_C = dot(C_, mps.C(ww));
                C_ = normalize(C_ ./ lambda_C);
            
                eta(w) = distance(AC_ , ...
                    repartition(multiplyleft(mps.AR(w), C_), rank(AC_)));
            end
            eta = max(eta, [], 'all');
        end
    end
    
    
    %% Option handling
    methods
        function alg = updatetols(alg, iter, eta)
            if alg.dynamical_tols
                alg.alg_eigs.Tol = between(alg.tol_min, eta * alg.eigs_tolfactor, ...
                    alg.tol_max / iter);
                alg.alg_canonical.Tol = between(alg.tol_min, ...
                    eta * alg.canonical_tolfactor, alg.tol_max / iter);
                alg.alg_environments.Tol = between(alg.tol_min, ...
                    eta * alg.environments_tolfactor, alg.tol_max / iter);
                
                if alg.verbosity > Verbosity.iter
                    fprintf('Updated subalgorithm tolerances: (%e,\t%e,\t%e)\n', ...
                        alg.alg_eigs.Tol, alg.alg_canonical.Tol, alg.alg_environments.Tol);
                end
            end
            
            if alg.dynamical_multiAC
                if eta < alg.tol_multiAC
                    alg.multiAC = 'sequential';
                else
                    alg.multiAC = 'parallel';
                end
            end
        end
    end
    
    %% Display
    methods (Access = private)
        function plot(alg, iter, mps, eta)
            if ~alg.doplot, return; end
            persistent axhistory axspectrum
            
            D = depth(mps);
            W = period(mps);
            
            if isempty(alg.progressfig) || ~isvalid(alg.progressfig) || iter == 1
                alg.progressfig = figure('Name', 'Vumps');
                axhistory = subplot(D + 1, W, 1:W);
                axspectrum = gobjects(D, W);
                for d = 1:D, for w = 1:W
                        axspectrum(d, w) = subplot(D + 1, W, w + d * W);
                    end, end
                linkaxes(axspectrum, 'y');
            end
           
            
            if isempty(axhistory.Children)
                semilogy(axhistory, iter, eta, '.', 'Color', colors(1), ...
                    'DisplayName', 'Errors', 'MarkerSize', 10);
                hold on
                ylim(axhistory, [alg.tol / 5, max([1 eta])]);
                yline(axhistory, alg.tol, '-', 'Convergence', ...
                    'LineWidth', 2);
                hold off
            else
                axhistory.Children(end).XData(end+1) = iter;
                axhistory.Children(end).YData(end+1) = eta;
            end
            
            plot_entanglementspectrum(mps, 1:period(mps), axspectrum);
            drawnow
        end
        
        function disp_init(alg)
            if alg.verbosity < Verbosity.conv, return; end
            fprintf('---- VUMPS ----\n');
        end
        
        function disp_iter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.iter, return; end
            
            s = settings;
            if abs(imag(lambda)) < eps(lambda)^(3/4) * abs(lambda)
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('Vumps %2d:\tE = %-0.4f\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t));
                    otherwise
                        fprintf('Vumps %4d:\tE = %-0.15f\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t, 's'));
                        
                end
            else
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('Vumps %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                    otherwise
                        fprintf('Vumps %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                        
                end
            end
        end
        
        function disp_conv(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.conv, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('Vumps converged %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                otherwise
                    fprintf('Vumps converged %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
        
        function disp_maxiter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.warn, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('Vumps max iterations %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                otherwise
                    fprintf('Vumps max iterations %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
    end
end

