classdef Vomps
    % Fixed point algorithm for maximizing overlap.
    
    %% Options
    properties
        tol = 1e-5
        miniter = 2
        maxiter = 10
        verbosity = Verbosity.iter
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
        alg_canonical = struct('Method', 'polar')
        alg_environments = struct
    end
    
    methods
        function v = Vomps(kwargs)
            arguments
                kwargs.?Vumps
            end
            
            fields = fieldnames(kwargs);
            if ~isempty(fields)
                for field = fields.'
                    v.(field{1}) = kwargs.(field{1});
                end
            end
            
            if ~isfield('alg_canonical', kwargs)
                v.alg_canonical.Tol = sqrt(v.tol_min * v.tol_max);
                v.alg_canonical.Verbosity = v.verbosity - 2;
            end
            
            if ~isfield('alg_environments', kwargs)
                v.alg_environments.Tol = sqrt(v.tol_min * v.tol_max);
                v.alg_environments.Verbosity = v.verbosity - 2;
            end
        end
        
        function [mps2, GL, GR] = approximate(alg, mpo, mps1, mps2)
            if period(mpo) ~= period(mps1) || period(mpo) ~= period(mps2)
                error('vumps:argerror', ...
                    'periodicitys should match: mpo (%d), mps1 (%d), mps2(%d)', ...
                    period(mpo), period(mps1), period(mps2));
            end
            
            t_total = tic;
            disp_init(alg);
            
            mps2 = canonicalize(mps2);
            [GL, GR] = environments(alg, mpo, mps1, mps2);
            
            for iter = 1:alg.maxiter
                t_iter = tic;

                AC = updateAC(alg, iter, mpo, mps1, GL, GR);
                C  = updateC (alg, iter, mpo, mps1, GL, GR);
                mps2 = updatemps(alg, iter, mps2, AC, C);
                
                [GL, GR, lambda] = environments(alg, mpo, mps1, mps2, GL, GR);
                eta = convergence(alg, mpo, mps1, mps2, GL, GR);
                
                if iter > alg.miniter && eta < alg.tol
                    disp_conv(alg, iter, lambda, eta, toc(t_total));
                    return
                end
                alg = updatetols(alg, iter, eta);
                disp_iter(alg, iter, lambda, eta, toc(t_iter));
            end
            
            disp_maxiter(alg, iter, lambda, eta, toc(t_total));
        end
    end
    
    %% Subroutines
    methods
        function AC = updateAC(alg, iter, mpo, mps, GL, GR)
            if strcmp(alg.multiAC, 'sequential')
                sites = mod1(iter, period(mps));
            else
                sites = 1:period(mps);
            end
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR, sites);
            ACs = arrayfun(@(x) x.AC(sites), mps, 'UniformOutput', false);
            AC = vertcat(ACs{:});
            for i = length(sites):-1:1
                for d = 1:depth(mpo)
                    AC{d, i} = H_AC{i}(d).apply(AC{prev(d, depth(mpo)), i});
                end
            end
        end
        
        function C = updateC(alg, iter, mpo, mps, GL, GR)
            if strcmp(alg.multiAC, 'sequential')
                sites = mod1(iter, period(mps));
            else
                sites = 1:period(mps);
            end
            
            H_C = C_hamiltonian(mpo, mps, GL, GR, sites);
            Cs = arrayfun(@(x) x.C(sites), mps, 'UniformOutput', false);
            C = vertcat(Cs{:});
            for i = length(sites):-1:1
                for d = 1:depth(mpo)
                    C{d, i} = H_C{i}(d).apply(C{prev(d, depth(mpo)), i});
                end
            end
        end
        
        function mps = updatemps(alg, iter, mps, AC, C)
            if strcmp(alg.multiAC, 'sequential')
                sites = mod1(iter, period(mps));
            else
                sites = 1:period(mps);
            end
            
            for d = size(AC, 1):-1:1
                for i = length(AC):-1:1
                    [Q_AC, ~] = leftorth(AC{d, i}, 'polar');
                    [Q_C, ~]  = leftorth(C{d, i}, 1, 2, 'polar');
                    mps(d).AL{sites(i)} = multiplyright(Q_AC, Q_C');
                end
            end
            
            kwargs = namedargs2cell(alg.alg_canonical);
            mps = canonicalize(mps, kwargs{:});
        end
        
        function [GL, GR, lambda] = environments(alg, mpo, mps1, mps2, GL, GR)
            arguments
                alg
                mpo
                mps1
                mps2
                GL = cell(1, period(mps1))
                GR = cell(1, period(mps1))
            end
            
            kwargs = namedargs2cell(alg.alg_environments);
            [GL, GR, lambda] = environments(mpo, mps1, mps2, GL, GR, ...
                kwargs{:});
        end
        
        function eta = convergence(alg, mpo, mps1, mps2, GL, GR)
            H_AC = AC_hamiltonian(mpo, mps1, GL, GR);
            H_C  = C_hamiltonian(mpo, mps1, GL, GR);
            eta = zeros(1, period(mps1));
            for w = 1:period(mps1)
                AC_ = apply(H_AC{w}, mps1.AC{w});
                lambda_AC = dot(AC_, mps2.AC{w});
                AC_ = normalize(AC_ ./ lambda_AC);
                
                ww = prev(w, period(mps1));
                C_ = apply(H_C{ww}, mps1.C{ww});
                lambda_C = dot(C_, mps2.C{ww});
                C_ = normalize(C_ ./ lambda_C);
            
                eta(w) = distance(AC_ , ...
                    repartition(multiplyleft(mps2.AR{w}, C_), rank(AC_)));
            end
            eta = max(eta, [], 'all');
        end
    end
    
    %% Option handling
    methods
        function alg = updatetols(alg, iter, eta)
            if alg.dynamical_tols
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
        function disp_init(alg)
            if alg.verbosity < Verbosity.conv, return; end
            fprintf('---- VOMPS ----\n');
        end
        
        function disp_iter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.iter, return; end
            
            s = settings;
            if abs(imag(lambda)) < eps(lambda)^(3/4) * abs(lambda)
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('Vomps %2d:\tE = %-0.4f\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t));
                    otherwise
                        fprintf('Vomps %4d:\tE = %-0.15f\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t, 's'));
                        
                end
            else
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('Vomps %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                    otherwise
                        fprintf('Vomps %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                        
                end
            end
        end
        
        function disp_conv(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.conv, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('Vomps converged %2d:\terror = %0.1e\t(%s)\n', ...
                        iter, eta, time2str(t));
                otherwise
                    fprintf('Vomps converged %4d:\terror = %0.4e\t(%s)\n', ...
                        iter, eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
        
        function disp_maxiter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.warn, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('Vomps max iterations %2d:\terror = %0.1e\t(%s)\n', ...
                        iter, eta, time2str(t));
                otherwise
                    fprintf('Vomps max iterations %4d:\terror = %0.4e\t(%s)\n', ...
                        iter, eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
    end
end

