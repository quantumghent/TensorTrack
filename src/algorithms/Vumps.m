classdef Vumps
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Options
    properties
        tol = 1e-10
        miniter = 5
        maxiter = 100
        verbosity = Verbosity.iter
        which = 'largestabs'
        
        
        KrylovDim = 20
        
        dynamical_tols = true
        tol_min                 = 1e-12
        tol_max                 = 1e-6
        eigs_tolfactor          = 1e-4
        canonical_tolfactor     = 1e-8
        environments_tolfactor  = 1e-4
    end
    
    properties (Access = private)
        alg_eigs = struct('MaxIter', 100, 'KrylovDim', 20)
        alg_canonical = struct('Method', 'polar')
        alg_environments = struct
    end
    
    
    %%
    methods
        function v = Vumps(kwargs)
            arguments
                kwargs.?Vumps
            end
            
            fields = fieldnames(kwargs);
            if ~isempty(fields)
                for field = fields.'
                    v.(field{1}) = kwargs.(field{1});
                end
            end
            
            if ~isfield('alg_eigs', kwargs)
                v.alg_eigs.Tol = sqrt(v.tol_min * v.tol_max);
                v.alg_eigs.Verbosity = v.verbosity - 2;
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
        
        function [mps, lambda, GL, GR] = fixedpoint(alg, mpo, mps)
            t_total = tic;
            disp_init(alg);
            
            mps = canonicalize(mps);
            [GL, GR] = environments(mpo, mps);
            
            for iter = 1:alg.maxiter
                t_iter = tic;
                AC = updateAC(alg, mpo, mps, GL, GR);
                C  = updateC (alg, mpo, mps, GL, GR);
                mps = updatemps(alg, AC, C);
                
                [GL, GR, lambda] = environments(alg, mpo, mps, GL, GR);
                eta = convergence(alg, mpo, mps, GL, GR);
                
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
        function AC = updateAC(alg, mpo, mps, GL, GR)
            kwargs = namedargs2cell(alg.alg_eigs);
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            for i = period(mps):-1:1
                [AC(i), ~] = eigsolve(H_AC{i}, mps.AC(i), 1, alg.which, kwargs{:});
            end
        end
        
        function C = updateC(alg, mpo, mps, GL, GR)
            kwargs = namedargs2cell(alg.alg_eigs);
            H_C = C_hamiltonian(mpo, mps, GL, GR);
            for i = period(mps):-1:1
                [C(i), ~] = eigsolve(H_C{i}, mps.C(i), 1, alg.which, kwargs{:});
            end
        end
        
        function mps = updatemps(alg, AC, C)
            for i = length(AC):-1:1
                [~, Q_AC] = rightorth(AC(i));
                [~, Q_C]  = rightorth(C(i), 1, 2);
                
                AR(i) = contract(conj(Q_C), [1 -1], Q_AC, [1 -2 -3], 'Rank', [1 2]);
            end
            mps = UniformMps([], AR, C, []);
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
            
                C_ = apply(H_C{w}, mps.C(w));
                lambda_C = dot(C_, mps.C(w));
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
            if ~alg.dynamical_tols, return; end
            
            alg.alg_eigs.Tol = between(alg.tol_min, eta * alg.eigs_tolfactor, ...
                alg.tol_max / iter);
            alg.alg_canonical.Tol = between(alg.tol_min, eta * alg.canonical_tolfactor, ...
                alg.tol_max / iter);
            alg.alg_environments.Tol = between(alg.tol_min, eta * alg.environments_tolfactor, ...
                alg.tol_max / iter);
            
            if alg.verbosity > Verbosity.iter
                fprintf('Updated subalgorithm tolerances: (%e,\t%e,\t%e)\n', ...
                    alg.alg_eigs.Tol, alg.alg_canonical.Tol, alg.alg_environments.Tol);
            end
        end
    end
    
    %% Display
    methods (Access = private)
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
                        fprintf('Iter %2d:\tE = %-0.4f\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t));
                    otherwise
                        fprintf('Iter %4d:\tE = %-0.15f\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t, 's'));
                        
                end
            else
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('Iter %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                    otherwise
                        fprintf('Iter %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                        
                end
            end
        end
        
        function disp_conv(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.conv, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('Conv %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                otherwise
                    fprintf('Conv %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
        
        function disp_maxiter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.warn, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('MaxIter %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                otherwise
                    fprintf('MaxIter %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                        iter, real(lambda), imag(lambda), eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
    end
end

