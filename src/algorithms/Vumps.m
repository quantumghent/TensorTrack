classdef Vumps < handle
    % `Variational fixed point algorithm for uniform matrix product states <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.045145>`_.
    %
    % Properties
    % ----------
    % tol : :class:`double`
    %   tolerance for convergence criterion, defaults to :code:`1e-10`.
    %
    % miniter : :class:`int`
    %   minimum number of iteration, defaults to :code:`5`.
    %
    % maxiter : :class:`int`
    %   maximum number of iteration, defaults to code:`100`.
    %
    % verbosity : :class:`.Verbosity`
    %   verbosity level of the algorithm, defaults to :code:`Verbosity.iter`.
    %
    % doplot : :class:`logical`
    %   plot progress, defaults to :code:`false`.
    %
    % which : :class:`char`
    %   eigenvalue selector (passed as the :code:`sigma` argument to :func:`.eigsolve`),
    %   defaults to :code:`'largestabs'`.
    %
    % dynamical_tols : :class:`logical`
    %   indicate whether or not to use a dynamical tolerance scaling for the algorithm's
    %   subroutines based on the current error measure, defaults to :code:`true`
    %
    % tol_min : :class:`double`
    %   smallest allowed convergence tolerance for soubroutines, defaults to :code:`1e-12`.
    %
    % tol_max : :class:`double`
    %   highest allowed convergence tolerance for soubroutines, defaults to :code:`1e-10`.
    %
    % eigs_tolfactor : :class:`double`
    %   relative scaling factor for determining the convergence tolerance of the local
    %   update solver subroutine based on the current error measure, defaults to
    %   :code:`1e-6`.
    %
    % canonical_tolfactor : :class:`double`
    %   relative scaling factor for determining the convergence tolerance of the
    %   canonicalization subroutine based on the current error measure, defaults to
    %   :code:`1e-8`.
    %
    % environments_tolfactor : :class:`double`
    %   relative scaling factor for determining the convergence tolerance of the environment
    %   solver subroutine based on the current error measure, defaults to :code:`1e-6`.
    %
    % multiAC : :class:`char`
    %   execution style for the local `AC` updates for a multi-site unit cell, options are:
    %
    %   - :code:`'parallel'`: (default) update all `AC` tensors simultaneously.
    %   - :code:`'sequential'`: update one `AC` tensor at a time, sweeping through the unit
    %     cell.
    %
    % dynamical_multiAC : :class:`logical`
    %   automatically switch from :code:`'sequential'` to :code:`'parallel'` if the error
    %   measure becomes small enough, defaults to :code:`false`.
    %
    % tol_multiAC : :class:`char`
    %   tolerance for automatically switching from :code:`'sequential'` to
    %   :code:`'parallel'` if the error measure falls below this value, defaults to
    %   :code:`Inf`.
    %
    % alg_eigs : :class:`.KrylovSchur` or :class:`.Arnoldi`
    %   algorithm used for the eigsolver subroutines, defaults to
    %   :code:`Arnoldi('MaxIter', 100, 'KrylovDim', 20)`.

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
        tol_max                 = 1e-10
        eigs_tolfactor          = 1e-6
        canonical_tolfactor     = 1e-8
        environments_tolfactor  = 1e-6
        
        multiAC = 'parallel'
        dynamical_multiAC = false;
        tol_multiAC = Inf

        doSave = false
        saveIterations = 1
        saveMethod = 'full'
        name = 'VUMPS'
        
        alg_eigs = Arnoldi('MaxIter', 100, 'KrylovDim', 20)
    end
    
    properties (Access = private)
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
                alg.alg_eigs.tol = sqrt(alg.tol_min * alg.tol_max);
                alg.alg_eigs.verbosity = alg.verbosity - 2;
            end
            
            if ~isfield('alg_canonical', kwargs)
                alg.alg_canonical.tol = sqrt(alg.tol_min * alg.tol_max);
                alg.alg_canonical.verbosity = alg.verbosity - 2;
            end
            
            if ~isfield('alg_environments', kwargs)
                alg.alg_environments.tol = sqrt(alg.tol_min * alg.tol_max);
                alg.alg_environments.verbosity = alg.verbosity - 2;
            end
            
            if isfield('verbosity', kwargs)
                alg.alg_eigs.verbosity = alg.verbosity - 2;
                alg.alg_environments.verbosity = alg.verbosity - 2;
            end
        end
        
        function [mps, lambda, GL, GR, eta] = fixedpoint(alg, mpo, mps)
            % Find the fixed point MPS of an infinite MPO, given an initial guess.
            %
            % Usage
            % -----
            % :code:`[mps, lambda, GL, GR, eta] = fixedpoint(alg, mpo, mps)`
            %
            % Arguments
            % ---------
            % alg : :class:`.Vumps`
            %   VUMPS algorithm.
            %
            % mpo : :class:`.InfMpo`
            %   matrix product operator.
            %
            % mps : :class:`.UniformMps`
            %   initial guess for MPS fixed point.
            %
            % Returns
            % -------
            % mps : :class:`.UniformMps`
            %   MPS fixed point.
            %
            % lambda : :class:`double`
            %   eigenvalue.
            %
            % GL : :class:`cell` of :class:`.MpsTensor`
            %   left environment tensors.
            %
            % GR : :class:`cell` of :class:`.MpsTensor`
            %   right environment tensors.
            %
            % eta : :class:`double`
            %   final error measure at convergence.
            
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

                if alg.doSave && mod(iter, alg.saveIterations) == 0
                    save_iteration(alg, mps, lambda, iter, eta, toc(t_total));
                end
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
                if alg.verbosity >= Verbosity.detail
                    fprintf('\nAC{%d} eigenvalue solver:\n------------------------\n', sites(i));
                end
                [AC{1, i}, ~] = eigsolve(alg.alg_eigs, @(x) H_AC{i}.apply(x), AC{1, i}, ...
                    1, alg.which);
                for d = 2:depth(mpo)
                    AC{d, i} = H_AC{i}(d).apply(AC{d-1, i});
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
                if alg.verbosity >= Verbosity.detail
                    fprintf('\nC{%d} eigenvalue solver:\n-----------------------\n', sites(i));
                end
                [C{1, i}, ~] = eigsolve(alg.alg_eigs, @(x) H_C{i}.apply(x), C{1, i}, 1, ...
                    alg.which);
                for d = 2:depth(mpo)
                    C{d, i} = H_C{i}(d).apply(C{d-1, i});
                end
            end
        end
        
        function mps = updatemps(alg, iter, mps, AC, C)
            if alg.verbosity >= Verbosity.detail
                fprintf('\nCanonicalize:\n------------------\n');
            end
            
            if strcmp(alg.multiAC, 'sequential')
                sites = mod1(iter, period(mps));
            else
                sites = 1:period(mps);
            end
            
            for d = size(AC, 1):-1:1
                for i = size(AC, 2):-1:1
                    [Q_AC, ~] = leftorth(AC{d, i}, 'polar');
                    [Q_C, ~]  = leftorth(C{d, i}, 1, 2, 'polar');
                    mps(d).AL{sites(i)} = multiplyright(Q_AC, Q_C');
                end
            end
            
            kwargs = namedargs2cell(alg.alg_canonical);
            mps = canonicalize(mps, kwargs{:});
        end
        
        function [GL, GR, lambda] = environments(alg, mpo, mps, GL, GR)
            arguments
                alg
                mpo
                mps
                GL = cell(depth(mpo), period(mps))
                GR = cell(depth(mpo), period(mps))
            end
            if alg.verbosity >= Verbosity.detail
                fprintf('\nEnvironments:\n------------------\n');
            end
            kwargs = namedargs2cell(alg.alg_environments);
            D = depth(mpo);
            lambda = zeros(D, 1);
            for d = 1:D
                dd = next(d, D);
                [GL(d, :), GR(d, :), lambda(d)] = environments(mpo.slice(d), ...
                    mps(d), mps(dd), GL(d, :), GR(d, :), kwargs{:});
            end
            lambda = prod(lambda);
        end
        
        function eta = convergence(alg, mpo, mps, GL, GR)
            % TODO: also implement galerkin error
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            H_C  = C_hamiltonian(mpo, mps, GL, GR);
            eta = zeros(depth(mpo), period(mps));
            for d = 1:depth(mpo)
                dd = next(d, depth(mpo));
                for w = 1:period(mps)
                    AC_ = apply(H_AC{w}(d), mps(d).AC{w});
                    lambda_AC = dot(AC_, mps(dd).AC{w});
                    AC_ = normalize(AC_ ./ lambda_AC);

                    ww = prev(w, period(mps));
                    C_ = apply(H_C{ww}(d), mps(d).C{ww});
                    lambda_C = dot(C_, mps(dd).C{ww});
                    C_ = normalize(C_ ./ lambda_C);

                    eta(dd, w) = distance(AC_ , ...
                        repartition(multiplyleft(mps(dd).AR{w}, C_), rank(AC_)));
                end
            end
            eta = max(eta, [], 'all');
        end
    end
    
    
    %% Option handling
    methods
        function alg = updatetols(alg, iter, eta)
            if alg.dynamical_tols
                alg.alg_eigs.tol = between(alg.tol_min, eta * alg.eigs_tolfactor / iter, ...
                    alg.tol_max);
                alg.alg_canonical.Tol = between(alg.tol_min, ...
                    eta * alg.canonical_tolfactor / iter, alg.tol_max);
                alg.alg_environments.Tol = between(alg.tol_min, ...
                    eta * alg.environments_tolfactor / iter, alg.tol_max);
                
                if alg.verbosity > Verbosity.iter
                    fprintf('Updated subalgorithm tolerances: (%e,\t%e,\t%e)\n', ...
                        alg.alg_eigs.tol, alg.alg_canonical.Tol, alg.alg_environments.Tol);
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
            
            plot_entanglementspectrum(mps, 1:D, 1:W, axspectrum);
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
                
        function save_iteration(alg, mps, lambda, iter, eta, t)
            fileName = alg.name;

            fileData = struct;
            fileData.mps            = mps;
            fileData.lambda         = lambda;
            fileData.iteration = iter;
            fileData.eta = eta;
            fileData.time = t;
            % save
            
            %if exist(fileName,'file')
            %    old_file=load(fileName);
            %    fileName_temp=[fileName(1:end-4),'_temp.mat'];
            %    save(fileName_temp, '-struct', 'old_file', '-v7.3');
            %    saved_temp=1;
            %else
            %    saved_temp=0;
            %end
            
            save(fileName, '-struct', 'fileData', '-v7.3');
            
            %if saved_temp
            %    delete(fileName_temp);
            %end
        end
    end
end

