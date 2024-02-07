classdef Dmrg
    % `Density Matrix Renormalisation Group algorithm <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.2863>`_ for marix product states.
    %
    % Properties
    % ----------
    % tol : :class:`double`
    %   tolerance for convergence criterion, defaults to :code:`1e-10`.
    %
    % miniter : :class:`int`
    %   minimum number of iteration, defaults to 5.
    %
    % maxiter : :class:`int`
    %   maximum number of iteration, defaults to 100.
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
    %   subroutines based on the current error measure, defaults to :code:`false`
    %
    % tol_min : :class:`double`
    %   smallest allowed convergence tolerance for soubroutines, defaults to :code:`1e-12`.
    %
    % tol_max : :class:`double`
    %   highest allowed convergence tolerance for soubroutines, defaults to :code:`1e-10`.
    %
    % eigs_tolfactor : :class:`double`
    %   relative scaling factor for determining the convergence tolerance of the eigensolver
    %   subroutine based on the current error measure, defaults to :code:`1e-6`
    %
    % sweepstyle : :class:`char`
    %   sweep style indicating how to sweep through the MPS at each iteration, options are:
    %
    %   - :code:`'f2f'`: (default) front-to-front, sweeping from site 1 to the end and back.
    %   - :code:`'b2b'`: back-to-back, sweeping from site N to the start and back.
    %   - :code:`'m2m'`: mid-to-mid, sweeping from the middle site to both ends and back.
    %
    % alg_eigs : :class:`.KrylovSchur` or :class:`.Arnoldi`
    %   algorithm used for the local eigsolver updates, defaults to
    %   :code:`KrylovSchur('MaxIter', 100, 'KrylovDim', 20)`.
    
    %% Options
    properties
        tol = 1e-10
        miniter = 5
        maxiter = 100
        verbosity = Verbosity.iter
        doplot = false
        which = 'largestabs'
        
        dynamical_tols = false
        tol_min                 = 1e-12
        tol_max                 = 1e-10
        eigs_tolfactor          = 1e-6
        
        sweepstyle {mustBeMember(sweepstyle, {'f2f', 'b2b', 'm2m'})} = 'f2f'
        
        doSave = false
        saveIterations = 1
        saveMethod = 'full'
        name = 'DMRG'
        
        alg_eigs = KrylovSchur('MaxIter', 100, 'KrylovDim', 20)
    end
    
    properties (Access = private)
        progressfig
    end
    
    methods
        function alg = Dmrg(kwargs)
            arguments
                kwargs.?Dmrg
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
        end
        
        function [mps, envs, eta] = fixedpoint(alg, mpo, mps, envs)
            % Find the fixed point MPS of a finite MPO, given an initial guess.
            %
            % Usage
            % -----
            % :code:`[mps, envs, eta] = fixedpoint(alg, mpo, mps, envs)`
            %
            % Arguments
            % ---------
            % alg : :class:`.Dmrg`
            %   DMRG algorithm.
            %
            % mpo : :class:`.FiniteMpo`
            %   matrix product operator.
            %
            % mps : :class:`.FiniteMps`
            %   initial guess for MPS fixed point.
            %
            % envs : :class:`.FiniteEnvironment`
            %   initial guess for the environments.
            %
            % Returns
            % -------
            % mps : :class:`.FiniteMps`
            %   MPS fixed point.
            %
            % envs : :class:`.FiniteEnvironment`
            %   corresponding environments.
            %
            % eta : :class:`double`
            %   final error measure at convergence.
            
            arguments
                alg
                mpo
                mps
                envs = initialize_envs(mpo)
            end
            
            if length(mpo) ~= length(mps)
                error('dmrg:argerror', ...
                    'length of mpo (%d) should be equal to that of the mps (%d)', ...
                    length(mpo), length(mps));
            end
            
            t_total = tic;
            disp_init(alg);
            
            for iter = 1:alg.maxiter
                t_iter = tic;
                
                order = sweeporder(alg, length(mps));
                eta     = zeros(1, length(mps));
                lambda  = zeros(1, length(mps));
                
                for pos = sweeporder(alg, length(mps))
                    [mps, envs, lambda(pos), eta(pos)] = ...
                        localupdate(alg, mps, mpo, envs, order(pos));
                end
                
                eta = norm(eta, Inf);
                lambda = sum(lambda) / length(lambda);
                
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
        
        function [mps, envs, lambda, eta] = localupdate(alg, mps, mpo, envs, pos)
            % update orthogonality center
            mps = movegaugecenter(mps, pos);
            envs = movegaugecenter(envs, mpo, mps, mps, pos);
            
            % compute update
            H_AC = AC_hamiltonian(mpo, mps, envs, pos);
            AC = mps.A{pos};
            [AC, lambda] = eigsolve(alg.alg_eigs, @(x) H_AC.apply(x), AC, 1, alg.which);
            
            % determine error
            phase = dot(AC, mps.A{pos});
            eta = distance(AC ./ sign(phase), mps.A{pos});
            
            % take step
            mps.A{pos} = AC;
            envs = invalidate(envs, pos);
        end
    end
   
    %% Option handling
    methods
        function alg = updatetols(alg, iter, eta)
            if alg.dynamical_tols
                alg.alg_eigs.tol = between(alg.tol_min, eta * alg.eigs_tolfactor / iter, ...
                    alg.tol_max);
                if alg.verbosity > Verbosity.iter
                    fprintf('Updated eigsolver tolerances: (%e)\n', alg.alg_eigs.tol);
                end
            end
        end
        
        function order = sweeporder(alg, L)
            switch alg.sweepstyle
                case 'f2f'
                    order = [1:L L-1:-1:2];
                case 'b2b'
                    order = [L:-1:1 2:L];
                case 'm2m'
                    order = circshift([2:L, L-1:-1:1], -floor(L / 2));
                otherwise
                    error('unknown sweepstyle %s', alg.sweepstyle);
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
            fprintf('---- %s ----\n', alg.name);
        end
        
        function disp_iter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.iter, return; end
            
            s = settings;
            if abs(imag(lambda)) < eps(lambda)^(3/4) * abs(lambda)
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('%s %2d:\tE = %-0.4f\terror = %0.1e\t(%s)\n', ...
                            alg.name, iter, real(lambda), eta, time2str(t));
                    otherwise
                        fprintf('%s %4d:\tE = %-0.15f\terror = %0.4e\t(%s)\n', ...
                            alg.name, iter, real(lambda), eta, time2str(t, 's'));
                        
                end
            else
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('%s %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                            alg.name, iter, real(lambda), imag(lambda), eta, time2str(t));
                    otherwise
                        fprintf('%s %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                            alg.name, iter, real(lambda), imag(lambda), eta, time2str(t));
                        
                end
            end
        end
        
        function disp_conv(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.conv, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('%s converged %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                        alg.name, iter, real(lambda), imag(lambda), eta, time2str(t));
                otherwise
                    fprintf('%s converged %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                        alg.name, iter, real(lambda), imag(lambda), eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
        
        function disp_maxiter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.warn, return; end
            s = settings;
            switch s.matlab.commandwindow.NumericFormat.ActiveValue
                case 'short'
                    fprintf('%s max iterations %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                        alg.name, iter, real(lambda), imag(lambda), eta, time2str(t));
                otherwise
                    fprintf('%s max iterations %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
                        alg.name, iter, real(lambda), imag(lambda), eta, time2str(t));
                    
            end
            fprintf('---------------\n');
        end
    end
end

