classdef IDmrg2
    % `Infinite Density Matrix Renormalization Group algorithm with 2-site updates <https://arxiv.org/abs/0804.2509>`_.
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
    %   highest allowed convergence tolerance for soubroutines, defaults to :code:`1e-6`.
    %
    % eigs_tolfactor : :class:`double`
    %   relative scaling factor for determining the convergence tolerance of the local
    %   update solver subroutine based on the current error measure, defaults to
    %   :code:`1e-7`.
    %
    % trunc : :class:`struct`
    %   truncation method for local 2-site update, see :meth:`.Tensor.tsvd` for details on
    %   truncation options.
    %
    % alg_eigs : :class:`.KrylovSchur` or :class:`.Arnoldi`
    %   algorithm used for the eigsolver subroutines, defaults to
    %   :code:`Arnoldi('MaxIter', 100, 'KrylovDim', 20)`.
    
    properties
        tol = 1e-10
        miniter = 5
        maxiter = 100
        verbosity   = Verbosity.iter
        doplot = false
        
        which = 'largestabs'
        
        dynamical_tols          = true
        tol_min                 = 1e-12
        tol_max                 = 1e-6
        eigs_tolfactor          = 1e-7
        
        trunc                   = {'TruncDim', 10}
        
        doSave = false
        saveIterations = false
        saveMethod = 'full'
        name = 'IDmrg2'
        
        alg_eigs = Arnoldi('MaxIter', 100, 'KrylovDim', 20)
    end
    
    properties (Access = private)
        % TODO: alg_canonical, alg_environments? cf Vumps
        progressfig
    end
    
    methods
        function alg = IDmrg2(kwargs)
            arguments
                kwargs.?IDmrg2
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
        
        function [mps, lambda, GL, GR] = fixedpoint(alg, mpo, mps)
            % Find the fixed point MPS of an infinite MPO, given an initial guess.
            %
            % Usage
            % -----
            % :code:`[mps, lambda, GL, GR] = fixedpoint(alg, mpo, mps)`
            %
            % Arguments
            % ---------
            % alg : :class:`.IDmrg2`
            %   IDMRG2 algorithm.
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

            if period(mpo) ~= period(mps)
                error('idmrg2:argerror', ...
                    'periodicity of mpo (%d) should be equal to that of the mps (%d)', ...
                    period(mpo), period(mps));
            end
            if period(mps) < 2
                error('idmrg2:argerror', ...
                    'IDmrg2 is only defined for periodicity > 1');
            end
            
            t_total = tic;
            disp_init(alg);
            
            [GL, GR] = environments(mpo, mps, mps);
            
            for iter = 1:alg.maxiter
                t_iter = tic;
                
                C_ = mps.C{end};
                kwargs = {};
                lambdas = zeros(1, period(mps));
                
                % sweep from left to right
                for pos = 1:period(mps)-1
                    AC2 = computeAC2(mps, 1, pos);
                    H = AC2_hamiltonian(mpo, mps, GL, GR, pos);
                    [AC2, lambdas(pos)] = ...
                        eigsolve(alg.alg_eigs, @(x) H{1}.apply(x), AC2, 1, alg.which);
                    [mps.AL{pos}.var, C, mps.AR{pos+1}.var] = ...
                        tsvd(AC2.var, [1 2], [3 4], alg.trunc{:});
                    
                    mps.C{pos} = normalize(C);
                    mps.AC{pos + 1} = computeAC(mps, 1, pos + 1, 'R');
                    
                    TL = transfermatrix(mpo, mps, mps, pos, 'Type', 'LL');
                    GL{pos + 1} = apply(TL, GL{pos});
                    TR = transfermatrix(mpo, mps, mps, pos + 1, 'Type', 'RR');
                    GR{pos + 1} = apply(TR.', GR{mod1(pos + 2, length(GR))});
                end
                
                % update edge
                AC2 = contract(mps.AC{end}, [-1 -2 1], inv(mps.C{end}), [1 2], ...
                    mps.AL{1}, [2 -3 3], mps.C{1}, [3 -4], 'Rank', [2 2]);
                H = AC2_hamiltonian(mpo, mps, GL, GR, period(mps));
                [AC2, lambdas(end)] = ...
                    eigsolve(alg.alg_eigs, @(x) H{1}.apply(x), AC2, 1, alg.which);
                
                [mps.AL{end}.var, C, mps.AR{1}.var] = ...
                    tsvd(AC2, [1 2], [3 4], alg.trunc{:});
                
                mps.C{end} = normalize(C);
                mps.AC{end} = computeAC(mps, 1, period(mps), 'L');
                mps.AC{1} = computeAC(mps, 1, 1, 'R');
                mps.AL{1} = multiplyright(mps.AC{1}, inv(mps.C{1}));
                
                TL = transfermatrix(mpo, mps, mps, period(mps), 'Type', 'LL');
                GL{1} = apply(TL, GL{end});
                TR = transfermatrix(mpo, mps, mps, 1, 'Type', 'RR');
                GR{1} = apply(TR.', GR{2});
                
                % sweep from right to left
                for pos = period(mps)-1:-1:1
                    AC2 = computeAC2(mps, 1, pos, 'L');
                    H = AC2_hamiltonian(mpo, mps, GL, GR, pos);
                    [AC2, lambdas(pos)] = ...
                        eigsolve(alg.alg_eigs, @(x) H{1}.apply(x), AC2, 1, alg.which);
                    
                    [mps.AL{pos}.var, C, mps.AR{pos + 1}.var] = ...
                       tsvd(AC2, [1 2], [3 4], alg.trunc{:});
                    mps.C{pos} = normalize(C);
                    mps.AC{pos} = computeAC(mps, 1, pos, 'L');
                    mps.AC{pos + 1} = computeAC(mps, 1, pos + 1, 'R');
                    
                    TL = transfermatrix(mpo, mps, mps, pos, 'Type', 'LL');
                    GL{pos + 1} = apply(TL, GL{pos});
                    TR = transfermatrix(mpo, mps, mps, pos + 1, 'Type', 'RR');
                    GR{pos + 1} = apply(TR.', GR{mod1(pos + 2, period(mps))});
                end
                
                % update edge
                AC2 = contract(mps.C{end-1}, [-1 1], mps.AR{end}, [1 -2 2], ...
                    inv(mps.C{end}), [2 3], mps.AC{1}, [3 -3 -4], 'Rank', [2 2]);
                H = AC2_hamiltonian(mpo, mps, GL, GR, period(mps));
                [AC2, lambdas(1)] = ...
                    eigsolve(alg.alg_eigs, @(x) H{1}.apply(x), AC2, 1, alg.which);
                
                [mps.AL{end}.var, C, mps.AR{1}.var] = ...
                    tsvd(AC2, [1 2], [3 4], alg.trunc{:});
                mps.C{end} = normalize(C);
                mps.AC{1} = computeAC(mps, 1, 1, 'R');
                
                mps.AR{end} = multiplyleft(multiplyright(mps.AL{end}, mps.C{end}), ...
                    inv(mps.C{end - 1}));
                
                
                TL = transfermatrix(mpo, mps, mps, period(mps), 'Type', 'LL');
                GL{1} = apply(TL, GL{end});
                TR = transfermatrix(mpo, mps, mps, 1, 'Type', 'RR');
                GR{1} = apply(TR.', GR{2});
                
                % error measure
                infspace = infimum(space(C_, 1), space(mps.C{end}, 1));
                e1 = C_.eye(space(mps.C{end}, 1), infspace);
                e2 = C_.eye(space(C_, 1), infspace);
                
                eta = distance(e2' * C_ * e2, e1' * mps.C{end} * e1);
                lambda = prod(sqrt(lambdas));
                
                if iter > alg.miniter && eta < alg.tol
                    mps = canonicalize(mps, 'Order', 'rl');
                    disp_conv(alg, iter, lambda, eta, toc(t_total));
                    return
                end
                
                alg = updatetols(alg, iter, eta);
                plot(alg, iter, mps, eta);
                disp_iter(alg, iter, lambda, eta, toc(t_iter));

                if alg.doSave && mod(iter, alg.saveIterations) == 0
                    save(alg, mps, lambdas);
                end
            end
            
            mps = canonicalize(mps, 'Order', 'rl');
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
                alg.alg_eigs.tol = between(alg.tol_min, eta * alg.eigs_tolfactor, ...
                    alg.tol_max / iter);
                
                if alg.verbosity > Verbosity.iter
                    fprintf('Updated subalgorithm tolerances: (%e)\n', ...
                        alg.alg_eigs.tol);
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
            fprintf('---- IDmrg2 ----\n');
        end
        
        function disp_iter(alg, iter, lambda, eta, t)
            if alg.verbosity < Verbosity.iter, return; end
            
            s = settings;
            if abs(imag(lambda)) < eps(lambda)^(3/4) * abs(lambda)
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('IDmrg2 %2d:\tE = %-0.4f\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t));
                    otherwise
                        fprintf('IDmrg2 %4d:\tE = %-0.15f\terror = %0.4e\t(%s)\n', ...
                            iter, real(lambda), eta, time2str(t, 's'));
                        
                end
            else
                switch s.matlab.commandwindow.NumericFormat.ActiveValue
                    case 'short'
                        fprintf('IDmrg2 %2d:\tE = %-0.4f %+0.4fi\terror = %0.1e\t(%s)\n', ...
                            iter, real(lambda), imag(lambda), eta, time2str(t));
                    otherwise
                        fprintf('IDmrg2 %4d:\tE = %-0.15f %+0.15fi\terror = %0.4e\t(%s)\n', ...
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

        function save(alg, mps, lambdas)
            fileName = alg.name;

            fileData = struct;
            fileData.mps            = mps;
            fileData.lambdas         = lambdas;

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

