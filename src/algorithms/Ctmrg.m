classdef Ctmrg
    % `Corner Transfer Matrix Renormalisation Group algorithm <https://arxiv.org/abs/cond-mat/9705072>`_ for PEPS.
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
    %   maximum number of iteration, defaults to :code:`100`.
    %
    % projectortype : :class:`char`
    %   projector scheme used in the algorithm, currently only supports a single default
    %   value.
    %
    % trunc : :class:`struct`
    %   specification of truncation options, see :meth:`.Tensor.tsvd` for details.
    %
    % verbosity : :class:`.Verbosity`
    %   verbosity level of the algorithm, defaults to :code:`Verbosity.iter`.
    %
    % doplot : :class:`logical`
    %   plot progress, defaults to :code:`false`.
    
    %% Options
    properties

        tol = 1e-10
        miniter = 5
        maxiter = 100

        projectortype = 'cornersboris'
        trunc

        verbosity = Verbosity.iter        
        doPlot = false
        
        doSave = false
        saveIterations = 1
        saveMethod = 'full'
        name = 'CTMRG'

    end
    
    methods
        function alg = Ctmrg(kwargs)
            arguments
                kwargs.?Ctmrg
            end
            
            fields = fieldnames(kwargs);
            if ~isempty(fields)
                for field = fields.'
                    alg.(field{1}) = kwargs.(field{1});
                end
            end
            
        end
        
        function [envs, new_norm, err] = fixedpoint(alg, peps_top, peps_bot, envs)
            % Find the fixed point CTMRG environment of an infinite PEPS overlap.
            %
            % Usage
            % -----
            % :code:`[envs, new_norm, err] = fixedpoint(alg, peps_top, peps_bot, envs)`
            %
            % Arguments
            % ---------
            % alg : :class:`.Ctmrg`
            %   CTMRG algorithm.
            %
            % peps_top : :class:`.UniformPeps`
            %   top-layer uniform PEPS, usually interpreted as the 'ket' in the overlap.
            %
            % peps_bot : :class:`.UniformPeps`
            %   bottom-layer uniform PEPS, usually interpreted as the 'bra' in the overlap, 
            %   defaults to the top layer state.
            %
            % envs : :class:`.CtmrgEnv`
            %   initial guess for the fixed point CTMRG environment.
            %
            % Returns
            % -------
            % envs : :class:`.CtmrgEnv`
            %  fixed point CTMRG environment.
            %
            % new_norm : :class:`double`
            %   corresponding norm.
            %
            % err : :class:`double`
            %   final error measure at convergence.
            
            arguments
                alg
                peps_top
                peps_bot = peps_top
                envs = CtmrgEnvironment(peps_top, peps_bot)
            end

            if isempty(alg.trunc)
               alg.trunc  = struct('TruncSpace', space(envs.corners{1,:,:},1));
            end

            err = Inf;
            iter = 1;
        
            t_total = tic;
            disp_init(alg);

            old_norm = contract_ctrmg(alg, peps_top, peps_bot, envs);

            eta1 = 1.0;
            while (err > alg.tol && iter <= alg.maxiter) || iter <= alg.miniter

                t_iter = tic;

                eta = 0.0;
                for i = 1:4
                    if isfield(alg.trunc,'doPlot') && alg.trunc.doPlot
                        subplot(2,2,i)
                    end
                    [envs, eta0] = left_move(alg, peps_top, peps_bot, envs);
                    eta = max(eta, eta0);
                    envs = rot90(envs);
                    peps_top = rot90(peps_top);
                    peps_bot = rot90(peps_bot);
                end
        
                new_norm = contract_ctrmg(alg, peps_top, peps_bot, envs);
        
                err = abs(old_norm - new_norm);
                deta = abs((eta1 - eta) / eta1);

                if alg.doPlot
                    plot(alg, iter, envs, err);
                end

                disp_iter(alg, iter, err, abs(new_norm), eta, deta, toc(t_iter));
                if alg.doSave && mod(iter, alg.saveIterations) == 0
                    save_iteration(alg, envs, err, new_norm, iter, eta, toc(t_total));
                end                

                if iter > alg.miniter && err < alg.tol
                    disp_conv(alg, err, abs(new_norm), eta, deta, toc(t_total));
                    return
                end
        
                old_norm = new_norm;
                eta1 = eta;
                iter = iter + 1;
            end
        
            disp_maxiter(alg, err, abs(new_norm), eta, deta, toc(t_total));
        end
        
        function [envs, eta] = left_move(alg, peps_top, peps_bot, envs)

            eta = 0.0;
            h = height(peps_top);
            w = width(peps_top);

            for j = 1:w
        
                [above_projs, below_projs, etaj] = get_projectors(alg, peps_top, peps_bot, envs, j);                
                eta = max(eta, etaj);

                %use the projectors to grow the corners/edges
                for i = 1:h

                    envs.corners{1,prev(i,h),j} = contract(envs.corners{1,prev(i,h),prev(j,w)}, [2 1], ...
                                                    envs.edges{1,prev(i,h),j}, [1 3 4 -2], ...
                                                    above_projs{prev(i,h)}, [-1 4 3 2], ...
                                                    'Rank', [1,1]);

                    envs.edges{4,i,j} = contract(   envs.edges{4,i,prev(j,w)}, [7 2 4 1], ...
                                                    peps_top.A{i,j}.var, [6 2 8 -2 3], ...
                                                    conj(peps_bot.A{i,j}.var), [6 4 9 -3 5], ...
                                                    below_projs{prev(i,h)}, [1 3 5 -4], ...
                                                    above_projs{i}, [-1 9 8 7], ...
                                                    'Rank', [3,1]);

                    envs.corners{4,next(i,h),j} = contract(envs.corners{4,next(i,h),prev(j,w)}, [1 2], ...
                                                    envs.edges{3,next(i,h),j}, [-1 3 4 1], ...
                                                    below_projs{i}, [2 3 4 -2], ...
                                                    'Rank', [1,1]);
        
                end
    
            envs.corners(1,:,j) = cellfun(@(x) x ./ norm(x), envs.corners(1,:,j), 'UniformOutput', false);
            envs.edges(4,:,j) = cellfun(@(x) x ./ norm(x), envs.edges(4,:,j), 'UniformOutput', false);
            envs.corners(4,:,j) = cellfun(@(x) x ./ norm(x), envs.corners(4,:,j), 'UniformOutput', false);
                
            end
        end

        function [above_projs, below_projs, etaj] = get_projectors(alg, peps_top, peps_bot, envs, j)

            etaj = 0;
            h = height(peps_top);
            w = width(peps_top);
            above_projs = cell(h,1);
            below_projs = cell(h,1);

            switch alg.projectortype

                case 'cornersboris'

                    for i = 1:h

                        %SW corner
                        Q1 = contract(  envs.edges{3, mod1(i+2,h), j}, [-1 5 3 1], ...
                                        envs.corners{4, mod1(i+2,h), prev(j,w)}, [1 2], ...
                                        envs.edges{4, next(i,h), prev(j,w)}, [2 6 4 -6], ...
                                        peps_top.A{next(i,h), j}.var, [7 6 5 -2 -5], ...
                                        conj(peps_bot.A{next(i,h), j}.var), [7 4 3 -3 -4], ...
                                        'Rank', [3,3]);
                        %NW corner
                        Q2 = contract(  envs.edges{4, i, prev(j,w)}, [-1 3 5 2], ...
                                        envs.corners{1, prev(i,h), prev(j,w)}, [2 1], ...
                                        envs.edges{1, prev(i,h), j}, [1 4 6 -6], ...
                                        peps_top.A{i, j}.var, [7 3 -2 -5 4], ...
                                        conj(peps_bot.A{i, j}.var), [7 5 -3 -4 6], ...
                                        'Rank', [3,3]);

                        trun = [fieldnames(alg.trunc),struct2cell(alg.trunc)]';
                        [U, S, V, etaij] = tsvd(contract(Q1,[-1,-2,-3,3,2,1],Q2,[1,2,3,-4,-5,-6],'Rank',[3,3]), trun{:});                        

                        etaj = max(etaj, etaij);            
                        isqS = inv(sqrtm(S));    
                        Q = isqS*U'*Q1;
                        P = Q2*V'*isqS;

                        %norm(contract(Q,[-1,1,2,3],P,[3,2,1,-2],'Rank',[1,1])-Tensor.eye(space(Q,1),space(Q,1)))
                        %norm(P*Q-Tensor.eye(space(Q,2:4)',space(Q,2:4)'))
    
                        above_projs{i} = Q;
                        below_projs{i} = P;                %should be ok for any arrow

                    end

            end


        end

        function total = contract_ctrmg(alg, peps_top, peps_bot, envs)

            total = 1.0;
            h = height(peps_top);
            w = width(peps_top);

            for i = 1:h
                for j = 1:w
                    total = total * contract(   envs.corners{1,prev(i,h),prev(j,w)}, [9 1], ...
                                                envs.edges{1,prev(i,h),j}, [1 4 7 2], ...
                                                envs.corners{2,prev(i,h),next(j,w)}, [2 3], ...
                                                envs.edges{2,i,next(j,w)}, [3 5 8 17], ...
                                                envs.corners{3,next(i,h),next(j,w)}, [17 16], ...
                                                envs.edges{3,next(i,h),j}, [16 14 15 13], ...
                                                envs.corners{4,next(i,h),prev(j,w)}, [13 12], ...
                                                envs.edges{4,i,prev(j,w)}, [12 10 11 9], ...
                                                peps_top.A{i,j}.var, [6 10 14 5 4], ...
                                                conj(peps_bot.A{i,j}.var), [6 11 15 8 7]);
                    
                    total = total * contract(   envs.corners{1,prev(i,h),prev(j,w)}, [4 1], ...
                                                envs.corners{2,prev(i,h),j}, [1 2], ...
                                                envs.corners{3,i,j}, [2 3], ...
                                                envs.corners{4,i,prev(j,w)}, [3 4]);
        
                    total = total / contract(   envs.corners{1,prev(i,h),prev(j,w)}, [2 1], ...
                                                envs.corners{2,prev(i,h),j}, [1 3], ...
                                                envs.edges{2,i,j}, [3 4 5 6], ...
                                                envs.corners{3,next(i,h),j}, [6 7], ...
                                                envs.corners{4,next(i,h),prev(j,w)}, [7 8], ...
                                                envs.edges{4,i,prev(j,w)}, [8 4 5 2]);

                    total = total / contract(   envs.corners{1,prev(i,h),prev(j,w)}, [1 3], ...
                                                envs.corners{4,i,prev(j,w)}, [2 1], ...
                                                envs.edges{1,prev(i,h),j}, [3 4 5 6], ...
                                                envs.corners{2,prev(i,h),next(j,w)}, [6 7], ...
                                                envs.corners{3,i,next(j,w)}, [7 8], ...
                                                envs.edges{3,i,j}, [8 4 5 2]);
                end
            end

        end
    
    end

    %% Display
    methods(Access = private)
        
        function plot(alg, iter, envs, eta)
            if ~alg.doPlot, return; end
            persistent axhistory axspectrum
            
            D = height(envs);
            W = width(envs);
            
            if iter == 1
                figure('Name', 'Ctmrg');
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
            plot_entanglementspectrum(UniformMps({MpsTensor(envs.edges{1,:,:})}), 1:D, 1:W, axspectrum);
            drawnow
        end
        
        function disp_init(alg)
            if alg.verbosity < Verbosity.conv, return; end
            fprintf('---- CTMRG ----\n');
        end

        function disp_iter(alg, iter, err, norm, eta, deta, t)
            if alg.verbosity < Verbosity.iter, return; end
            fprintf('iteration: %4d\t\terror: %.2e\t\tnorm: %.10e\t\teta: %.2e\t\tdeta: %.2e\t(%s)\n', ...
                iter, err, norm, eta, deta, time2str(t));
        end
        
        function disp_conv(alg, err, norm, eta, deta, t)
            if alg.verbosity < Verbosity.conv, return; end
            fprintf('CTMRG converged \t\terror: %.2e\t\tnorm: %.10e\t\teta: %.2e\t\tdeta: %.2e\t(%s)\n', ...
                err, norm, eta, deta, time2str(t));
            fprintf('---------------\n');
        end
        
        function disp_maxiter(alg, err, norm, eta, deta, t)
            if alg.verbosity < Verbosity.warn, return; end
            fprintf('CTMRG max iterations \terror: %.2e\t\tnorm: %.10e\t\teta: %.2e\t\tdeta: %.2e\t(%s)\n', ...
                err, norm, eta, deta, time2str(t));
            fprintf('---------------\n');
        end

        function save_iteration(alg, envs, err, norm, iter, eta, t)
            fileName = alg.name;
            fileData = struct;
            fileData.envs = envs;
            fileData.err = err;
            fileData.norm = norm;
            fileData.iteration = iter;
            fileData.eta = eta;
            fileData.time = t;            
            save(fileName, '-struct', 'fileData', '-v7.3');
        end
    end
      
end

