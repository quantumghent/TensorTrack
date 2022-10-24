classdef IDmrg2
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
        trunc                   = {'TruncDim', 10}
    end
    
    properties (Access = private)
        alg_eigs = struct('MaxIter', 100, 'KrylovDim', 20)
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
                alg.alg_eigs.Tol = sqrt(alg.tol_min * alg.tol_max);
                alg.alg_eigs.Verbosity = alg.verbosity - 2;
            end
        end
        
        function [mps, lambda, GL, GR] = fixedpoint(alg, mpo, mps)
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
                
                C_ = mps.C(end);
                kwargs = {};
                lambdas = zeros(1, period(mps));
                
                % sweep from left to right
                for pos = 1:period(mps)-1
                    AC2 = contract(mps.AC(pos), [-1 -2 1], ...
                        mps.AR(pos + 1), [1 -3 -4], ...
                        'Rank', [2 2]);
                    H = AC2_hamiltonian(mpo, mps, GL, GR, pos);
                    [AC2, lambdas(pos)] = ...
                        eigsolve(H{1}, AC2, 1, alg.which, kwargs{:});
                    [mps.AL(pos), C, mps.AR(pos + 1), delta] = ...
                        tsvd(AC2, [1 2], [3 4], alg.trunc{:});
                    
                    mps.C(pos) = normalize(C);
                    mps.AC(pos + 1) = multiplyleft(mps.AR(pos + 1), mps.C(pos));
                    
                    TL = transfermatrix(mpo, mps, mps, pos, 'Type', 'LL');
                    GL{pos + 1} = apply(TL, GL{pos}) / sqrt(lambdas(pos));
                    TR = transfermatrix(mpo, mps, mps, pos + 1, 'Type', 'RR');
                    GR{pos + 1} = apply(TR.', GR{mod1(pos + 2, length(GR))}) / ...
                        sqrt(lambdas(pos));
                end
                
                % update edge
                AC2 = contract(mps.AC(end), [-1 -2 1], inv(mps.C(end)), [1 2], ...
                    mps.AL(1), [2 -3 3], mps.C(1), [3 -4], 'Rank', [2 2]);
                H = AC2_hamiltonian(mpo, mps, GL, GR, period(mps));
                [AC2, lambdas(end)] = eigsolve(H{1}, AC2, 1, alg.which, kwargs{:});
                
                [mps.AL(end), C, mps.AR(1)] = ...
                    tsvd(AC2, [1 2], [3 4], alg.trunc{:});
                mps.C(end) = normalize(C);
                mps.AC(end) = multiplyright(mps.AL(end), mps.C(end));
                mps.AC(1) = multiplyleft(mps.AR(1), mps.C(end));
                mps.AL(1) = multiplyright(mps.AC(1), inv(mps.C(1)));
                
                TL = transfermatrix(mpo, mps, mps, period(mps), 'Type', 'LL');
                GL{1} = apply(TL, GL{end}) / sqrt(lambdas(end));
                TR = transfermatrix(mpo, mps, mps, 1, 'Type', 'RR');
                GR{1} = apply(TR.', GR{2}) / sqrt(lambdas(end));
                
                % sweep from right to left
                for pos = period(mps)-1:-1:1
                    AC2 = contract(mps.AL(pos), [-1 -2 1], ...
                        mps.AC(pos + 1), [1 -3 -4], 'Rank', [2 2]);
                    H = AC2_hamiltonian(mpo, mps, GL, GR, pos);
                    [AC2, lambdas(pos)] = eigsolve(H{1}, AC2, 1, alg.which, kwargs{:});
                    
                    [mps.AL(pos), C, mps.AR(pos + 1)] = ...
                        tsvd(AC2, [1 2], [3 4], alg.trunc{:});
                    mps.C(pos) = normalize(C);
                    mps.AC(pos) = multiplyright(mps.AL(pos), mps.C(pos));
                    mps.AC(pos + 1) = multiplyleft(mps.AR(pos + 1), mps.C(pos));
                    
                    TL = transfermatrix(mpo, mps, mps, pos, 'Type', 'LL');
                    GL{pos + 1} = apply(TL, GL{pos}) / sqrt(lambdas(pos));
                    TR = transfermatrix(mpo, mps, mps, pos + 1, 'Type', 'RR');
                    GR{pos + 1} = apply(TR.', GR{mod1(pos + 2, period(mps))}) / ...
                        sqrt(lambdas(pos));
                end
                
                % update edge
                AC2 = contract(mps.C(end-1), [-1 1], mps.AR(end), [1 -2 2], ...
                    inv(mps.C(end)), [2 3], mps.AC(1), [3 -3 -4], 'Rank', [2 2]);
                H = AC2_hamiltonian(mpo, mps, GL, GR, period(mps));
                [AC2, lambdas(1)] = eigsolve(H{1}, AC2, 1, alg.which, kwargs{:});
                
                [mps.AL(end), C, mps.AR(1)] = ...
                    tsvd(AC2, [1 2], [3 4], alg.trunc{:});
                mps.C(end) = normalize(C);
                mps.AC(1) = multiplyleft(mps.AR(1), mps.C(end));
                mps.AR(end) = multiplyleft(multiplyright(mps.AL(end), mps.C(end)), ...
                    inv(mps.C(end - 1)));
                
                TL = transfermatrix(mpo, mps, mps, period(mps), 'Type', 'LL');
                GL{1} = apply(TL, GL{end}) / sqrt(lambdas(end));
                TR = transfermatrix(mpo, mps, mps, 1, 'Type', 'RR');
                GR{1} = apply(TR.', GR{2}) / sqrt(lambdas(end));
                
                % error measure
                infspace = infimum(space(C_, 1), space(mps.C(end), 1));
                e1 = C_.eye(space(mps.C(end), 1), infspace);
                e2 = C_.eye(space(C_, 1), infspace);
                
                eta = distance(e2' * C_ * e2, e1' * mps.C(end) * e1);
                lambda = prod(sqrt(lambdas));
                
                if iter > alg.miniter && eta < alg.tol
                    mps = canonicalize(mps, 'Order', 'rl');
                    disp_conv(alg, iter, lambda, eta, toc(t_total));
                    return
                end
                
                alg = updatetols(alg, iter, eta);
                disp_iter(alg, iter, lambda, eta, toc(t_iter));
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
    end
end

