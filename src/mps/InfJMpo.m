classdef InfJMpo < InfMpo
    % Infinite Mpo with a Jordan block structure
    
    methods
        function mpo = InfJMpo(varargin)
            mpo@InfMpo(varargin{:});
            if nargin > 0
                assert(istriu(mpo.O{1}));
                assert(iseye(mpo.O{1}(1, 1, 1, 1)) && iseye(mpo.O{1}(end, 1, end, 1)));
                assert(isconnected(mpo));
            end
        end
        
        function bool = isconnected(mpo)
            bool = true;
        end
        
        function [GL, lambda] = leftenvironment(mpo, mps1, mps2, GL, linopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GL = cell(1, period(mps1))
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            linkwargs = namedargs2cell(linopts);
            L = period(mps1);
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'LL');
            fp_left = fixedpoint(mpo, mps1, 'l_LL', 1);
            fp_right = fixedpoint(mpo, mps1, 'r_LL', L);
            
            if isempty(GL) || isempty(GL{1}) % initialize environment
                GL = cell(1, period(mps1));
                GL{1} = repartition(MpsTensor(SparseTensor.zeros(domain(T), [])), ...
                    rank(fp_left));
                GL{1}.var(1) = fp_left;
            end
            
            for i = 2:size(GL{1}, 2)
                rhs = apply(slice(T, i, 1:i-1), slice(GL{1}, 1, 1:i-1, 1));
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GL{1}.var(i) = rhs;
                else
                    if iseye(T, i)
                        lambda = overlap(rhs, fp_right);
                        rhs = rhs - lambda * fp_left;
                    end
                    [GL{1}.var(i), flag, relres, iter] = linsolve(@(x) x - apply(Tdiag, x), rhs, GL{1}.var(i), ...
                        linkwargs{:});
                    if iseye(T, i)
                        GL{1}.var(i) = GL{1}.var(i) - overlap(GL{1}.var(i), fp_right) * fp_left;
                    end
                    
                    if flag == 1 && linopts.Verbosity >= Verbosity.warn
                        fprintf('GL(%d) reached maxiter, error = %3e\n', i, relres);
                    elseif flag == 2 && linopts.Verbosity >= Verbosity.warn
                        fprintf('GL(%d) stagnated, error = %3e\n', i, relres);
                    elseif linopts.Verbosity >= Verbosity.conv
                        fprintf('GL(%d) converged after %d iterations, error = %3e\n', i, iter, relres);
                    end
                end
            end
            
            if nnz(GL{1}.var) == numel(GL{1}.var)
                GL{1}.var = full(GL{1}.var);
            end
            
            for w = 1:period(mps1)-1
                T = transfermatrix(mpo, mps1, mps2, w, 'Type', 'LL');
                GL{next(w, period(mps1))} = apply(T, GL{w});
            end
        end
        
        function [GR, lambda] = rightenvironment(mpo, mps1, mps2, GR, linopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GR = cell(1, period(mps1))
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            linkwargs = namedargs2cell(linopts);
            L = period(mpo);
            
            fp_left = fixedpoint(mpo, mps1, 'l_RR', 1);
            fp_right = repartition(fixedpoint(mpo, mps1, 'r_RR', L), [3 0]);
            
            T = transfermatrix(mpo, mps1, mps2, 'Type', 'RR').';
            N = size(T(1).O{1}, 2);
            
            if isempty(GR) || isempty(GR{1})
                GR = cell(1, period(mps1));
                GR{1} = MpsTensor(SparseTensor.zeros(domain(T), []));
                GR{1}.var(1, N, 1) = fp_right;
            end
            
            for i = N-1:-1:1
                rhs = apply(slice(T, i, i+1:N), slice(GR{1}, 1, i+1:N, 1));
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GR{1}.var(i) = rhs;
                else
                    if iseye(T, i)
                        lambda = overlap(rhs, fp_left);
                        rhs = rhs - lambda * fp_right;
                    end
                    [GR{1}.var(i), flag, relres, iter] = ...
                        linsolve(@(x) x - apply(Tdiag, x), rhs, GR{1}.var(i), linkwargs{:});
                    if iseye(T, i)
                        GR{1}.var(i) = GR{1}.var(i) - ...
                            overlap(GR{1}.var(i), fp_left) * fp_right;
                    end
                    
                    if flag == 1 && linopts.Verbosity >= Verbosity.warn
                        fprintf('GR(%d) reached maxiter, error = %3e\n', i, relres);
                    elseif flag == 2 && linopts.Verbosity >= Verbosity.warn
                        fprintf('GR(%d) stagnated, error = %3e\n', i, relres);
                    elseif linopts.Verbosity >= Verbosity.conv
                        fprintf('GR(%d) converged after %d iterations, error = %3e\n', i, iter, relres);
                    end
                end
            end
            
            if nnz(GR{1}.var) == numel(GR{1}.var)
                GR{1}.var = full(GR{1}.var);
            end
            
            for w = period(mps1):-1:2
                T = transfermatrix(mpo, mps1, mps2, w, 'Type', 'RR').';
                GR{w} = apply(T, GR{next(w, period(mps1))});
            end
        end
        
        function GBL = leftquasienvironment(mpo, qp, GL, GR, linopts)
            arguments
                mpo
                qp
                GL
                GR
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(qp))^(3/5)
            end
            
            linkwargs = namedargs2cell(linopts);
            expP = exp(-1i * qp.p);
            L = period(mpo);
            
            needsRegularization = istrivial(qp);
            if needsRegularization || true
                fp_left  = fixedpoint(mpo, qp, 'l_RL_0', 1);
                fp_right = fixedpoint(mpo, qp, 'r_RL_1', L);
            end
            
            T = transfermatrix(mpo, qp, qp, 'Type', 'RL');
            TB = transfermatrix(mpo, qp, qp, 'Type', 'BL');
            
            % initialize and precompute GL * TB
            GBL = cell(size(GL));
            GBL{1} = MpsTensor(SparseTensor.zeros(domain(T), auxspace(qp, 1)), 1);
            for w = 1:L
                GBL{next(w, L)} = ...
                    (apply(TB(w), GL{w}) + apply(T(w), GBL{w})) * (expP^(1/L));
            end
            
            N = size(GBL{1}, 2);
            for i = 2:N % GBL{1}(1) = 0 because of quasiparticle gauge
                rhs = apply(slice(T, i, 1:i-1), slice(GBL{1}, 1, 1:i-1, 1, 1)) * expP;
                rhs = rhs + slice(GBL{1}, i);
                
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GBL{1}.var(i) = rhs;
                else
                    if needsRegularization && iseye(T, i)
                        lambda = overlap(rhs, fp_right);
                        fp_left = repartition(fp_left, rank(rhs));
                        rhs = rhs - lambda * fp_left;
                        H_effective = @(x) x - expP * ...
                            apply_regularized(Tdiag, fp_left, fp_right, x);
                    else
                        H_effective = @(x) x - expP * apply(Tdiag, x);
                    end
                    
                    [GBL{1}.var(i), flag, relres, iter, resvec] = ...
                        linsolve(H_effective, rhs, [], linkwargs{:});
                    
                    if flag == 1 && linopts.Verbosity >= Verbosity.warn
                        fprintf('GBL(%d) reached maxiter, error = %3e\n', i, relres);
                    elseif flag == 2 && linopts.Verbosity >= Verbosity.warn
                        fprintf('GBL(%d) stagnated, error = %3e\n', i, relres);
                    elseif linopts.Verbosity >= Verbosity.conv
                        fprintf('GBL(%d) converged after %.1f iterations, error = %3e\n', i, iter, relres);
                    end
                end
            end
            
            if nnz(GBL{1}) == numel(GBL{1}.var)
                GBL{1}.var = full(GBL{1}.var);
            end
            
            GBL{1} = MpsTensor(GBL{1}, 1);
            for w = 1:L-1
                GBL{next(w, L)} = expP^(1 / L) * ...
                    (apply(TB(w), GL{w}) + apply(T(w), GBL{w}));
            end
        end
        
        function GBR = rightquasienvironment(mpo, qp, GL, GR, linopts)
            arguments
                mpo
                qp
                GL
                GR
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(qp))^(3/5)
            end
            
            linkwargs = namedargs2cell(linopts);
            expP = exp(+1i*qp.p);
            L = period(mpo);
            
            needsRegularization = istrivial(qp);
            if needsRegularization || true
                fp_left  = fixedpoint(mpo, qp, 'l_LR_1', 1);
                fp_right = repartition(fixedpoint(mpo, qp, 'r_LR_0', L), [3, 1]);
            end
            
            T = transfermatrix(mpo, qp, qp, 'Type', 'LR').';
            TB = transfermatrix(mpo, qp, qp, 'Type', 'BR').';
            
            GBR = cell(size(GR));
            GBR{1} = MpsTensor(SparseTensor.zeros(domain(T), auxspace(qp, 1)), 1);
            for w = L:-1:1
                ww = next(w, L);
                TB_w = transfermatrix(mpo, qp, qp, w, 'Type', 'BR').';
                T_w = transfermatrix(mpo, qp, qp, w, 'Type', 'LR').';
                GBR{w} = expP^(1/L) * (apply(TB_w, GR{ww}) + apply(T_w, GBR{ww}));
            end
            
            N = size(GBR{1}, 2);
            for i = N:-1:1
                if i == N
                    rhs = slice(GBR{1}, 1, i, 1, 1);
                else
                    rhs = apply(slice(T, i, i+1:N), slice(GBR{1}, 1, i+1:N, 1, 1)) * expP;
                    rhs = rhs + slice(GBR{1}, 1, i, 1, 1);
                end
                
                Tdiag = slice(T, i, i);
                if iszero(Tdiag)
                    GBR{1}.var(i) = rhs;
                else
                    if needsRegularization && iseye(T, i)
                        lambda = overlap(rhs, fp_left);
                        fp_right = repartition(fp_right, rank(rhs));
                        rhs = rhs - lambda * fp_right;
                        H_effective = @(x) x - expP * ...
                            apply_regularized(Tdiag, fp_right, fp_left, x);
                    else
                        H_effective = @(x) x - expP * apply(Tdiag, x);
                    end
                    
                    [GBR{1}.var(i), flag, relres, iter, resvec] = ...
                        linsolve(H_effective, rhs, [], linkwargs{:});
                    
                    if flag == 1 && linopts.Verbosity >= Verbosity.warn
                        fprintf('GBR(%d) reached maxiter, error = %3e\n', i, relres);
                    elseif flag == 2 && linopts.Verbosity >= Verbosity.warn
                        fprintf('GBR(%d) stagnated, error = %3e\n', i, relres);
                    elseif linopts.Verbosity >= Verbosity.conv
                        fprintf('GBR(%d) converged after %.1f iterations, error = %3e\n', i, iter, relres);
                    end
                end
            end
            
            if nnz(GBR{1}) == numel(GBR{1})
                GBR{1} = full(GBR{1});
            end
            
            GBR{1} = MpsTensor(GBR{1}, 1);
            for w = L:-1:2
                ww = next(w, L);
                TB_w = transfermatrix(mpo, qp, qp, w, 'Type', 'BR').';
                T_w = transfermatrix(mpo, qp, qp, w, 'Type', 'LR').';
                GBR{w} = expP^(1/L) * (apply(TB_w, GR{ww}) + apply(T_w, GBR{ww}));
            end
        end
        
        function [GL, GR, lambda] = environments(mpo, mps1, mps2, GL, GR, linopts)
            arguments
                mpo
                mps1
                mps2 = mps1
                GL = cell(1, period(mps1))
                GR = cell(1, period(mps1))
                linopts.Algorithm = 'bicgstab'
                linopts.MaxIter = 500
                linopts.Verbosity = Verbosity.warn
                linopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            kwargs = namedargs2cell(linopts);
            [GL, lambdaL] = leftenvironment(mpo, mps1, mps2, GL, kwargs{:});
            [GR, lambdaR] = rightenvironment(mpo, mps1, mps2, GR, kwargs{:});
            lambda = (lambdaL + lambdaR) / 2;
            if ~isapprox(lambdaL, lambdaR, 'AbsTol', eps^(1/3), 'RelTol', eps(lambda)^(1/3))
                warning('lambdas disagree (%e, %e)', lambdaL, lambdaR);
            end
        end
        
        function mpo = renormalize(mpo, lambda)
            mpo = mpo - lambda;
        end
        
        function mpo = plus(a, b)
            if isa(a, 'InfJMpo') && isnumeric(b)
                if period(a) > 1 && isscalar(b)
                    b = repmat(b, 1, period(a));
                end
                
                for i = 1:period(a)
                    a.O{i}(1, 1, end, 1) = a.O{i}(1, 1, end, 1) + b(i);
                end
                mpo = a;
                
            elseif isnumeric(a) && isa(b, 'InfJMpo')
                mpo = b + a;
            end
        end
        
        function mpo = minus(a, b)
            mpo = a + (-b);
        end
        
        function mpo = mtimes(mpo, b)
            if isnumeric(mpo) || isnumeric(b)
                mpo = mpo .* b;
                return
            end
            
            mpo = [mpo; b];
        end
        
        function finitempo = open_boundary_conditions(mpo, L)
            Os = repmat(mpo.O, 1, L);
            
            Os{1}   = Os{1}(1, :, :, :);
            Os{end} = Os{end}(:, :, end, :);
            
            rspace = subspaces(rightvspace(mpo, period(mpo)), size(Os{end}, 3));
            rightedge = MpsTensor(Tensor.eye([one(rspace) rspace'], one(rspace)));
            
            lspace = subspaces(leftvspace(mpo, 1), size(Os{1}, 1));
            leftedge = MpsTensor(Tensor.eye([one(lspace) lspace], one(lspace))');
            
            finitempo = FiniteMpo(leftedge, Os, rightedge);
        end
    end
    
    methods (Static)
        function mpo = twosite(Htwosite, Honesite, kwargs)
            arguments
                Htwosite
                Honesite = []
                kwargs.Trunc
            end
            newkwargs = namedargs2cell(kwargs);
            local_ops = MpoTensor.decompose_local_operator(Htwosite, newkwargs{:});
            L = local_ops{1};
            R = local_ops{2};
            
            P = pspace(L);
            assert(P == pspace(R), 'operators:spacemismatch', ...
                sprintf('incompatible physical spaces %s and %s', pspace(L), pspace(R)));
            
            total_leftvspace = [leftvspace(L) leftvspace(R) rightvspace(R)'];
            total_rightvspace = [leftvspace(L)' rightvspace(L) rightvspace(R)];
            assert(all(conj(total_rightvspace) == total_leftvspace))
            
            cod = SumSpace(total_leftvspace, P);
            dom = SumSpace(P, conj(total_rightvspace));
            
            O = MpoTensor.zeros(cod, dom);
            O(1, 1, 1, 1) = 1;
            O(3, 1, 3, 1) = 1;
            O(1, 1, 2, 1) = L;
            O(2, 1, 3, 1) = R;
            
            if ~isempty(Honesite)
                local_op = MpoTensor.decompose_local_operator(Honesite, newkwargs{:});
                assert(leftvspace(local_op{1}) == subspaces(leftvspace(O), 1) && ...
                    rightvspace(local_op{1}) == subspaces(rightvspace(O), 3) && ...
                    pspace(local_op{1}) == subspaces(pspace(O), 1), ...
                    'operators:spacemismatch', ...
                    'onesite operator incompatible with twosite operator.');
                O(1, 1, 3, 1) = local_op{1};
            end
            
            mpo = InfJMpo(O);
        end
    end
end
