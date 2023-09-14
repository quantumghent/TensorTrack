classdef FiniteMpo
    % Finite Matrix product operators
    
    properties
        L MpsTensor
        O
        R MpsTensor
    end
    
    methods
        function mpo = FiniteMpo(L, O, R)
            if nargin == 0, return; end
            if ~iscell(O), O = {O}; end
            mpo.L = L;
            mpo.O = O;
            mpo.R = R;
        end
        
        function v = apply(mpo, v)
            N = length(mpo);
            for d = 1:depth(mpo)
                if N == 0
                    v = applytransfer(mpo(d).L, mpo(d).R, v);
                elseif N == 1
                    v = applychannel(mpo(d).O{1}, mpo(d).L, mpo(d).R, v);
                else
                    v = applympo(mpo(d).O{:}, mpo(d).L, mpo(d).R, v);
                end
            end
        end
        
        function v = apply_regularized(mpo, fp1, fp2, v)
            v = apply(mpo, v);
            v = v - overlap(v, fp2) * repartition(fp1, rank(v));
        end
        
        function varargout = eigsolve(mpo, v0, howmany, sigma, options)
            arguments
                mpo
                v0 = []
                howmany = 1
                sigma = 'largestabs'
                
                options.Tol = eps(underlyingType(mpo))^(3/4)
                options.Algorithm = 'Arnoldi'
                options.MaxIter = 100
                options.KrylovDim = 20
                options.IsSymmetric logical = false
                options.DeflateDim
                options.ReOrth = 2
                options.NoBuild
                options.Verbosity = 0
            end
            
            if isempty(v0), v0 = initialize_fixedpoint(mpo(1)); end
            
            kwargs = namedargs2cell(options);
            [varargout{1:nargout}] = ...
                eigsolve(@(x) mpo.apply(x), v0, howmany, sigma, kwargs{:});
        end
        
        function v = initialize_fixedpoint(mpo)
            % Initialize a dense tensor for the fixedpoint of a :class:`FiniteMPO`.
            
            N = prod(cellfun(@(x) size(x, 4), mpo.O));
            for i = N:-1:1
                v(i) = MpsTensor.randnc(domain(slice(mpo, i, 1:N)), []);
            end
        end
        
        function mps = initialize_mps(mpo, kwargs)
            arguments
                mpo
                kwargs.MaxVspace
            end
            
            pspaces = arrayfun(@(x) pspace(mpo, x), 1:length(mpo), 'UniformOutput', false);
            
            vspacefirst = rightvspace(mpo.L)';
            vspacelast  = leftvspace(mpo.R);
            
            newkwargs = namedargs2cell(kwargs);
            mps = FiniteMps.randnc(pspaces{:}, 'LeftVspace', vspacefirst, ...
                'RightVspace', vspacelast, newkwargs{:});
            mps = normalize(mps);
        end
        
        function envs = initialize_envs(mpo)
            arguments
                mpo
            end
            
            GL = cell(1, length(mpo) + 1);
            GR = cell(1, length(mpo) + 1);
            
            GL{1} = mpo.L;
            GR{end} = mpo.R;
            
            envs = FiniteEnvironment(GL, GR);
        end
        
        function T = transfermatrix(mpo, mps1, mps2, sites)
        arguments
                mpo
                mps1
                mps2 = mps1
                sites = 1:length(mps1)
            end
            
            assert(all(diff(sites) == 1), 'sites must be neighbouring and increasing.');
            A1 = mps1.A(sites);
            A2 = mps2.A(sites);
            O = mpo.O(sites);                                                               %#ok<PROPLC>
            T = FiniteMpo.mps_channel_operator(A1, O, A2);                                  %#ok<PROPLC>
        end
        
        function H = AC_hamiltonian(mpo, mps, envs, pos)
            envs = movegaugecenter(envs, mpo, mps, mps, pos);
            H = FiniteMpo(envs.GL{pos}, mpo.O(pos), envs.GR{pos + 1});
        end
        
        function s = pspace(mpo, x)
            s = pspace(mpo.O{x});
        end
        
        function s = domain(mpo)
            sO = flip(cellfun(@(x) domainspace(x), mpo(1).O, 'UniformOutput', false));
            s = [leftvspace(mpo(1).R) [sO{:}] rightvspace(mpo(1).L)]';
        end
        
        function s = codomain(mpo)
            sO = cellfun(@(x) codomainspace(x), mpo(1).O, 'UniformOutput', false);
            s = [leftvspace(mpo(end).L) [sO{:}] rightvspace(mpo(end).R)];
        end
        
        function d = depth(mpo)
            d = builtin('length', mpo);
        end
        
        function l = length(mpo)
            l = length(mpo(1).O);
        end
        
        function type = underlyingType(mpo)
            type = underlyingType(mpo(1).L);
        end
        
        function mpo = ctranspose(mpo)
            if depth(mpo) > 1
                mpo = flip(mpo);
            end
            
            for d = 1:depth(mpo)
                mpo(d).L = tpermute(mpo(d).L', ...
                    [1, flip(2:nspaces(mpo(d).L)-1), nspaces(mpo(d).L)]);
                mpo(d).O = cellfun(@ctranspose, mpo(d).O, ...
                    'UniformOutput', false);
                mpo(d).R = tpermute(mpo(d).R', ...
                    [1, flip(2:nspaces(mpo(d).R)-1), nspaces(mpo(d).R)]);
            end
        end
        
        function mpo = transpose(mpo)
            if depth(mpo) > 1
                mpo = flip(mpo);
            end
            
            for d = 1:depth(mpo)
                [mpo(d).L, mpo(d).R] = swapvars(mpo(d).L, mpo(d).R);
                mpo(d).O = cellfun(@transpose, ...
                    fliplr(mpo(d).O), 'UniformOutput', false);
            end
        end
        
        function mpo = slice(mpo, i, j)
            if isempty(mpo(1).O), return; end
            if strcmp(i, ':')
                i = 1:size(mpo(end).O, 2);
            end
            if strcmp(j, ':')
                j = 1:size(mpo(1).O{1}, 4);
            end
            assert(all(1 <= i) && all(i <= size(mpo(end).O{1}, 2)));
            assert(all(1 <= j) && all(j <= size(mpo(1).O{1}, 4)));
            
            mpo(end).O{1} = mpo(end).O{1}(:, i, :, :);
            mpo(1).O{1} = mpo(1).O{1}(:, :, :, j);
        end
        
        function bool = iszero(mpo)
            if isempty(mpo(1).O)
                bool = false;
                return
            end
            for i = 1:depth(mpo)
                bool = any(cellfun(@nnz, mpo(i).O) == 0);
                if bool, return; end
            end
        end
        
        function bool = iseye(mpo, i)
            if isempty(mpo(1).O)
                bool = true;
                return
            end
            for d = 1:depth(mpo)
                bool = all(cellfun(@(x) iseye(x(:,i,:,i)), mpo(d).O));
                if ~bool, return; end
            end
        end
%         
%         function v = applyleft(mpo, v)
%             
%         end
%         
%         function mpo1 = plus(mpo1, mpo2)
%             
%         end
%         
%         function v = applyright(mpo, v)
%             arguments
%                 mpo
%                 v MpsTensor
%             end
%             
%             assert(depth(mpo) == 1, 'mps:TBA', 'Not implemented yet.');
%             assert(v.plegs == width(mpo), 'mps:ArgError', 'Incompatible sizes.');
%             
%             w = width(mpo);
%             
%             mpopart = cell(2, w);
%             for i = 1:w
%                 mpopart{1, i} = mpo.O{i};
%                 mpopart{2, i} = [2 * i, 2 * (i + 1) + 1, -(1 + i), 2 * i + 1];
%             end
%             
%             v = MpsTensor(contract(...
%                 v, [1, 2:2:(2 * w), 2 * (w + 1), -(1:v.alegs) - (w + 2)], ...
%                 mpo.L, [-1 3 1], ...
%                 mpo.R, [-(w + 2), 2 * (w + 1) + 1, 2 * (w + 1)], ...
%                 mpopart{:}, 'Rank', rank(v)));
%         end
%         
        function t = Tensor(mpo)
            assert(depth(mpo) == 1, 'not implemented for 1 < depth');
            W = length(mpo);
            if W == 0
                t = contracttransfer(mpo.L, mpo.R);
            else
                t = contractmpo(mpo.O{:}, mpo.L, mpo.R);
            end
        end
        
    end
    
    methods (Static)
        function mpo = randnc(pspaces, vspaces)
            assert(length(pspaces) == length(vspaces) + 1);
            
            L = MpsTensor(Tensor.randnc([pspaces(1) vspaces(1)'], pspaces(1)));
            O = cell(1, length(pspaces)-2);
            for i = 1:length(O)
                O{i} = MpoTensor(Tensor.randnc([vspaces(i) pspaces(i+1)], ...
                    [pspaces(i+1) vspaces(i+1)]));
            end
            R = MpsTensor(Tensor.randnc([pspaces(end)' vspaces(end)], pspaces(end)'));
            mpo = FiniteMpo(L, O, R);
        end
        
        function T = mps_channel_operator(Atop, O, Abot)
            arguments
                Atop    % top mps tensors
                O       % cell of mpotensors (unrotated)
                Abot    % bottom mps tensors (unconjugated)
            end
            
            assert(length(Atop) == length(O) && length(O) == length(Abot), ...
                'FiniteMpo:argerror', 'unmatching sizes');
            
            for i = numel(Atop):-1:1
                atop = Atop{i};
                o = rot90(O{i});
                twistinds = find(isdual(space(atop, 2:(nspaces(atop) - 1 - atop.alegs))));
                abot = twist(Abot{i}, twistinds + 1)'; % IMPORTANT: ' needs to be outside of twist
                
                assert(isequal(pspace(abot)', leftvspace(o)));
                assert(isequal(rightvspace(o)', pspace(atop)));
                
                T(i) = FiniteMpo(abot, {o}, atop);
            end
        end
    end
end

