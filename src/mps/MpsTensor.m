classdef (InferiorClasses = {?Tensor, ?SparseTensor}) MpsTensor < AbstractTensor
    % Generic mps tensor objects that have a notion of virtual, physical and auxiliary legs.
    
    properties
        var
        plegs
        alegs = 0
    end
    
    
    %% Constructors
    methods
        function A = MpsTensor(tensor, alegs)
            arguments
                tensor = []
                alegs = 0
            end
            
            if iscell(tensor)
                for i = length(tensor):-1:1
                    A(i) = MpsTensor(tensor{i}, alegs);
                end
                return
            end
            
            % TODO: copy constructor
            
            if ~isempty(tensor)
                A.var = tensor;
                A.plegs = nspaces(tensor) - alegs - 2;
                A.alegs = alegs;
            end
        end
    end
    
    methods (Static)
        function A = new(varargin, kwargs)
            arguments (Repeating)
                varargin
            end
            arguments
                kwargs.alegs = 0
            end
            
            A = MpsTensor(Tensor.new(varargin{:}), kwargs.alegs);
        end
        
        function A = randnc(varargin, kwargs)
            arguments (Repeating)
                varargin
            end
            arguments
                kwargs.alegs = 0
            end
            
            A = MpsTensor(Tensor.randnc(varargin{:}), kwargs.alegs);
        end
    end
    
    
    %% Properties
    methods
%         function varargout = size(A, varargin)
%             [varargout{1:nargout}] = size(A.var, varargin{:});
%         end
        
%         function n = numel(A)
%             n = numel(A.var);
%         end
        
%         function l = length(A)
%             l = length(A.var);
%         end
        
        function A = cat(dim, varargin)
            ismpstensor = cellfun(@(x) isa(x, 'MpsTensor'), varargin);
            i = find(ismpstensor, 1);
            A = varargin{i};
            for j = 1:i-1
                B = varargin{j};
                if ismpstensor(j)
                    assert(B.alegs == A.alegs && B.plegs == A.plegs);
                    A.var = cat(dim, B.var, A.var);
                else
                    A.var = cat(dim, B, A.var);
                end
            end
        end
        
        function bool = isscalar(~)
            bool = false;
        end
        
        function A = horzcat(varargin)
            A = cat(2, varargin{:});
        end
        
        function A = vertcat(varargin)
            A = cat(1, varargin{:});
        end
        
        function n = nnz(A)
            n = nnz(A.var);
        end
        
        function A = slice(A, varargin)
            s = substruct('()', varargin);
            A.var = subsref(A.var, s);
        end
        
        function tdst = insert_onespace(tsrc, varargin)
            % insert a trivial space at position i.
            tdst = MpsTensor(insert_onespace(tsrc.var, varargin{:}), tsrc.alegs + 1);
        end
    end
    
    %% Spaces
    methods
        function s = space(A, varargin)
            s = space(A.var, varargin{:});
        end
        
        function n = nspaces(A)
            n = nspaces(A.var);
        end
        
        function cod = codomain(A)
            cod = A.var.codomain;
        end
        
        function dom = domain(A)
            dom = A.var.domain;
        end
        
        function r = rank(A)
            r = rank(A.var);
        end
        
        function s = pspace(A)
            % The physical space of an :class:`MpsTensor`.
            s = space(A, 1 + (1:A.plegs));
        end
        
        function s = leftvspace(A)
            % The left virtual space of an :class:`MpsTensor`.
            s = space(A.var(1), 1);
        end
        
        function s = rightvspace(A)
            % The right virtual space of an :class:`MpsTensor`.
            s = space(A.var(1), nspaces(A.var(1)) - A.alegs);
        end
    end
    
    
    
    %% Linear Algebra
    methods
        
        function varargout = matrixblocks(t)
            [varargout{1:nargout}] = matrixblocks(t.var);
        end
        
        function d = dot(A, B)
            if isa(A, 'MpsTensor')
                A = A.var;
            end
            if isa(B, 'MpsTensor')
                B = B.var;
            end
            d = dot(A, B);
        end
        
        function d = overlap(A, B)
            arguments
                A
                B MpsTensor
            end
            
            inds = 1:nspaces(A);
            
            d = contract(A, inds, B, [fliplr(inds(1:end-B.alegs)) inds(end-B.alegs+1:end)]);
        end
        
        function t = mrdivide(t1, t2)
            if isa(t1, 'MpsTensor')
                t = t1;
                t1 = t1.var;
            else
                t = t2;
            end
            if isa(t2, 'MpsTensor')
                t2 = t2.var;
            end            
            
            t.var = t1 * inv(t2);
        end
        
        function C = mtimes(A, B)
            
            if isscalar(A) || isscalar(B)
                C = A .* B;
                return
            end
            
            szA = size(A);
            szB = size(B);
            if szA(2) ~= szB(1)
                error('mtimes:dimagree', ...
                    'incompatible dimensions (%d) (%d)', szA(2), szB(1));
            end

            if ~(isnumeric(A) || isnumeric(B))
                error('not implemented when neither input is numeric.');
            end
            
            % TODO: fairly hacky for now but it works; should be revisited properly
            
            if isnumeric(A)
                al = repmat([B(1, :).alegs], szA(1), 1);
                B = reshape([B.var], szB);
            else
                al = repmat([A(:, 1).alegs], 1, szB(2));
                A = reshape([A.var], szA);
            end
            Cv = A * B;
            
            for i = szA(1):-1:1
                for j = szB(2):-1:1
                    C(i, j) = MpsTensor(Cv(i, j), al(i, j));
                end
            end


%             if isnumeric(A)
%                 al = repmat([B(1, :).alegs], szA(1), 1);
%                 %B = reshape([B.var], szB);
%             else
%                 al = repmat([A(:, 1).alegs], 1, szB(2));
%                 %A = reshape([A.var], szA); % this does not zork becquse it qlso needs to hqndle the cqse zhen the vqrs the,selves qre qlreqdy qrrqys; you should reqlly do the ,qnuql loop
%             end
%             
%             for i = szA(1):-1:1
%                 for j = szB(2):-1:1
%                     C(i, j) = MpsTensor(A(i, 1) * B(1, j), al(i, j));
%                     for k = 2:szA(2)
%                         C(i, j) = C(i, j) + ...
%                              MpsTensor(A(i, k).var * B(k, j), al(i, j));
%                     end
%                 end
%             end

        end
        
        function C = times(A, B)
            if isscalar(A)
                C = B;
                for i = 1:numel(C)
                    C(i).var = A .* B(i).var;
                end
                return
            end
            if isscalar(B)
                C = A;
                for i = 1:numel(C)
                    C(i).var = A(i).var .* B;
                end
                return
            end
            
            error('not implemented when neither input is scalar.');
        end
        
        function A = repartition(A, varargin)
            A.var = repartition(A.var, varargin{:});
        end
        
        function A = plus(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    alegs = varargin{i}.alegs;
                    varargin{i} = varargin{i}.var;
                end
            end
            A = MpsTensor(plus(varargin{:}), alegs);
        end
        
        function A = minus(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    alegs = varargin{i}.alegs;
                    varargin{i} = varargin{i}.var;
                end
            end
            A = MpsTensor(minus(varargin{:}), alegs);
        end
        
        function A = rdivide(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    alegs = varargin{i}.alegs;
                    varargin{i} = varargin{i}.var;
                end
            end
            A = MpsTensor(rdivide(varargin{:}), alegs);
        end
        
        function A = ldivide(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    alegs = varargin{i}.alegs;
                    varargin{i} = varargin{i}.var;
                end
            end
            A = MpsTensor(ldivide(varargin{:}), alegs);
        end
        
        function n = norm(A, varargin)
            n = norm(A.var, varargin{:});
        end
        
        function d = distance(A, B)
            d = norm(A - B);
        end

        
        function [A, n] = normalize(A)
            [A.var, n] = normalize(A.var);
        end
        
        function A = conj(A)
            A.var = conj(A.var);
        end

        function A = twist(A, varargin)
            A.var = twist(A.var, varargin{:});
        end
        
        function A = twistdual(A, varargin)
            A.var = twistdual(A.var, varargin{:});
        end
        
        function t = ctranspose(t)
            % Compute the adjoint of a tensor. This is defined as swapping the codomain and
            % domain, while computing the adjoint of the matrix blocks.
            %
            % Usage
            % -----
            % :code:`t = ctranspose(t)`
            % :code:`t = t'`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   adjoint tensor.
            
            for i = 1:numel(t)
                t(i).var = t(i).var';
            end
            
            t = permute(t, ndims(t):-1:1);
        end
        
        function A = tpermute(A, varargin)
            A.var = tpermute(A.var, varargin{:});
        end
        
        function C = tensorprod(varargin)
            for i = 1:2
                if isa(varargin{i}, 'MpsTensor')
                    varargin{i} = varargin{i}.var;
                end
            end
            C = tensorprod(varargin{:});
        end
        
        function T = transfermatrix(A, B)
            arguments
                A MpsTensor
                B MpsTensor = A
            end
            
            B = twist(B, [isdual(space(B, 1:B.plegs+1)) ~isdual(space(B, B.plegs+2))]);
            T = FiniteMpo(B', {}, A);
        end
        
        function C = initializeC(A)
            % Initialize a set of gauge tensors for the given mpstensors.
            arguments (Repeating)
                A MpsTensor
            end
            C = cell(size(A));
            
            for i = length(A):-1:1
                C{i} = A{i}.var.eye(rightvspace(A{i})', leftvspace(A{next(i, length(A))}));
            end
        end
        
        function A = expand(A, leftvspace, rightvspace, noisefactor)
            % Expand a tensor to the given virtual spaces.
            arguments
                A
                leftvspace
                rightvspace
                noisefactor = 1e-3
            end
            
            spaces = space(A);
            spaces(nspaces(A) - A.alegs) = rightvspace;
            spaces(1) = conj(leftvspace);
            r = rank(A);
            A.var = embed(A.var, ...
                noisefactor * ...
                normalize(A.var.randnc(spaces(1:r(1)), spaces(r(1)+1:end)')));
        end
        
        function type = underlyingType(A)
            type = underlyingType(A.var);
        end
        
        function disp(t)
            fprintf('%s with %d plegs and %d alegs:\n', class(t), t.plegs, t.alegs);
            disp(t.var);
        end
    end
    
    
    %% Factorizations
    methods
        function [A, L] = leftorth(A, alg)
            arguments
                A
                alg = 'polar'
            end
            
            if A.alegs == 0
                [A.var, L] = leftorth(A.var, 1:nspaces(A)-1, nspaces(A), alg);
                if isdual(space(L, 1)) == isdual(space(L, 2))
                    L.codomain = conj(L.codomain);
                    L = twist(L, 1);
                    A.var.domain = conj(A.var.domain);
                end
            else
                [A.var, L] = leftorth(A.var, [1:A.plegs+1 A.plegs+3], A.plegs+2, alg);
                A.var = permute(A.var, [1:A.plegs+1 A.plegs+3 A.plegs+2], rank(A));
            end
        end
        
        function [R, A] = rightorth(A, alg)
            arguments
                A
                alg = 'rqpos'
            end
            
            [R, A.var] = rightorth(A.var, 1, 2:nspaces(A), alg);
            if isdual(space(R, 1)) == isdual(space(R, 2))
                R.domain = conj(R.domain);
                R = twist(R, 2);
                A.var.codomain = conj(A.var.codomain);
            end
        end
        
        function varargout = tsvd(t, varargin)
            [varargout{1:nargout}] = tsvd(t.var, varargin{:});
        end
        
        function A = leftnull(A, alg)
            arguments
                A
                alg = 'svd'
            end
            
            if numel(A) > 1
                for i = numel(A):-1:1
                    A(i) = leftnull(A(i), alg);
                end
                return
            end
            
            p = 1:nspaces(A);
            p1 = p(1:(end - A.alegs - 1));
            p2 = p((end - A.alegs):end);
            A.var = leftnull(A.var, p1, p2, alg);
        end
        
        function A = rightnull(A, alg)
            arguments
                A
                alg = 'svd'
            end
            
            if numel(A) > 1
                for i = numel(A):-1:1
                    A(i) = rightnull(A(i), alg);
                end
                return
            end
            
            A.var = rightnull(A.var, 1, 2:nspaces(A), alg);
        end
    end
    
    %% Contractions
    methods
        function v = applytransfer(L, R, v)
            arguments
                L MpsTensor
                R MpsTensor
                v = []
            end
            
            if isempty(v)
                v = tracetransfer(L, R);
                return
            end
            
            auxlegs_v = nspaces(v) - 2;
            auxlegs_l = L.alegs;
            auxlegs_r = R.alegs;
            newrank = rank(v); newrank(2) = newrank(2) + auxlegs_l + auxlegs_r;
            
            v = contract( ...
                v, [1, R.plegs+2, (-(1:auxlegs_v) - 2 - auxlegs_l)], ...
                L, [-1, (1:L.plegs)+1, 1 (-(1:auxlegs_l) - (1+L.plegs))], ...
                R, [R.plegs+2, (R.plegs:-1:1)+1, -2, (-(1:auxlegs_r) - (2+R.plegs) - auxlegs_l - auxlegs_v)], ...
                'Rank', newrank);
        end
        
        function v = contracttransfer(L, R)
            arguments
                L MpsTensor
                R MpsTensor
            end
            
            auxlegs_l = L.alegs;
            auxlegs_r = R.alegs;
            assert(R.plegs == L.plegs);
            plegs = L.plegs; %#ok<PROPLC>
            
            v = contract(...
                L, [-1, 1:plegs, -4, (-(1:auxlegs_l) - 4)], ...
                R, [-3, flip(1:plegs), -2 (-(1:auxlegs_r) - 4 - auxlegs_l)], ...
                'Rank', [2, 2] + [0, auxlegs_l + auxlegs_r]); %#ok<PROPLC>
        end
        
        function v = tracetransfer(L, R)
            arguments
                L MpsTensor
                R MpsTensor
            end
            
            auxlegs_l = L.alegs;
            auxlegs_r = R.alegs;
            assert(R.plegs == L.plegs);
            plegs = L.plegs; %#ok<PROPLC>
            newrank = [1 1 + auxlegs_l + auxlegs_r];
            
            v = contract(...
                L, [-1 1:(plegs + 1) (-(1:auxlegs_l) - 2)], ...
                R, [flip(1:(plegs + 1)) -2 (-(1:auxlegs_r) - 2 - auxlegs_l)], ...
                'Rank', newrank); %#ok<PROPLC>
        end
        
        function A = multiplyleft(A, C)
            % Multiply a gauge matrix from the left.
            assert(nspaces(C) == 2, 'mps:tba', 'not implemented yet');
            A.var = contract(C, [-1 1], ...
                A.var, [1 -((1:A.plegs)+1) -((1:A.alegs+1)+A.plegs+1)], 'Rank', rank(A));
        end
        
        function A = multiplyright(A, C)
            % Multiply a gauge matrix from the right.
            r = rank(A);    npA = A.plegs;  naA = A.alegs; 
            nC = nspaces(C);    naC = nC - 2;
            r(2) = r(2) + naC;
            A.var = contract(A.var, [-(1:npA + 1) 1 -((1:naA) + npA + 2)], ...
                C, [1 -(2 + npA + (1:(nC - 1)))], 'Rank', r);
            A.alegs = A.alegs + naC;
        end
    end
    
    %%
    methods
        function bool = isisometry(t, varargin)
            bool = isisometry(t.var, varargin{:});
        end
    end
    
    
    %% Converters
    methods
        function t = Tensor(A)
            t = full(A.var);
        end
        
        function t = SparseTensor(A)
            t = reshape([A.var], size(A));
            t = sparse(t);
        end
    end
    
    
    %% Solvers
    methods
        function v = vectorize(t, type)
            % Collect all parameters in a vector, weighted to reproduce the correct
            % inproduct.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % type : 'real' or 'complex'
            %   optionally specify if complex entries should be seen as 1 or 2 parameters.
            %   Defaults to 'complex', with complex parameters.
            %
            % Returns
            % -------
            % v : numeric
            %   real or complex vector containing the parameters of the tensor.
            
            arguments
                t
                type = 'complex'
            end
            
            v = vectorize(t.var, type);
        end
        
        function t = devectorize(v, t, type)
            % Collect all parameters from a vector, and insert into a tensor.
            %
            % Arguments
            % ---------
            % v : numeric
            %   real or complex vector containing the parameters of the tensor.
            %
            % t : :class:`Tensor`
            %   input tensor.
            %
            % type : 'real' or 'complex'
            %   optionally specify if complex entries should be seen as 1 or 2 parameters.
            %   Defaults to 'complex', with complex parameters.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor, filled with the parameters.
            
            arguments
                v
                t
                type = 'complex'
            end
            
            t.var = devectorize(v, t.var, type);
        end
    end
    
    
    %%
    methods (Static)
        function local_tensors = decompose_local_state(psi, kwargs)
            % convert a tensor into a product of local operators.
            %
            % Usage
            % -----
            % :code:`local_operators = MpoTensor.decompose_local_operator(H, kwargs)`.
            %
            % Arguments
            % ---------
            % H : :class:`AbstractTensor`
            %   tensor representing a local operator on N sites.
            %
            % Keyword Arguments
            % -----------------
            % 'Trunc' : cell
            %   optional truncation method for the decomposition. See also
            %   :meth:`Tensor.tsvd`
            arguments
                psi
                kwargs.Trunc = {'TruncBelow', 1e-14}
            end
            
            L = nspaces(psi);
            assert(L >= 3, 'argerror', ...
                sprintf('state must have at least 3 legs. (%d)', L));
            
            local_tensors = cell(1, L - 2);
            for i = 1:length(local_tensors)-1
                [u, s, psi] = tsvd(psi, [1 2], [3:nspaces(psi)], kwargs.Trunc{:});
                local_tensors{i} = multiplyright(MpsTensor(u), s);
            end
            local_tensors{end} = MpsTensor(repartition(psi, [2 1]));
        end
    end
    
    methods (Hidden)
        function tdst = zerosLike(t, varargin)
            tdst = repmat(0 * t, varargin{:});
        end
    end

end

