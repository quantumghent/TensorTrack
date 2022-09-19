classdef Tensor
    %TENSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        codomain
        domain
    end
    
    properties
        var
    end
    
    %% Constructors
    methods
        function t = Tensor(varargin)
            % Create a tensor object.
            %
            % Usage
            % -----
            % :code:`t = Tensor(array)`
            %
            % :code:`t = Tensor(codomain, domain)`
            %
            % Arguments
            % ---------
            % array : numeric
            %   numeric input array to convert to a :class:`Tensor`
            %
            % codomain, domain : :class:`AbstractSpace`
            %   spaces that define the structure of the output tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   new tensor object.
            
            if nargin == 0, return; end
            
            % t = Tensor(tensor)
            if nargin == 1 && isa(varargin{1}, 'Tensor')
                t.codomain = varargin{1}.codomain;
                t.domain = varargin{1}.domain;
                t.var = varargin{1}.var;
                return
            end
            
            % t = Tensor(array)
            if nargin == 1 && isnumeric(varargin{1})
                t.domain = [];
                t.codomain = CartesianSpace.new(size(varargin{1}));
                t.var = fill_tensor_data(AbstractBlock.new(t.codomain, t.domain), ...
                    varargin);
                return;
            end
                
            % t = Tensor(codomain, domain)
            if nargin == 2
                codomain = varargin{1};
                domain = varargin{2};
                
                assert(~isempty(codomain) || ~isempty(domain), ...
                    'Cannot create (0,0) tensors.');
                t.domain    = domain;
                t.codomain  = codomain;
                
                t.var = AbstractBlock.new(codomain, domain);
                return
            end
            
            error('tensors:SyntaxError', 'Invalid constructor syntax.');
        end
        
        function t = fill_matrix(t, data, charges)
            % Fill the matrix blocks of a tensor.
            %
            % Usage
            % -----
            % :code:`t = fill_matrix(t, matrices, charges)`
            % 
            % :code:`t = fill_matrix(t, fun, charges)`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor to fill into.
            %
            % matrices : cell or numeric
            %   list of matrices or single matrix to fill with.
            %
            % fun : :class:`function_handle`
            %   function of signature :code:`fun(dims, charge)` to fill with.
            % 
            % charges : :class:`AbstractCharge`
            %   optional list of charges to identify the matrix blocks.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   filled tensor.
            
            arguments
                t
                data
                charges = [];
            end
            
            if isnumeric(data), data = {data}; end
            
            if iscell(data)
                t.var = fill_matrix_data(t.var, data, charges);
            else
                t.var = fill_matrix_fun(t.var, data, charges);
            end
        end
        
        function t = fill_tensor(t, data, trees)
            % Fill the tensor blocks of a tensor.
            %
            % Usage
            % -----
            % :code:`t = fill_tensor(t, tensors, trees)`
            % 
            % :code:`t = fill_tensor(t, fun, trees)`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor to fill into.
            %
            % tensors : cell or numeric
            %   list of tensors or single tensor to fill with.
            %
            % fun : :class:`function_handle`
            %   function of signature :code:`fun(dims, trees)` to fill with.
            % 
            % trees : :class:`FusionTree`
            %   optional list of fusion trees to identify the tensor blocks.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   filled tensor.
            
            arguments
                t
                data
                trees = []
            end
            
            if isnumeric(data), data = {data}; end
            
            if iscell(data)
                t.var = fill_tensor_data(t.var, data, trees);
            else
                t.var = fill_tensor_fun(t.var, data, trees);
            end
        end
        
        function t = similar(fun, tensors, indices, kwargs)
            % Create a tensor based on the indices of other tensors.
            %
            % Usage
            % -----
            % :code:`t = similar(fun, tensors, indices, kwargs)`
            %
            % Arguments
            % ---------
            % fun : :class:`function_handle`
            %   function to fill the tensor data with. This should have signature
            %   based on keyword Mode. If left empty, this defaults to a random complex
            %   tensor.
            %
            % Repeating Arguments
            % -------------------
            % tensors : :class:`Tensor`
            %   input tensors used to copy legs.
            %
            % indices : int
            %   array of which indices to copy for each input tensor.
            %
            % Keyword Arguments
            % -----------------
            % Rank : (1, 2) int
            %   rank of the output tensor, by default this is :code:`[nspaces(t) 0]`.
            %
            % Conj : logical
            %   flag to indicate whether the space should be equal to the input space, or
            %   fit onto the input space. This can be either an array of size(tensors), or a
            %   scalar, in which case it applies to all tensors.
            %
            % Mode : 'tensor' or 'matrix'
            %   method of filling the tensor data. By default this is matrix, where the
            %   function should be of signature :code:`fun(dims, charge)`, for 'tensor' this
            %   should be of signature :code:`fun(dims, tree)`.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor.
            %
            % Examples
            % --------
            % :code:`t = similar([], mpsbar, 1, mpo, 4, mps, 1, 'Conj', true)` creates a 
            % left mps environment tensor.
            
            arguments
                fun = []
            end
            
            arguments (Repeating)
                tensors
                indices
            end
            
            arguments
                kwargs.Rank (1, 2) = [sum(cellfun(@length, indices)) 0]
                kwargs.Conj = false(size(tensors))
                kwargs.Mode = 'matrix'
            end
            
            if isempty(fun)
                fun = @(dims, charge) randnc(dims);
            end
            
            ntotal = sum(cellfun('length', indices));
            ctr = 0;
            for i = length(indices):-1:1
                ctr = ctr + length(indices{i});
                if kwargs.Conj(i)
                    sp(ntotal - ctr + (1:length(indices{i}))) = ...
                        conj(space(tensors{i}, indices{i}));
                else
                    sp(ntotal - ctr + (1:length(indices{i}))) = ...
                        space(tensors{i}, indices{i});
                end
            end
            
            assert(sum(kwargs.Rank) == ntotal, 'tensors:ArgumentError', ...
                'Invalid rank specified.');
            t = Tensor.new(fun, sp(1:kwargs.Rank(1)), ...
                sp((1:kwargs.Rank(2)) + kwargs.Rank(1))', 'Mode', kwargs.Mode);
        end
    end
    
    methods (Static)
        function t = new(fun, varargin, kwargs)
            % Create a tensor with data using a function handle.
            %
            % Usage
            % -----
            % :code:`Tensor.new(fun, dims)`
            %
            % :code:`Tensor.new(fun, dims, arrows)`
            %
            % :code:`Tensor.new(fun, codomain, domain)`
            %
            % :code:`Tensor.new(fun, tensor)`
            %
            % :code:`Tensor.new(..., 'Rank', r, 'Mode', mode)`
            %
            % Arguments
            % ---------
            % fun : :class:`function_handle`
            %   function of signature :code:`fun(dims, id)` where id is determined by Mode.
            %   If this is left empty, the tensor data will be uninitialized.
            %
            % dims : int
            %   list of dimensions for non-symmetric tensors.
            %
            % arrows : logical
            %   optional list of arrows for tensor legs.
            %
            % tensor : :class:`Tensor`
            %   input tensor to copy structure.
            %
            % Repeating Arguments
            % -------------------
            % charges : cell
            %   list of charges for each tensor index.
            %
            % degeneracies : cell
            %   list of degeneracies for each tensor index.
            %
            % arrow : logical
            %   arrow for each tensor index.
            %
            % Keyword Arguments
            % -----------------
            % Rank : int (1, 2)
            %   rank of the constructed tensor. By default this is [nspaces(t) 0].
            %
            % Mode : 'matrix' or 'tensor'
            %   method of filling the resulting tensor. When this is 'matrix' (default),
            %   the function signature is :code:`fun(dims, charge)`, while for 'tensor' the
            %   signature should be :code:`fun(dims, tree)`.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor.
            
            arguments
                fun = []
            end
            
            arguments (Repeating)
                varargin
            end
            
            arguments
                kwargs.Rank
                kwargs.Mode = 'matrix'
            end
            
            % Parse special constructors
            if ~isempty(varargin{1}) && isnumeric(varargin{1})
                if length(varargin) == 1
                    spaces = CartesianSpace.new(varargin{1});
                elseif isnumeric(varargin{2}) || islogical(varargin{2})
                    assert(length(varargin{1}) == length(varargin{2}))
                    spaces = ComplexSpace.new(varargin{1}, varargin{2});
                end
                
                if ~isfield(kwargs, 'Rank')
                    kwargs.Rank = [length(spaces) 0];
                else
                    assert(sum(kwargs.Rank) == length(spaces), 'tensors:ArgumentError', ...
                        'User supplied invalid rank.');
                end
                
                codomain = spaces(1:kwargs.Rank(1));
                domain = spaces(kwargs.Rank(1) + (1:kwargs.Rank(2)))';
                
            elseif length(varargin) == 1 && isa(varargin{1}, 'Tensor')
                codomain = varargin{1}.codomain;
                domain = varargin{1}.domain;
                
                if isfield(kwargs, 'Rank')
                    assert(isequal(kwargs.Rank, [length(codomain) length(domain)]), ...
                        'tensors:ArgumentError', 'Rank and spaces incompatible.');
                end
                
            elseif isa(varargin{1}, 'AbstractSpace') || isa(varargin{2}, 'AbstractSpace')
                codomain = varargin{1};
                domain = varargin{2};
                
                if isfield(kwargs, 'Rank')
                    assert(isequal(kwargs.Rank, [length(codomain) length(domain)]), ...
                        'tensors:ArgumentError', 'Rank and spaces incompatible.');
                end
            end
            
            t = Tensor(codomain, domain);
            
            % Fill tensor
            if ~isempty(fun)
                if strcmp(kwargs.Mode, 'matrix')
                    t = fill_matrix(t, fun);
                elseif strcmp(kwargs.Mode, 'tensor')
                    t = fill_tensor(t, fun);
                else
                    error('tensors:ArgumentError', 'Unknown options for Mode');
                end
            end
        end

        function t = zeros(varargin)
            t = Tensor.new(@(dims, charge) zeros(dims), varargin{:});
        end
        
        function t = ones(varargin)
            t = Tensor.new(@(dims, charge) ones(dims), varargin{:});
        end
        
        function t = eye(varargin)
            t = Tensor.new(@(dims, charge) eye(dims), varargin{:});
        end
        
        function t = rand(varargin)
            t = Tensor.new(@(dims, charge) rand(dims), varargin{:});
        end
        
        function t = randc(varargin)
            t = Tensor.new(@(dims, charge) randc(dims), varargin{:});
        end
        
        function t = randnc(varargin)
            t = Tensor.new(@(dims, charge) randnc(dims), varargin{:});
        end
    end
    
    
    %% Structure
    methods
        function n = indin(t)
            n = length(t.domain);
        end
        
        function n = indout(t)
            n = length(t.codomain);
        end
        
        function n = nspaces(t)
            n = length(t.domain) + length(t.codomain);
        end
        
        function r = rank(t, i)
            r = [length(t.codomain) length(t.domain)];
            if nargin > 1
                r = r(i);
            end
        end
        
        function sz = dims(t, inds)
            sz = dims([t.codomain, t.domain']);
            if nargin > 1
                sz = sz(inds);
            end
        end
        
        function sp = space(t, inds)
            sp = [t.codomain t.domain'];
            if nargin > 1
                sp = sp(inds);
            end
        end
        
        function f = fusiontrees(t)
            f = fusiontrees(t.codomain, t.domain);
        end
        
        function b = tensorblocks(t)
            b = tensorblocks(t.var);
        end
        
        function varargout = matrixblocks(t)
            [varargout{1:nargout}] = matrixblocks(t.var);
        end
    end
    
    
    %% Comparison
    methods
        function d = distance(A, B)
            % Compute the Euclidean distance between two tensors.
            %
            % Arguments
            % ---------
            % A, B : :class:`Tensor`
            %
            % Returns
            % -------
            % d : numeric
            %   Euclidean distance, defined as the norm of the distance.
            
            assert(isequal(size(A), size(B)) || isscalar(A) || isscalar(B), ...
                'tensors:SizeError', 'Incompatible sizes for vectorized function.');
            
            % make everything a vector
            A = arrayfun(@repartition, A);
            B = arrayfun(@repartition, B);
            
            d = norm(A - B);
        end
    end
        
    
    %% Linear algebra
    methods
        function t = conj(t)
            % Compute the conjugate of a tensor. This is defined as taking an element-wise
            % conjugate, but implemented as taking the adjoint and reversing the order of
            % the indices.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   conjugate tensor.
            
            t = permute(t', nspaces(t):-1:1, rank(t));
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
            
            [t.codomain, t.domain] = swapvars(t.codomain, t.domain);
            t.var = t.var';
        end
        
        function d = dot(t1, t2)
            % Compute the scalar dot product of two tensors. This is defined as the overlap 
            % of the two tensors, which therefore must have equal domain and codomain. This
            % function is sesquilinear in its arguments.
            %
            % Arguments
            % ---------
            % t1, t2 : :class:`Tensor`
            %   tensors of equal structure.
            %
            % Returns
            % -------
            % d : double
            %   scalar dot product of the two tensors.
            
            assert(isequal(t1.domain, t2.domain) && isequal(t1.codomain, t2.codomain), ...
                'tensors:SpaceMismatch', ...
                'Dot product only defined for tensors of equal structure.');

            d = 0;
            [mblocks1, mcharges] = matrixblocks(t1.var);
            mblocks2 = matrixblocks(t2.var);
            qdims = qdim(mcharges);
            for i = 1:length(t1.var)
                d = d + qdims(i) * sum(conj(mblocks1{i}) .* mblocks2{i}, 'all');
            end
        end
        
        function t1 = minus(t1, t2)
            % Compute the difference between two tensors.
            %
            % Arguments
            % ---------
            % t1, t2 : :class:`Tensor` or numeric
            %   input tensors, scalars are interpreted as scalar * eye.
            %
            % Returns
            % -------
            % t1 : :class:`Tensor`
            %   output tensor
            
            if isnumeric(t1), t1 = t1 + (-t2); return; end
            assert(isequal(size(t1), size(t2)), 'Incompatible sizes for vectorized minus.');
            
            for i = 1:numel(t1)
                if isnumeric(t2(i))
                    assert(isequal(t1(i).codomain, t1(i).domain), 'tensors:NonSquare', ...
                        'Subtraction with scalars only defined for square tensors.');
                    I = t1(i).eye(t1(i).codomain, t1(i).domain);
                    t1(i).var = t1(i).var - I.var .* t2(i);
                else
                    assert(isequal(t1(i).domain, t2(i).domain) && ...
                        isequal(t1(i).codomain, t2(i).codomain), 'tensors:SpaceMismatch', ...
                        'Cannot subtract tensors of different structures.');
                    t1(i).var = t1(i).var - t2(i).var;
                end
            end
        end
        
        function t = mldivide(t1, t2)
            % Left division of tensors.
            %
            % Usage
            % -----
            % :code:`t = mldivide(t1, t2)`
            %
            % :code:`t = t1 \ t2`
            %
            % Arguments
            % ---------
            % t1, t2 : :class:`Tensor` or numeric
            %   input tensor or scalar.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor.
            
            t = inv(t1) * t2;
        end
            
        
        function t = mrdivide(t1, t2)
            % Right division of tensors.
            %
            % Usage
            % -----
            % :code:`t = mrdivide(t1, t2)`
            %
            % :code:`t1 = t1 / t2`
            %
            % Arguments
            % ---------
            % t1, t2 : :class:`Tensor` or numeric
            %   input tensor or scalar.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor.
            
            t = t1 * inv(t2);
        end
        
        function C = mtimes(A, B)
            % Compose tensors, interpreted as maps from domain to codomain.
            %
            % Usage
            % -----
            % :code:`A = mtimes(A, B)`
            %
            % :code:`A = A * B`
            %
            % Arguments
            % ---------
            % A, B : :class:`Tensor`
            %   input tensors, satisfying A.domain = B.codomain.
            %
            % Returns
            % -------
            % C : :class:`Tensor`
            %   output tensor.
            
            szA = size(A);
            szB = size(B);
            assert(szA(2) == szB(1));
            
            % Scalar multiplications
            if isnumeric(A) || isnumeric(B)
                for i = szA(1):-1:1
                    for j = szB(2):-1:1
                        C(i,j) = A(i, 1) .* B(1, j);
                        for k = 2:szA(2)
                            C(i,j) = C(i,j) + A(i, k) .* B(k, j);
                        end
                    end
                end
                return
            end
            
            % Tensor multiplications
            for i = szA(1):-1:1
                for j = szB(2):-1:1
                    C(i,j) = tensorprod(A(i, 1), B(1, j), ...
                        rank(A(i, 1), 1) + (1:rank(A(i, 1), 2)), ...
                        rank(B(1, j), 1):-1:1);
                    for k = 2:szA(2)
                        C(i,j) = C(i,j) + tensorprod(A(i, k), B(k, j), ...
                            rank(A(i, k), 1) + (1:rank(A(i, k), 2)), ...
                            rank(B(k, j), 1):-1:1);
                    end
                end
            end
        end
        
        function n = norm(t, p)
            % Compute the matrix p-norm of a tensor.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor, considered as a matrix from domain to codomain.
            %
            % p : 1, 2, inf or 'fro'
            %   type of norm to compute
            %
            % Returns
            % -------
            % n : :class:`double`
            %   matrix norm of the tensor.
            
            arguments
                t
                p = 'fro'
            end
            
            switch p
                case {1, 2, Inf}
                    n = 0;
                    for i = 1:numel(t)
                        mblocks = matrixblocks(t(i).var);
                        for j = 1:length(mblocks)
                            n = max(norm(mblocks{j}, p), n);
                        end
                    end
                    
                case 'fro'
                    n = 0;
                    for i = 1:numel(t)
                        [mblocks, mcharges] = matrixblocks(t(i).var);
                        qdims = qdim(mcharges);
                        for j = 1:length(mblocks)
                            n = n + qdims(j) * sum(conj(mblocks{j}) .* mblocks{j}, 'all');
                        end
                    end
                    n = sqrt(n);
                    
                otherwise
                    error('tensors:ArgumentError', 'Invalid norm selected.');
            end
        end
        
        function t = normalize(t)
            for i = 1:numel(t)
                t(i) = t(i) .* (1 / norm(t(i)));
            end
        end
        
        function t1 = plus(t1, t2)
            if isnumeric(t1), t1 = t2 + t1; return; end
            assert(isequal(size(t1), size(t2)), 'Incompatible sizes for vectorized plus.');
            
            for i = 1:numel(t1)
                if isnumeric(t2(i))
                    assert(isequal(t1(i).codomain, t1(i).domain), 'tensors:NonSquare', ...
                        'Addition with scalars only defined for square tensors.');
                    I = t1(i).eye(t1(i).codomain, t1(i).domain);
                    t1(i).var = t1(i).var + I.var .* t2(i);
                else
                    assert(isequal(t1(i).domain, t2(i).domain) && ...
                        isequal(t1(i).codomain, t2(i).codomain), 'tensors:SpaceMismatch', ...
                        'Cannot add tensors of different structures.');
                    t1(i).var = t1(i).var + t2(i).var;
                end
            end
        end
        
        function t = permute(t, p, r)
            % Permute the spaces of a tensor.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % p : (1, :) int
            %   permutation vector, by default a trivial permutation.
            %
            % r : (1, 2) int
            %   rank of the output tensor, by default equal to the rank of the input tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   permuted tensor with desired rank.
            
            arguments
                t
                p = 1:nspaces(t)
                r = rank(t)
            end
            
            if isempty(p), p = 1:nspaces(t); end
            if isempty(r), r = rank(t); end
            
            assert(length(p) == nspaces(t), 'Invalid permutation.');
            assert(length(p) == sum(r), 'Invalid new rank.');
            
            if (all(p == 1:nspaces(t)) && all(rank(t) == r)), return; end
            
            persistent cache
            if isempty(cache), cache = LRU; end
            
            if Options.CacheEnabled()
                key = GetMD5({GetMD5_helper(t.codomain), GetMD5_helper(t.domain), p, r}, ...
                    'Array', 'hex');
                med = get(cache, key);
                if isempty(med)
                    med = struct;
                    med.structure = similar(@(x,charge) uninit(x), t, p, 'Rank', r);
                    med.map = permute(fusiontrees(t), p, r);
                    cache = set(cache, key, med);
                end
                
                tdst = med.structure;
                map = med.map;
            else
                tdst = similar(@(x,charge) uninit(x), t, p, 'Rank', r);
                map = permute(fusiontrees(t), p, r);
            end
            
            tdst.var = axpby(1, t.var, 0, tdst.var, p, map);
            t = tdst;
        end
        
        function t = repartition(t, r)
            % Change the rank of a tensor.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % r : (1, 2) int
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   repartitioned tensor with desired rank.
            
            arguments
                t
                r (1,2) = [nspaces(t) 0]
            end
            
            assert(sum(r) == sum(rank(t)), 'tensors:ValueError', 'Invalid new rank.');
            t = permute(t, 1:nspaces(t), r);
        end
        
        function t = rdivide(t, a)
            % Scalar division of a tensor and a scalar.
            %
            % Usage
            % -----
            % :code:`t = rdivide(t, a)`
            %
            % :code:`t = t ./ a`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            % uminus
            % a : numeric
            %   input scalar.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor.
            
            t.var = rdivide(t.var, a);
        end
        
        function C = tensorprod(A, B, dimA, dimB, ca, cb, options)
            % Compute the contraction of two tensors through the selected spaces.
            %
            % Arguments
            % ---------
            % A, B : :class:`Tensor`
            %   input tensors, must satisfy space(A, dimA) = conj(space(B, dimB)).
            %
            % dimA, dimB : (1, :) int
            %   selected indices to contract.
            %
            % Keyword Arguments
            % -----------------
            % NumDimensionsA : int
            %   number of spaces of A, to satisfy builtin tensorprod syntax.
            %
            % Returns
            % -------
            % C : :class:`Tensor` or numeric
            %   output tensor, with the uncontracted spaces of A as codomain, and the
            %   uncontracted spaces of B as domain, or output scalar, if no uncontracted
            %   spaces remain.
            
            arguments
                A
                B
                dimA
                dimB
                ca = false
                cb = false
                options.NumDimensionsA
            end
            
            uncA = 1:nspaces(A);    uncA(dimA) = [];
            iA = [uncA dimA];       rA = [length(uncA) length(dimA)];
            
            uncB = 1:nspaces(B);    uncB(dimB) = [];
            iB = [flip(dimB) uncB]; rB = [length(dimB) length(uncB)];
            
            if ca
                A = A';
                idx = length(iA):-1:1;
                iA = idx(iA);
            end
            
            if cb
                B = B';
                idx = length(iB):-1:1;
                iB = idx(iB);
            end
            
            persistent cache
            if isempty(cache), cache = LRU; end
            
            if Options.CacheEnabled()
                key = GetMD5({GetMD5_helper(A.codomain), GetMD5_helper(A.domain), ...
                    GetMD5_helper(B.codomain), GetMD5_helper(B.domain), ...
                    dimA, dimB, ca, cb}, 'Array', 'hex');
                med = get(cache, key);
                if isempty(med)
                    med = struct;
                    A_ = similar(@(x,charge) uninit(x), A, iA, 'Rank', rA);
                    med.varA = A_.var;
                    med.mapA = permute(fusiontrees(A), iA, rA);
                    
                    B_ = similar(@(x,charge) uninit(x), B, iB, 'Rank', rB);
                    med.varB = B_.var;
                    med.mapB = permute(fusiontrees(B), iB, rB);
                    
                    assert(isequal(A_.domain, B_.codomain), 'tensors:SpaceMismatch', ...
                        'Contracted spaces incompatible.');
                    if ~isempty(A_.codomain) || ~isempty(B_.domain) 
                        med.C = Tensor(A_.codomain, B_.domain);
                    else
                        med.C = [];
                    end
                    cache = set(cache, key, med);
                end
                
                varA = axpby(1, A.var, 0, med.varA, iA, med.mapA);
                varB = axpby(1, B.var, 0, med.varB, iB, med.mapB);
                C = med.C;
                if isempty(C)
                    Ablocks = matrixblocks(varA);
                    Bblocks = matrixblocks(varB);
                    C = horzcat(Ablocks{:}) * vertcat(Bblocks{:});
                else
                    C.var = mul(C.var, varA, varB);
                end
                return
            end
            
            A = permute(A, iA, rA);
            B = permute(B, iB, rB);
            
            assert(isequal(A.domain, B.codomain), 'tensors:SpaceMismatch', ...
                'Contracted spaces incompatible.');
            
            % scalar result
            if isempty(uncA) && isempty(uncB)
                Ablocks = matrixblocks(A);
                Bblocks = matrixblocks(B);
                C = horzcat(Ablocks{:}) * vertcat(Bblocks{:});
                return
            end
            
            % tensor results
            C = Tensor(A.codomain, B.domain);
            C.var = mul(C.var, A.var, B.var);
        end
        
        function t = times(t, a)
            % Scalar product of a tensor and a scalar.
            %
            % Usage
            % -----
            % :code:`t = times(t, a)`
            % :code:`t uminus= t .* a`
            % :code:`t = a .* t`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            % 
            % a : numeric
            %   input scalar.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor.
            
            if isnumeric(t), [t, a] = swapvars(t, a); end
            t.var = times(t.var, a);
        end
        
        function tr = trace(t)
            % Compute the matrix trace of a tensor.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor, considered as a matrix from domain to codomain.
            %
            % Returns
            % -------
            % tr : double
            %   matrix trace of the tensor.
            
            tr = 0;
            for i = 1:numel(t)
                [mblocks, mcharges] = matrixblocks(t(i).var);
                qdims = qdim(mcharges);
                for j = 1:length(mblocks)
                    tr = tr + qdims(j) * trace(mblocks{j});
                end
            end
        end
        
        function t = transpose(t, p, r)
            % Compute the transpose of a tensor. This is defined as rotating the domain to
            % the codomain and vice versa, while cyclicly permuting the tensor blocks.
            % Currently not implemented.
            %
            % Usage
            % -----
            % :code:`t = transpose(t, p, rank)`
            % :code:`t = t.'`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % p : (1, :) int
            %   permutation vector, which must be cyclic. By default this is no permutation.
            %
            % r : (1, 2) int
            %   rank of the output tensor, by default equal to the rank of the input tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   transposed output tensor.
            
            error('tensors:TBA', 'This method has not been implemented.');
        end
        
        function t = uplus(t)
            
        end
        
        function t = uminus(t)
            t.var = -t.var;
        end
    end
    
    
    %% Factorizations
    methods
        function [V, D, W] = eig(A)
            % Compute the eigenvalues and eigenvectors of a square tensor.
            %
            % Usage
            % -----
            % :code:`D = eig(A)`
            %
            % :code:`[V, D] = eig(A)`
            %
            % :code:`[V, D, W] = eig(A)`
            %
            % Arguments
            % ---------
            % A : :class:`Tensor`
            %   square input tensor.
            %
            % Returns
            % -------
            % D : (:,:) :class:`Tensor`
            %   diagonal matrix of eigenvalues.
            %
            % V : (1,:) :class:`Tensor`
            %   row vector of right eigenvectors such that A * V = V * D.
            %
            % W : (1,:) :class:`Tensor`
            %   row vector of left eigenvectors such that W' * A = D * W'.
            
            assert(isequal(A.codomain, A.domain), 'tensors:ArgumentError', ...
                'Input should be square.');
            
            dims = struct;
            [mblocks, dims.charges] = matrixblocks(A);
            Ds = cell(size(mblocks));
            dims.degeneracies = zeros(size(mblocks));
            
            if nargout > 2
                Vs = cell(size(mblocks));
                Ws = cell(size(mblocks));
                for i = 1:length(mblocks)
                    [Vs{i}, Ds{i}, Ws{i}] = eig(mblocks{i});
                    dims.degeneracies(i) = size(Ds{i}, 1);
                end
            elseif nargout > 1
                Vs = cell(size(mblocks));
                for i = 1:length(mblocks)
                    [Vs{i}, Ds{i}] = eig(mblocks{i});
                    dims.degeneracies(i) = size(Ds{i}, 1);
                end
                
            else
                for i = 1:length(mblocks)
                    Ds{i} = diag(eig(mblocks{i}));
                    dims.degeneracies(i) = size(Ds{i}, 1);
                end
            end
            
            space = A.domain.new(dims, false);
            D = A.zeros(space, space);
            D.var = fill_matrix_data(D.var, Ds, dims.charges);
            
            if nargout <= 1
                V = D;
                return
            end
            if nargout > 1
                V = A.eye(A.domain, space);
                V.var = fill_matrix_data(V.var, Vs, dims.charges);
            end
            if nargout > 2
                W = A.eye(A.domain, space);
                W.var = fill_matrix_data(W.var, Ws, dims.charges);
            end
        end
        
        function [Q, R] = leftorth(t, p1, p2, alg)
            % Factorize a tensor into an orthonormal basis `Q` and remainder `R`, such that
            % :code:`permute(t, [p1 p2], [length(p1) length(p2)]) = Q * R`.
            %
            % Usage
            % -----
            % :code:`[Q, R] = leftorth(t, p1, p2, alg)`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor to factorize.
            %
            % p1, p2 : int
            %   partition of left and right indices, by default this is the partition of the
            %   input tensor.
            %
            % alg : char or string
            %   selection of algorithms for the decomposition:
            %   
            %   - 'qr' produces an upper triangular remainder R
            %   - 'qrpos' corrects the diagonal elements of R to be positive.
            %   - 'ql' produces a lower triangular remainder R
            %   - 'qlpos' corrects the diagonal elements of R to be positive.
            %   - 'polar' produces a Hermitian and positive semidefinite R.
            %   - 'svd' uses a singular value decomposition.
            %
            % Returns
            % -------
            % Q : :class:`Tensor`
            %   Orthonormal basis tensor
            %
            % R : :class:`Tensor`
            %   Remainder tensor, depends on selected algorithm.
            
            arguments
                t
                p1 = 1:t.rank(1)
                p2 = t.rank(1) + (1:t.rank(2))
                alg {mustBeMember(alg, {'qr', 'qrpos', 'ql', 'qlpos', 'polar', 'svd'})} ...
                    = 'qrpos'
            end
            
            if isempty(p1), p1 = 1:rank(t, 1); end
            if isempty(p2), p2 = rank(t, 1) + (1:rank(t,2)); end
            
            t = permute(t, [p1 p2], [length(p1) length(p2)]);
            
            dims = struct;
            [mblocks, dims.charges] = matrixblocks(t);
            
            Qs = cell(size(mblocks));
            Rs = cell(size(mblocks));
            dims.degeneracies = zeros(size(mblocks));
            
            for i = 1:length(mblocks)
                [Qs{i}, Rs{i}] = leftorth(mblocks{i}, alg);
                dims.degeneracies(i) = size(Qs{i}, 2);
            end
            
            V = t.codomain.new(dims, false);
            
            if strcmp(alg, 'polar')
                assert(isequal(V, prod(t.domain)));
                W = t.domain;
            elseif length(p1) == 1 && V == t.codomain
                W = t.codomain;
            elseif length(p2) == 1 && V == t.domain
                W = t.domain;
            else
                W = V;
            end
            
            Q = t.eye(t.codomain, W);
            Q.var = fill_matrix_data(Q.var, Qs, dims.charges);
            
            R = t.zeros(W, t.domain);
            R.var = fill_matrix_data(R.var, Rs, dims.charges);
        end
        
        function [R, Q] = rightorth(t, p1, p2, alg)
            % Factorize a tensor into an orthonormal basis `Q` and remainder `L`, such that
            % :code:`permute(t, [p1 p2], [length(p1) length(p2)]) = L * Q`.
            %
            % Usage
            % -----
            % :code:`[R, Q] = rightorth(t, p1, p2, alg)`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor to factorize.
            %
            % p1, p2 : int
            %   partition of left and right indices, by default this is the partition of the
            %   input tensor.
            %
            % alg : char or string
            %   selection of algorithms for the decomposition:
            %   
            %   - 'rq' produces an upper triangular remainder R
            %   - 'rqpos' corrects the diagonal elements of R to be positive.
            %   - 'lq' produces a lower triangular remainder R
            %   - 'lqpos' corrects the diagonal elements of R to be positive.
            %   - 'polar' produces a Hermitian and positive semidefinite R.
            %   - 'svd' uses a singular value decomposition.
            %
            % Returns
            % -------
            % R : :class:`Tensor`
            %   Remainder tensor, depends on selected algorithm.
            %
            % Q : :class:`Tensor`
            %   Orthonormal basis tensor.
            
            arguments
                t
                p1 = 1:t.rank(1)
                p2 = t.rank(1) + (1:t.rank(2))
                alg {mustBeMember(alg, {'rq', 'rqpos', 'lq', 'lqpos', 'polar', 'svd'})} ...
                    = 'rqpos'
            end

            if isempty(p1), p1 = 1:rank(t, 1); end
            if isempty(p2), p2 = rank(t, 1) + (1:rank(t,2)); end
            
            t = permute(t, [p1 p2], [length(p1) length(p2)]);
            
            dims = struct;
            [mblocks, dims.charges] = matrixblocks(t);
            Qs = cell(size(mblocks));
            Rs = cell(size(mblocks));
            dims.degeneracies = zeros(size(mblocks));
            
            for i = 1:length(mblocks)
                [Rs{i}, Qs{i}] = rightorth(mblocks{i}, alg);
                dims.degeneracies(i) = size(Qs{i}, 1);
            end
            
            V = t.domain.new(dims, false);
            
            if strcmp(alg, 'polar')
                assert(isequal(V, prod(t.codomain)));
                W = t.codomain;
            elseif length(p1) == 1 && V == t.codomain
                W = t.codomain;
            elseif length(p2) == 1 && V == t.domain
                W = t.domain;
            else
                W = V;
            end
            
            Q = t.eye(W, t.domain);
            Q.var = fill_matrix_data(Q.var, Qs, dims.charges);
            
            R = t.zeros(t.codomain, W);
            R.var = fill_matrix_data(R.var, Rs, dims.charges);
        end
        
        function N = leftnull(t, p1, p2, alg, atol)
            % Compute the left nullspace of a tensor, such that
            % :code:`N' * permute(t, [p1 p2], [length(p1) length(p2)]) = 0`.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor to compute the nullspace.
            %
            % p1, p2 : int
            %   partition of left and right indices, by default this is the partition of the
            %   input tensor.
            %
            % alg : char or string
            %   selection of algorithms for the nullspace:
            %
            %   - 'svd'
            %   - 'qr'
            %
            % Returns
            % -------
            % N : :class:`Tensor`
            %   orthogonal basis for the left nullspace.
            
            arguments
                t
                p1 = 1:t.rank(1)
                p2 = t.rank(1) + (1:t.rank(2))
                alg {mustBeMember(alg, {'svd', 'qr'})} = 'svd'
                atol = norm(t) * eps(underlyingType(t))
            end
            
            if isempty(p1), p1 = 1:rank(t, 1); end
            if isempty(p2), p2 = rank(t, 1) + (1:rank(t, 2)); end
            
            t = permute(t, [p1 p2], [length(p1) length(p2)]);
            
            dims = struct;
            [mblocks, dims.charges] = matrixblocks(t);
            Ns = cell(size(mblocks));
            dims.degeneracies = zeros(size(mblocks));
            
            for i = 1:length(mblocks)
                Ns{i} = leftnull(mblocks{i}, alg, atol);
                dims.degeneracies(i) = size(Ns{i}, 2);
            end
            
            N = t.eye(t.codomain, t.codomain.new(dims, false));
            N.var = fill_matrix_data(N.var, Ns, dims.charges);
        end
        
        function N = rightnull(t, p1, p2, alg, atol)
             % Compute the right nullspace of a tensor, such that
            % :code:`permute(t, [p1 p2], [length(p1) length(p2)]) * N = 0`.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor to compute the nullspace.
            %
            % p1, p2 : int
            %   partition of left and right indices, by default this is the partition of the
            %   input tensor.
            %
            % alg : char or string
            %   selection of algorithms for the nullspace:
            %
            %   - 'svd'
            %   - 'lq'
            %
            % Returns
            % -------
            % N : :class:`Tensor`
            %   orthogonal basis for the right nullspace.
            
            arguments
                t
                p1 = 1:t.rank(1)
                p2 = t.rank(1) + (1:t.rank(2))
                alg {mustBeMember(alg, {'svd', 'lq'})} = 'svd'
                atol = norm(t) * eps(underlyingType(t))
            end
            
            if isempty(p1), p1 = 1:rank(t, 1); end
            if isempty(p2), p2 = rank(t, 1) + (1:rank(t, 2)); end
            
            t = permute(t, [p1 p2], [length(p1) length(p2)]);
            
            dims = struct;
            [mblocks, dims.charges] = matrixblocks(t);
            Ns = cell(size(mblocks));
            dims.degeneracies = zeros(size(mblocks));
            
            for i = 1:length(mblocks)
                Ns{i} = rightnull(mblocks{i}, alg, atol);
                dims.degeneracies(i) = size(Ns{i}, 1);
            end
            
            N = Tensor.eye(t.domain.new(dims, false), t.domain);
            N.var = fill_matrix_data(N.var, Ns, dims.charges);
        end
        
        function [U, S, V, eta] = tsvd(t, p1, p2, trunc)
            % Compute the singular value decomposition of a tensor. This computes left and
            % right isometries U and V, and a non-negative diagonal tensor S such that
            % norm(permute(t, [p1 p2], [length(p1) length(p2)]) - U * S * V) = 0
            % Additionally, the dimension of S can be truncated in such a way to minimize
            % this norm, which gives the truncation error eta.
            %
            % Usage
            % -----
            % [U, S, V] = tsvd(t, p1, p2)
            % [U, S, V, eta] = tsvd(t, p1, p2, trunc, tol)
            % S = tsvd(t, ...)
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % p1, p2 : int
            %   partition of left and right indices, by default this is the partition of the
            %   input tensor.
            %
            % Keyword Arguments
            % -----------------
            % TruncDim : int
            %   truncate such that the size of S is not larger than this value.
            %
            % TruncBelow : numeric
            %   truncate such that there are no singular values below this value.
            %
            % TruncSpace : :class:`AbstractSpace`
            %   truncate such that the space of S is smaller than this value.
            %
            % Returns
            % -------
            % U, S, V : :class:`Tensor`
            %   left isometry U, non-negative diagonal S and right isometry V that satisfy
            %   U * S * V = permute(t, [p1 p2], [length(p1) length(p2)]).
            %
            % eta : numeric
            %   truncation error.
            
            arguments
                t
                p1 = 1:t.rank(1)
                p2 = t.rank(1) + (1:t.rank(2))
                trunc.TruncDim
                trunc.TruncBelow
                trunc.TruncSpace
            end
            
            t = permute(t, [p1 p2], [length(p1) length(p2)]);
            
            dims = struct;
            [mblocks, dims.charges] = matrixblocks(t);
            Us = cell(size(mblocks));
            Ss = cell(size(mblocks));
            Vs = cell(size(mblocks));
            dims.degeneracies = zeros(size(mblocks));
            
            doTrunc = ~isempty(fieldnames(trunc));
            if doTrunc, eta = 0; end
            for i = 1:length(mblocks)
                if doTrunc
                    [Us{i}, Ss{i}, Vs{i}] = svd(mblocks{i}, 'econ');
                    dims.degeneracies(i) = size(Ss{i}, 1);
                else
                    [Us{i}, Ss{i}, Vs{i}] = svd(mblocks{i});
                end
                Vs{i} = Vs{i}';
            end
            
            if isfield(trunc, 'TruncBelow')
                for i = 1:length(mblocks)
                    s = diag(Ss{i});
                    dims.degeneracies(i) = sum(s > trunc.TruncBelow);
                    eta = eta + sum(s(dims.degeneracies(i) + 1:end));
                    Us{i} = Us{i}(:, 1:dims.degeneracies(i));
                    Ss{i} = diag(s(1:dims.degeneracies(i)));
                    Vs{i} = Vs{i}(1:dims.degeneracies(i), :);
                end
            end
            if isfield(trunc, 'TruncDim')
                for i = 1:length(mblocks)
                    dims.degeneracies(i) = min(dims.degeneracies(i), trunc.TruncDim);
                    s = diag(Ss{i});
                    eta = eta + sum(s(dims.degeneracies(i) + 1:end));
                    Us{i} = Us{i}(:, 1:dims.degeneracies(i));
                    Ss{i} = diag(s(1:dims.degeneracies(i)));
                    Vs{i} = Vs{i}(1:1:dims.degeneracies(i), :);
                end
            end
            if isfield(trunc, 'TruncSpace')
                error('TBA');
            end
            
            if ~doTrunc
                W1 = prod(t.codomain);
                W2 = prod(t.domain);
                
                U = Tensor.eye(t.codomain, W1);
                S = Tensor.zeros(W1, W2);
                V = Tensor.eye(W2, t.domain);
            else
                W = t.domain.new(dims, false);
                
                U = Tensor.eye(t.codomain, W);
                S = Tensor.zeros(W, W);
                V = Tensor.eye(W, t.domain);
            end
            
            U.var = fill_matrix_data(U.var, Us, dims.charges);
            S.var = fill_matrix_data(S.var, Ss, dims.charges);
            V.var = fill_matrix_data(V.var, Vs, dims.charges);
        end
    end
    
    
    %% Matrix functions
    methods
        function t = expm(t)
            % Compute the matrix exponential of a square tensor. This is done via a scaling
            % and squaring algorithm with a Pade approximation, block-wise.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor.
            
            assert(isequal(t.codomain, t.domain), 'tensors:ArgumentError', ...
                'Input should be square.');
            
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                mblocks{i} = expm(mblocks{i});
            end
            t.var = fill_matrix_data(t.var, mblocks);
        end
        
        function t = inv(t)
            % Compute the matrix inverse of a square tensor, such that t * inv(t) = I.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   output tensor.
            
            assert(isequal(t.codomain, t.domain), 'tensors:ArgumentError', ...
                'Input should be square.');
            
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                mblocks{i} = inv(mblocks{i});
            end
            t.var = fill_matrix_data(t.var, mblocks);
        end
        
        function A = mpower(X, Y)
            % Raise a square tensor to a scalar power, or a scalar to a square tensor power.
            %
            % Usage
            % -----
            % :class:`A = x^Y`
            %
            % :class:`A = X^y`
            %
            % Arguments
            % ---------
            % X, Y : :class:`Tensor`
            %   Square input tensor.
            %
            % x, y : numeric
            %   Input scalars.
            %
            % Returns
            % -------
            % A : :class:`Tensor`
            %   Output tensor.
            
            % tensor to a scalar power
            if isnumeric(Y) && isscalar(Y)  
                assert(isequal(X.codomain, X.domain), 'tensors:ArgumentError', ...
                    'Input tensor should be square.');
                mblocks = matrixblocks(X);
                for i = 1:length(mblocks)
                    mblocks{i} = mblocks{i}^Y;
                end
                A = X;
                A.var = fill_matrix_data(A.var, mblocks);
                return
            end
            
            % scalar to a tensor power
            if isnumeric(X) && isscalar(X)
                assert(isequal(Y.codomain, Y.domain), 'tensors:ArgumentError', ...
                    'Input tensor should be square.');
                mblocks = matrixblocks(Y);
                for i = 1:length(mblocks)
                    mblocks{i} = X^mblocks{i};
                end
                A = Y;
                A.var = fill_matrix_data(A.var, mblocks);
                return
            end
            
            error('tensors:ArgumentError', ...
                'At least one of the arguments should be scalar.');
        end
        
        function [A, resnorm] = sqrtm(A)
            % Compute the principal square root of a square tensor. This is the unique root
            % for which every eigenvalue has nonnegative real part. If `t` is singular, then
            % the result may not exist.
            %
            % Arguments
            % ---------
            % A : :class:`Tensor`
            %   input tensor.
            %
            % Returns
            % -------
            % X : :class:`Tensor`
            %   principal square root of the input, which has X^2 = A.
            %
            % resnorm : numeric
            %   the relative residual, norm(A - X^2, 1) / norm(A, 1). If this argument is
            %   returned, no warning is printed if exact singularity is detected.
            
            assert(isequal(A.codomain, A.domain), 'tensors:ArgumentError', ...
                'Input should be square.');
            
            mblocks = matrixblocks(A);
            if nargout > 1
                resnorms = zeros(size(mblocks));
                for i = 1:length(mblocks)
                    [mblocks{i}, resnorms(i)] = sqrtm(mblocks{i});
                end
                resnorm = sum(resnorms);
            else
                for i = 1:length(mblocks)
                    mblocks{i} = sqrtm(mblocks{i});
                end
            end
            A.var = fill_matrix_data(A.var, mblocks);
        end
    end
    
    
    %% 
    methods
        function bool = isposdef(t)
            % Test if a tensor is a positive-definite map. Generally, a Hermitian matrix `M`
            % is positive-definite if the real number `z' * M * z`  is positive for every
            % nonzero complex column vector `z`.
            % This is equivalent to any of the following conditions:
            %
            % - M is Hermitian and all eigenvalues are real and positive.
            % - M is congruent with a diagonal matrix with positive real entries.
            % - There exists an invertible B such that `M = B' * B`.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % Returns
            % -------
            % bool : logical
            %   true if `t` is positive definite.
            
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                if ~isposdef(mblocks{i})
                    bool = false;
                    return
                end
            end
            bool = true;
        end
        
        function bool = isisometry(t, side, tol)
            % Test if a tensor is an isometric map. Generally, a matrix `M` is left or right
            % isometric if `M' * M = I` or `M * M' = I`.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % side : char
            %   either 'left', 'right' or 'both' (default).
            %
            % Keyword Arguments
            % -----------------
            % AbsTol, RelTol : numeric
            %   `norm(M * M' - eye(size(M))) < max(AbsTol, RelTol * norm(M))`.
            %   By default `AbsTol = 0` and `RelTol = eps`.
            %
            % Returns
            % -------
            % bool : logical
            %   true if t is isometric.
            
            arguments
                t
                side {mustBeMember(side, {'left', 'right', 'both'})} = 'both'
                tol.RelTol = sqrt(eps(underlyingType(t)))
                tol.AbsTol = 0
            end
            
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                if ~isisometry(mblocks{i}, side, 'RelTol', tol.RelTol, 'AbsTol', tol.AbsTol)
                    bool = false;
                    return
                end
            end
            bool = true;
        end
        
        function bool = istriu(t)
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                if ~istriu(mblocks{i})
                    bool = false;
                    return
                end
            end
            bool = true;
        end
        
        function bool = istril(t)
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                if ~istril(mblocks{i})
                    bool = false;
                    return
                end
            end
            bool = true;
        end
        
        function bool = isdiag(t)
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                if ~isdiag(mblocks{i})
                    bool = false;
                    return
                end
            end
            bool = true;
        end
        
        function bool = isreal(t)
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                if ~isreal(mblocks{i})
                    bool = false;
                    return
                end
            end
        end
        
        function r = cond(t, p)
            % Condition number with respect to inversion. This is defined as
            % :math:`||t|| * ||t^{-1}||` in the p-norm. For well conditioned
            % tensors, the value is near 1.0, while for badly conditioned tensors the value
            % diverges.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % p : 1, 2, inf or 'fro'
            %   kind of norm to use. Default is the 2-norm.
            %
            % Returns
            % -------
            % r : :class:`double`
            %   Condition number.
            
            arguments
                t
                p = 2
            end
            
            r = norm(t, p) * norm(inv(t), p);
        end
    end
    
    
    %% Solvers
    methods
        function varargout = linsolve(A, b, x0, M1, M2, options)
            % Find a solution for a linear system `A(x) = b` or `A * x = b`.
            %
            % Arguments
            % ---------
            % A : operator
            %   either a function handle implementing or an object that supports
            %   right multiplication.
            %
            % b : :class:`Tensor`
            %   right-hand side of the equation, interpreted as vector.
            %
            % x0 : :class:`Tensor`
            %   optional initial guess for the solution.
            %
            % M1, M2 : operator
            %   preconditioner M = M1 or M = M1 * M2 to effectively solve the system A *
            %   inv(M) * y = b with y = M * x.
            %   M is either a function handle implementing or an object that supports
            %   left division.
            %
            % Keyword Arguments
            % -----------------
            % Tol : numeric
            %   specifies the tolerance of the method, by default this is the square root of
            %   eps.
            %
            % Algorithm : char
            %   specifies the algorithm used. Can be either one of the following:
            %
            %   - 'bicgstab'
            %   - 'bicgstabl'
            %   - 'gmres'
            %   - 'pcg'
            %
            % MaxIter : int
            %   Maximum number of iterations.
            %
            % Restart : int
            %   For 'gmres', amount of iterations after which to restart.
            %
            % Verbosity : int
            %   Level of output information, by default nothing is printed if `flag` is
            %   returned, otherwise only warnings are given.
            %
            %   - 0 : no information
            %   - 1 : information at failure
            %   - 2 : information at convergence
            %
            % Returns
            % -------
            % x : :class:`Tensor`
            %   solution vector.
            %
            % flag : int
            %   a convergence flag:
            %
            %   - 0 : linsolve converged to the desired tolerance.
            %   - 1 : linsolve reached the maximum iterations without convergence.
            %   - 2 : linsolve preconditioner was ill-conditioned.
            %   - 3 : linsolve stagnated.
            %   - 4 : one of the scalar quantities calculated became too large or too small.
            %
            % relres : numeric
            %   relative residual, norm(b - A * x) / norm(b).
            %
            % iter : int
            %   iteration number at which x was computed.
            %
            % resvec : numeric
            %   vector of estimated residual norms at each part of the iteration.
            
            arguments
                A
                b
                x0 = []
                M1 = []
                M2 = []
                
                options.Tol = eps(underlyingType(b)) ^ (3/4)
                options.Algorithm {mustBeMember(options.Algorithm, ...
                    {'pcg', 'gmres', 'bicgstab', 'bicgstabl'})} = 'gmres'
                options.MaxIter = 400
                options.Restart = 30
                options.Verbosity = 0
            end
            
            % Convert input objects to vectors
            b_vec = vectorize(b);
            b_sz = size(b_vec);
            
            if ~isempty(x0)
                x0_vec = vectorize(x0);
            else
                x0_vec = [];
            end
            
            % Convert input operators to handle vectors
            if isa(A, 'function_handle')
                A_fun = @(x) vectorize(A(devectorize(x, b)));
            else
                A_fun = @(x) vectorize(A * devectorize(x, b));
            end
            
            if isempty(M1)
                M1_fun = [];
            elseif isa(M1, 'function_handle')
                M1_fun = @(x) vectorize(M1(devectorize(x, b)));
            else
                M1_fun = @(x) vectorize(M1 \ devectorize(x, b));
            end
            
            if isempty(M2)
                M2_fun = [];
            elseif isa(M2, 'function_handle')
                M2_fun = @(x) vectorize(M2(devectorize(x, b)));
            else
                M2_fun = @(x) vectorize(M2 \ devectorize(x, b));
            end
            
            % Sanity check on parameters
            options.Restart = min(options.Restart, b_sz(1));
            if options.Tol < eps(underlyingType(b))^0.9
                warning('Requested tolerance might be too strict.');
            end
            
            % Apply MATLAB implementation
            switch options.Algorithm
                case 'bicgstab'
                    [varargout{1:nargout}] = bicgstab(A_fun, b_vec, ...
                        options.Tol, options.MaxIter, M1_fun, M2_fun, x0_vec);
                case 'bicgstabl'
                    [varargout{1:nargout}] = bicgstabl(A_fun, b_vec, ...
                        options.Tol, options.MaxIter, M1_fun, M2_fun, x0_vec);
                case 'gmres'
                    options.MaxIter = min(b_sz(1), options.MaxIter);
                    [varargout{1:nargout}] = gmres(A_fun, b_vec, ...
                        options.Restart, options.Tol, options.MaxIter, ...
                        M1_fun, M2_fun, x0_vec);
                case 'pcg'
                    [varargout{1:nargout}] = pcg(A_fun, b_vec, ...
                        options.Tol, options.MaxIter, M1_fun, M2_fun, x0_vec);
            end
            
            
            % Convert output
            varargout{1} = devectorize(varargout{1}, b);
        end
        
        function varargout = eigsolve(A, x0, howmany, sigma, options)
            % Find a few eigenvalues and eigenvectors of an operator.
            %
            % Usage
            % -----
            % :code:`[V, D, flag] = eigsolve(A, x0, howmany, sigma, kwargs)`
            % :code:`D = eigsolve(A, x0, ...)`
            %
            % Arguments
            % ---------
            % A : :class:`Tensor` or function_handle
            %   A square tensormap interpreted as matrix.
            %   A function handle which implements one of the following, depending on sigma:
            %
            %   - A \ x, if `sigma` is 0 or 'smallestabs'
            %   - (A - sigma * I) \ x, if sigma is a nonzero scalar
            %   - A * x, for all other cases
            %
            % x0 : :class:`Tensor`
            %   initial guess for the eigenvector. If A is a :class:`Tensor`, this defaults
            %   to a random complex :class:`Tensor`, for function handles this is a required
            %   argument.
            %
            % howmany : int
            %   amount of eigenvalues and eigenvectors that should be computed. By default
            %   this is 1, and this should not be larger than the total dimension of A.
            %
            % sigma : 'char' or numeric
            %   selector for the eigenvalues, should be either one of the following:
            %
            %   - 'largestabs', 'largestreal', 'largestimag' : default, eigenvalues of
            %       largest magnitude, real part or imaginary part.
            %   - 'smallestabs', 'smallestreal', 'smallestimag' : eigenvalues of smallest
            %       magnitude, real part or imaginary part.
            %   - numeric : eigenvalues closest to sigma.
            %
            % Keyword Arguments
            % -----------------
            % Tol : numeric
            %   tolerance of the algorithm.
            %
            % Algorithm : char
            %   choice of algorithm. Currently only 'eigs' is available, which leverages the
            %   default Matlab eigs.
            %
            % MaxIter : int
            %   maximum number of iterations, 100 by default.
            %
            % KrylovDim : int
            %   number of vectors kept in the Krylov subspace.
            %
            % IsSymmetric : logical
            %   flag to speed up the algorithm if the operator is symmetric, false by
            %   default.
            %
            % Verbosity : int
            %   Level of output information, by default nothing is printed if `flag` is
            %   returned, otherwise only warnings are given.
            %
            %   - 0 : no information
            %   - 1 : information at failure
            %   - 2 : information at convergence
            %   - 3 : information at every iteration
            %
            % Returns
            % -------
            % V : (1, howmany) :class:`Tensor`
            %   vector of eigenvectors.
            %
            % D : numeric
            %   vector of eigenvalues if only a single output argument is asked, diagonal
            %   matrix of eigenvalues otherwise.
            %
            % flag : int
            %   if flag = 0 then all eigenvalues are converged, otherwise not.
            
            arguments
                A
                x0 = A.randnc(A.domain, [])
                howmany = 1
                sigma = 'largestabs'
                
                options.Tol = 1e-12
                options.Algorithm = 'eigs'
                options.MaxIter = 100
                options.KrylovDim = 20
                options.IsSymmetric logical = false
                options.Verbosity = 0
            end
            
            assert(isnumeric(sigma) || ismember(sigma, {'largestabs', 'smallestabs', ...
                'largestreal', 'smallestreal', 'bothendsreal', ...
                'largestimag', 'smallestimag', 'bothendsimag'}), ...
                'tensors:ArgumentError', 'Invalid choice of eigenvalue selector.');
            
            x0_vec = vectorize(x0);
            sz = size(x0_vec);
            
            if isa(A, 'function_handle')
                A_fun = @(x) vectorize(A(devectorize(x, x0)));
            else
                A_fun = @(x) vectorize(A * devectorize(x, x0));
            end
            
            options.KrylovDim = min(sz(1), options.KrylovDim);
            
            [varargout{1:nargout}] = eigs(A_fun, sz(1), howmany, sigma, ...
                'Tolerance', options.Tol, 'MaxIterations', options.MaxIter, ...
                'SubspaceDimension', options.KrylovDim, 'IsFunctionSymmetric', ...
                options.IsSymmetric, 'StartVector', x0_vec);
            if nargout > 1
                for i = howmany:-1:1
                    V(:, i) = devectorize(varargout{1}(:, i), x0);
                end
                varargout{1} = V;
            end
        end
        
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
    
    
    %% Converters
    methods
        function a = double(t)
            % Convert tensor to array of double.
            %
            % Arguments
            % ---------
            % t : Tensor
            %
            % Returns
            % -------
            % a : double
            
            trees = fusiontrees(t);
            blocks = tensorblocks(t);
            if isempty(trees)
                a = blocks{1};
                return
            end
            
            % Initialize output
            a = zeros(dims(t));
            s = [t.codomain flip(t.domain)];
            dimsizes = cell(size(s));
            for i = 1:length(s)
                dimsizes{i} = qdim(charges(s(i))) .* degeneracies(s(i));
            end
            a_cell = mat2cell(a, dimsizes{:});
            
            % Locate non-empty blocks
            [~, locb] = ismember(trees.uncoupled, ...
                charges(s).', 'rows');
            
            % Fill output
            if fusionstyle(trees) == FusionStyle.Unique
                a_cell(locb) = blocks;
            else
                tree_array = fusiontensor(trees);
                for i = 1:length(locb)
                    tree_double = tree_array{i};
                    tree_size = size(tree_double, 1:nspaces(t));
                    block_size = size(blocks{i}, 1:nspaces(t));
                    a_cell{locb(i)} = a_cell{locb(i)} + ...
                        reshape(contract(tree_double, -(1:2:2*nspaces(t)), ...
                        blocks{i}, -(2:2:2*nspaces(t))), tree_size .* block_size);
                end
            end
            
            a = cell2mat(a_cell);
        end
    end
    
    
    %% Utility
    methods
        function s = GetMD5_helper(t)
            s = {t.codomain t.domain};
        end
        
        function type = underlyingType(t)
            type = underlyingType(t(1).var);
        end
        
        function disp(t)
            if isscalar(t)
                r = t.rank;
                fprintf('Rank (%d, %d) %s:\n\n', r(1), r(2), class(t));
                s = space(t);
                for i = 1:length(s)
                    fprintf('%d.\t', i);
                    disp(s(i));
                    fprintf('\n');
                end
            else
                builtin('disp', t);
            end
        end
    end
end
