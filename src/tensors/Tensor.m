classdef Tensor < AbstractTensor
    % Tensor - Base implementation of a dense tensor array with optional symmetries.
    
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
                for i = numel(varargin{1}):-1:1
                    t(i).codomain = varargin{1}(i).codomain;
                    t(i).domain = varargin{1}(i).domain;
                    t(i).var = varargin{1}(i).var;
                end
                t = reshape(t, size(varargin{1}));
                return
            end
            
            if nargin == 1 && isa(varargin{1}, 'SparseTensor')
                t = full(varargin{1});
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
                
                if isa(codomain, 'SumSpace') || isa(domain, 'SumSpace')
                    if isempty(codomain)
                        sz = flip(nsubspaces(domain));
                    elseif isempty(domain)
                        sz = nsubspaces(codomain);
                    else
                        if ~isa(codomain, 'SumSpace'), codomain = SumSpace(codomain); end
                        if ~isa(domain, 'SumSpace'), domain = SumSpace(domain); end
                        sz = [nsubspaces(codomain) flip(nsubspaces(domain))];
                    end
                    subs = ind2sub_(sz, 1:prod(sz));
                    for i = size(subs, 1):-1:1
                        [cod, dom] = slice(codomain, domain, subs(i,:));
                        t(i).codomain = cod;
                        t(i).domain = dom;
                        t(i).var = AbstractBlock.new(cod, dom);
                    end
                else
                    t.domain    = domain;
                    t.codomain  = codomain;
                    t.var = AbstractBlock.new(codomain, domain);
                end
                
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
                if isempty(charges)
                    t.var = fill_matrix_data(t.var, data);
                else
                    t.var = fill_matrix_data(t.var, data, charges);
                end
            else
                if isempty(charges)
                    t.var = fill_matrix_fun(t.var, data);
                else
                    t.var = fill_matrix_fun(t.var, data, charges);
                end
            end
        end
        
        function t = fill_tensor(t, data)
            % Fill the tensor blocks of a tensor.
            %
            % Usage
            % -----
            % :code:`t = fill_tensor(t, tensors)`
            %
            % :code:`t = fill_tensor(t, fun)`
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
            % Returns
            % -------
            % t : :class:`Tensor`
            %   filled tensor.
            
            arguments
                t
                data
            end
            
            if isnumeric(data), data = {data}; end
            
            if iscell(data)
                t.var = fill_tensor_data(t.var, data);
            else
                [tmp, trees] = tensorblocks(t);
                for i = 1:length(tmp)
                    tmp{i} = data(size(tmp{i}), trees(i));
                end
                t.var = fill_tensor_data(t.var, tmp);
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
                kwargs.Mode {mustBeMember(kwargs.Mode, {'matrix', 'tensor'})} ...
                    = 'matrix'
            end
            
            % Parse special constructors
            if ~isempty(varargin{1}) && isnumeric(varargin{1})
                if length(varargin) == 1
                    spaces = CartesianSpace.new(varargin{1});
                elseif isnumeric(varargin{2}) || islogical(varargin{2})
                    assert(length(varargin{1}) == length(varargin{2}))
                    args = [num2cell(varargin{1}); num2cell(varargin{2})];
                    spaces = ComplexSpace.new(args{:});
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
                switch kwargs.Mode
                    case 'matrix'
                        t = fill_matrix(t, fun);
                    case 'tensor'
                        t = fill_tensor(t, fun);
                end
            end
        end
        
        function t = zeros(varargin)
            t = Tensor.new(@zeros, varargin{:});
        end
        
        function t = ones(varargin)
            t = Tensor.new(@ones, varargin{:});
        end
        
        function t = eye(varargin)
            t = Tensor.new(@eye, varargin{:});
        end
        
        function t = rand(varargin)
            t = Tensor.new(@rand, varargin{:});
        end
        
        function t = randn(varargin)
            t = Tensor.new(@randn, varargin{:});
        end
        
        function t = randc(varargin)
            t = Tensor.new(@randc, varargin{:});
        end
        
        function t = randnc(varargin)
            t = Tensor.new(@randnc, varargin{:});
        end
    end
    
    methods (Hidden)
        function tdst = zerosLike(t, varargin)
            tdst = repmat(0 * t, varargin{:});
        end
    end
    
    
    %% Structure
    methods
        function n = indin(t)
            n = length(t(1).domain);
        end
        
        function n = indout(t)
            n = length(t(1).codomain);
        end
        
        function n = nspaces(t)
            n = length(t(1).domain) + length(t(1).codomain);
        end
        
        function r = rank(t, i)
            r = [length(t(1).codomain) length(t(1).domain)];
            if nargin > 1
                r = r(i);
            end
        end
        
        function sp = space(t, inds)
            if isscalar(t)
                sp = [t.codomain t.domain'];
            else
                [cod, dom] = deduce_spaces(t);
                sp = [cod, dom'];
            end
            if nargin > 1
                sp = sp(inds);
            end
        end
        
        function f = fusiontrees(t)
            f = fusiontrees(t.codomain, t.domain);
        end
        
        function [b, f] = tensorblocks(t)
            b = tensorblocks(t.var);
            if nargout > 1, f = fusiontrees(t); end
        end
        
        function varargout = matrixblocks(t)
            [varargout{1:nargout}] = matrixblocks(t.var);
        end
        
        function style = braidingstyle(t)
            style = braidingstyle(t.codomain, t.domain);
        end
        
        function style = fusionstyle(t)
            style = fusionstyle(t.codomain, t.domain);
        end
        
        function t = full(t)
        end
        
        function tdst = insert_onespace(tsrc, i, dual)
            % insert a trivial space at position i.
            arguments
                tsrc
                i = nspaces(tsrc) + 1
                dual = false
            end
            
            spaces = insertone(space(tsrc), i, dual);
            data = matrixblocks(tsrc);
            
            r = rank(tsrc);
            if i <= r(1)
                r(1) = r(1) + 1;
            else
                r(2) = r(2) + 1;
            end
            tdst = Tensor(spaces(1:r(1)), spaces(r(1)+1:end)');
            tdst.var = fill_matrix_data(tdst.var, data);
        end
        
        function tdst = embed(tsrc, tdst)
            % embed a tensor in a different tensor.
            
            assert(isequal(rank(tsrc), rank(tdst)), 'tensors:argerror', ...
                'tensors must have the same rank');
            
            [bsrc, fsrc] = tensorblocks(tsrc);
            [bdst, fdst] = tensorblocks(tdst);
            
            [lia, locb] = ismember(fsrc, fdst);
            nsp = nspaces(tdst);
            
            for i = find(lia).'
                sz = min(size(bsrc{i}, 1:nsp), size(bdst{locb(i)}, 1:nsp));
                inds = arrayfun(@(x) 1:x, sz, 'UniformOutput', false);
                bdst{locb(i)}(inds{:}) = bsrc{i}(inds{:});
            end
            
            tdst.var = fill_tensor_data(tdst.var, bdst);
        end
    end
    
    
    %% Comparison
    methods
        
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
            
            for i = 1:numel(t)
                t(i) = tpermute(t(i)', nspaces(t(i)):-1:1, rank(t(i)));
            end
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
                [t(i).codomain, t(i).domain] = swapvars(t(i).codomain, t(i).domain);
                t(i).var = t(i).var';
            end
            
            t = permute(t, ndims(t):-1:1);
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
            
            assert(isequal(size(t1), size(t2)), 'tensors:dimerror', ...
                'input tensors must have the same size.');
            
            d = 0;
            for i = 1:numel(t1)
                assert(isequal(t1(i).domain, t2(i).domain) && ...
                    isequal(t1(i).codomain, t2(i).codomain), ...
                    'tensors:SpaceMismatch', ...
                    'dot product only defined for tensors of equal structure.');
                [mblocks1, mcharges] = matrixblocks(t1(i).var);
                mblocks2 = matrixblocks(t2(i).var);
                qdims = qdim(mcharges);
                for j = 1:length(mblocks1)
                    d = d + qdims(j) * sum(conj(mblocks1{j}) .* mblocks2{j}, 'all');
                end
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
                elseif iszero(t2(i))
                    continue;
                elseif iszero(t1(i))
                    t1(i) = -t2(i);
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
            
            if isnumeric(A)
                if issparse(A) && ~all(any(A, 2))
                    [cod, dom] = deduce_spaces(B);
                    cod2 = cell(size(cod));
                    for i = flip(1:length(cod))
                        tmp = subspaces(cod(i));
                        cod2{i} = tmp(1);
                    end
                    C = SparseTensor.zeros(SumSpace(cod2{:}), dom);
                else
                    C = Tensor.empty(szA(1), szB(2), 0);
                end
                
                for i = 1:szA(1)
                    for j = 1:szB(2)
                        C(i, j) = sum(A(i, :).' .* B(:, j), 'all');
                    end
                end
                return
            end
            
            if isnumeric(B)
                if issparse(B) && ~all(any(B, 1))
                    [cod, dom] = deduce_spaces(A);
                    dom2 = cell(size(dom));
                    for i = 1:length(dom)
                        tmp = subspaces(dom(i));
                        dom2{i} = tmp(1);
                    end
                    C = SparseTensor.zeros(cod, SumSpace(dom2{:}));
                end
                C(szA(1), szB(2)) = Tensor();
                for i = flip(1:szA(1))
                    for j = flip(1:szB(2))
                        C(i, j) = sum(A(i, :) .* B(:, j).', 'all');
                    end
                end
                return
            end
            
            % Tensor multiplications
            for i = szA(1):-1:1
                for j = szB(2):-1:1
                    C(i, j) = sum(reshape(A(i, :), [], 1) .* B(:, j));
                end
            end
        end
        
        function n = nnz(A)
            n = numel(A);
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
                        if iszero(t(i)), continue; end
                        mblocks = matrixblocks(t(i).var);
                        for j = 1:length(mblocks)
                            n = max(norm(mblocks{j}, p), n);
                        end
                    end
                    
                case 'fro'
                    n = 0;
                    for i = 1:numel(t)
                        if iszero(t(i)), continue; end
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
        
        function [t, n] = normalize(t)
            n = zeros(size(t));
            for i = 1:numel(t)
                n(i) = norm(t(i));
                t(i) = t(i) .* (1 / n(i));
            end
        end
        
        function t1 = plus(t1, t2)
            if isnumeric(t1), t1 = t2 + t1; return; end
            if ~isequal(size(t1), size(t2))
                error('plus:dimagree', ...
                    'Incompatible sizes for vectorized plus. (%s) (%s)', ...
                    dim2str(size(t1)), dim2str(size(t2)));
            end
            
            for i = 1:numel(t1)
                if isnumeric(t2(i))
                    assert(isequal(t1(i).codomain, t1(i).domain), 'tensors:NonSquare', ...
                        'Addition with scalars only defined for square tensors.');
                    I = t1(i).eye(t1(i).codomain, t1(i).domain);
                    t1(i).var = t1(i).var + I.var .* t2(i);
                elseif iszero(t2(i))
                    continue;
                elseif iszero(t1(i))
                    t1(i) = t2(i);
                else
                    if ~isequal(t1(i).domain, t2(i).domain) || ...
                            ~isequal(t1(i).codomain, t2(i).codomain)
                        t1_str = sprintf('%s\t<-\t%s', join(string(t1(i).codomain)), ...
                            join(string(t1(i).domain)));
                        t2_str = sprintf('%s\t<-\t%s', join(string(t2(i).codomain)), ...
                            join(string(t2(i).domain)));
                        error('tensors:spacemismatch', ...
                            'Added spaces incompatible.\n%s\n%s', t1_str, t2_str);
                    end
                    t1(i).var = t1(i).var + t2(i).var;
                end
            end
        end
        
        function t = tpermute(t, p, r)
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
                p = []
                r = []
            end
            
            if ~isscalar(t)
                for i = 1:numel(t)
                    t(i) = tpermute(t(i), p, r);
                end
                t = permute(t, p);
                return
            end
            
            if isempty(p), p = 1:nspaces(t); end
            if isempty(r), r = rank(t); end
            
            assert(length(p) == nspaces(t), 'Invalid permutation.');
            assert(length(p) == sum(r), 'Invalid new rank.');
            
            if (all(p == 1:nspaces(t)) && all(rank(t) == r)), return; end
            
            global cache
            if isempty(cache), cache = LRU; end
            
            if Options.CacheEnabled()
                key = GetMD5({t.codomain, t.domain, p, r}, ...
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
            
            for i = 1:numel(t)
                t(i) = tpermute(t(i), 1:nspaces(t(i)), r);
                assert(sum(r) == sum(rank(t(i))), ...
                    'tensors:ValueError', 'Invalid new rank.');
            end
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
            
            if isscalar(a)
                for i = 1:numel(t)
                    t(i).var = rdivide(t(i).var, a);
                end
                return
            end
            
            error('undefined');
        end
        
        function C = sum(A, dim)
            arguments
                A
                dim = []
            end
            
            if isscalar(A), C = A; return; end
            
            if isempty(dim), dim = find(size(A) ~= 1, 1); end
            
            if strcmp(dim, 'all')
                C = A(1);
                for i = 2:numel(A)
                    C = C + A(i);
                end
                return
            end
            
            if ismatrix(A)
                if dim == 1
                    C = A(1, :);
                    for i = 2:size(A, 1)
                        C = C + A(i, :);
                    end
                    return
                end
                
                if dim == 2
                    C = A(:, 1);
                    for i = 2:size(A, 2)
                        C = C + A(:, i);
                    end
                    return
                end
            end
            
            error('TBA');
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
                options.NumDimensionsA = ndims(A)
            end
            
            if ~isscalar(A) || ~isscalar(B)
                szA = size(A, 1:options.NumDimensionsA);
                szB = size(B, 1:max(ndims(B), max(dimB)));
                
                assert(all(szA(dimA) == szB(dimB)), 'tensors:SizeMismatch', ...
                    'Invalid contraction sizes.');
                
                uncA = 1:length(szA); uncA(dimA) = [];
                uncB = 1:length(szB); uncB(dimB) = [];
                
                if isempty(uncA)
                    if isempty(uncB)
                        szC = [1 1];
                    elseif length(uncB) == 1
                        szC = [1 szB(uncB)];
                    else
                        szC = szB(uncB);
                    end
                elseif isempty(uncB)
                    if length(uncA) == 1
                        szC = [szA(uncA) 1];
                    else
                        szC = szA(uncA);
                    end
                else
                    szC = [szA(uncA) szB(uncB)];
                end
                
                A = reshape(permute(A, [uncA dimA]), prod(szA(uncA)), prod(szA(dimA)));
                B = reshape(permute(B, [dimB, uncB]), prod(szB(dimB)), prod(szB(uncB)));
                
                for i = prod(szA(uncA)):-1:1
                    for j = prod(szB(uncB)):-1:1
                        C(i,j) = tensorprod(A(i,1), B(1,j), dimA, dimB, ca, cb, ...
                            'NumDimensionsA', options.NumDimensionsA);
                        for k = 2:prod(szA(dimA))
                            C(i,j) = C(i,j) + ...
                                tensorprod(A(i,k), B(k,j), dimA, dimB, ca, cb, ...
                                'NumDimensionsA', options.NumDimensionsA);
                        end
                    end
                end
                C = reshape(C, szC);
                return
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
            
            global cache
            if isempty(cache), cache = LRU; end
            
            if Options.CacheEnabled()
                key = GetMD5({A.codomain, A.domain, ...
                    B.codomain, B.domain, ...
                    dimA, dimB, ca, cb}, 'Array', 'hex');
                med = get(cache, key);
                if isempty(med)
                    med = struct;
                    
                    A_ = similar(@(x,charge) uninit(x), A, iA, 'Rank', rA);
                    med.varA = A_.var;
                    [med.mapA, f] = permute(fusiontrees(A), iA, rA);
                    for i = rA(1) + (1:rA(2))
                        if ~isdual(space(A_, i))
                            [c2, f] = twist(f, i);
                            med.mapA = med.mapA * c2;
                        end
                    end
                    
                    B_ = similar(@(x,charge) uninit(x), B, iB, 'Rank', rB);
                    med.varB = B_.var;
                    med.mapB = permute(fusiontrees(B), iB, rB);
                    
                    if ~isequal(A_.domain, B_.codomain)
                        error('tensors:SpaceMismatch', ...
                            'Contracted spaces incompatible.\n%s\n%s', ...
                            join(string(A_.domain), '    '), ...
                            join(string(B_.codomain), '    '));
                    end
                    if ~isempty(A_.codomain) || ~isempty(B_.domain)
                        med.C = Tensor.zeros(A_.codomain, B_.domain);
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
                    C = sum(cellfun(@mtimes, Ablocks, Bblocks), 'all'); % horzcat(Ablocks{:}) * vertcat(Bblocks{:});
                else
                    C.var = mul(C.var, varA, varB);
                end
                return
            end
            
            A = tpermute(A, iA, rA);
            for i = rA(1) + (1:rA(2))
                if ~isdual(space(A, i)), A = twist(A, i); end
            end
            %             A = twist(A, [false(1, length(A.codomain)) ~isdual(A.domain')]);
            B = tpermute(B, iB, rB);
            
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
        
        function C = times(A, B)
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
            
            if isscalar(A) && ~isscalar(B)
                A = repmat(A, size(B));
            elseif isscalar(B) && ~isscalar(A)
                B = repmat(B, size(A));
            end
            
            nd = max(ndims(A), ndims(B));
            assert(isequal(size(A, 1:nd), size(B, 1:nd)), ...
                'times:dimagree', 'incompatible dimensions.');
            
            if isnumeric(A)
                if issparse(A)
                    [cod, dom] = deduce_spaces(B);
                    C = SparseTensor.zeros(cod, dom);
                    I = find(A);
                    if isempty(I), return; end
                    C(I) = full(A(I)) .* B(I);
                    return
                end
                
                C = B;
                for i = 1:numel(C)
                    C(i).var = C(i).var .* A(i);
                end
                return
            end
            
            if isnumeric(B)
                C = B .* A;
                return
            end
            
            for i = numel(A):-1:1
                assert(isequal(A(i).domain, B(i).codomain), 'tensors:SpaceMismatch', ...
                    'Multiplied spaces incompatible.');
                if ~isempty(A(i).codomain) || ~isempty(B(i).domain)
                    C(i) = Tensor.zeros(A(i).codomain, B(i).domain);
                    C(i).var = mul(C(i).var, A(i).var, B(i).var);
                else
                    Ablocks = matrixblocks(A(i).var);
                    Bblocks = matrixblocks(B(i).var);
                    C(i) = horzcat(Ablocks{:}) * vertcat(Bblocks{:});
                end
            end
            C = reshape(C, size(A));
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
            if nargin < 2
                p = circshift(1:nspaces(t), length(t.domain));
            else
                assert(iscircperm(p));
            end
            
            if nargin < 3
                r = flip(rank(t));
            end
            
            t = tpermute(t, p, r);
        end
        
        function t = twist(t, i, inv)
            % Twist the spaces of a tensor.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % i : (1, :) int or logical
            %   indices to twist.
            %
            % inv : logical
            %   flag to indicate inverse twisting.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   twisted tensor with desired rank.
            
            arguments
                t
                i
                inv = false
            end
            
            if isempty(i) || ~any(i) || istwistless(braidingstyle(t(1)))
                return
            end
            
            if numel(t) > 1
                for i = 1:numel(t)
                    t(i) = twist(t(i), i, inv);
                end
                return
            end
            
            global cache
            if isempty(cache), cache = LRU; end
            
            if Options.CacheEnabled()
                key = GetMD5({t.codomain, t.domain, i, inv}, ...
                    'Array', 'hex');
                med = get(cache, key);
                if isempty(med)
                    med = twist(fusiontrees(t), i, inv);
                    cache = set(cache, key, med);
                end
            else
                med = twist(fusiontrees(t), i, inv);
            end
            
            t.var = axpby(1, t.var, 0, t.var, 1:nspaces(t), med);
        end
        
        function t = twistdual(t, i, inv)
            % Twist the spaces of a tensor if they are dual.
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % i : (1, :) int or logical
            %   indices to twist.
            %
            % inv : logical
            %   flag to indicate inverse twisting.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   twisted tensor with desired rank.
            arguments
                t
                i
                inv = false
            end
            
            for j = 1:numel(t)
                i_dual = i(isdual(space(t(j), i)));
                t(j) = twist(t(j), i_dual, inv);
            end
        end
        
        function t = uplus(t)
            
        end
        
        function t = uminus(t)
            for i = 1:numel(t)
                t(i).var = -t(i).var;
            end
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
            % :code:`tpermute(t, [p1 p2], [length(p1) length(p2)]) = Q * R`.
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
            if isempty(p2), p2 = rank(t, 1) + (1:rank(t, 2)); end
            
            t = tpermute(t, [p1 p2], [length(p1) length(p2)]);
            
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
                assert(isequal(V, prod(t.domain)), ...
                    'linalg:polar', 'polar decomposition should lead to square R.');
                W = V;
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
            % :code:`tpermute(t, [p1 p2], [length(p1) length(p2)]) = L * Q`.
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
            
            t = tpermute(t, [p1 p2], [length(p1) length(p2)]);
            
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
                assert(isequal(V, prod(t.codomain)), ...
                    'linalg:polar', 'polar decomposition should lead to square R.');
                W = V;
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
            % :code:`N' * tpermute(t, [p1 p2], [length(p1) length(p2)]) = 0`.
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
            
            [U, S] = tsvd(t, p1, p2);
            
            dims = struct;
            [Sblocks, c1] = matrixblocks(S);
            [Ublocks, c2] = matrixblocks(U);
            
            lia = ismember_sorted(c2, c1);
            ctr = 0;
            for i = 1:length(Ublocks)
                if ~lia(i), continue; end
                ctr = ctr + 1;
                u = Ublocks{i};
                s = Sblocks{ctr};
                if isvector(s)
                    diags = s(1);
                else
                    diags = diag(s);
                end
                r = sum(diags > atol);
                Ublocks{i} = u(:, (r + 1):size(u, 1));
            end
            
            dims.degeneracies = cellfun(@(x) size(x, 2), Ublocks);
            dims.charges = c2;
            
            mask = dims.degeneracies > 0;
            dims.charges = c2(mask);
            dims.degeneracies = dims.degeneracies(mask);
            Ns = Ublocks(mask);
            
            W = t.codomain.new(dims, false);
            N = t.eye(U.codomain, W);
            N.var = fill_matrix_data(N.var, Ns, dims.charges);
        end
        
        function N = rightnull(t, p1, p2, alg, atol)
            % Compute the right nullspace of a tensor, such that
            % :code:`tpermute(t, [p1 p2], [length(p1) length(p2)]) * N = 0`.
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
            
            [~, S, V] = tsvd(t, p1, p2);
            
            dims = struct;
            [Sblocks, c1] = matrixblocks(S);
            [Vblocks, c2] = matrixblocks(V);
            
            lia = ismember_sorted(c2, c1);
            ctr = 0;
            for i = 1:length(Vblocks)
                if ~lia(i), continue; end
                ctr = ctr + 1;
                v = Vblocks{i};
                s = Sblocks{ctr};
                if isvector(s)
                    diags = s(1);
                else
                    diags = diag(s);
                end
                r = sum(diags > atol);
                Vblocks{i} = v((r + 1):size(v, 2), :);
            end
            
            dims.degeneracies = cellfun(@(x) size(x, 1), Vblocks);
            dims.charges = c2;
            
            mask = dims.degeneracies > 0;
            dims.charges = c2(mask);
            dims.degeneracies = dims.degeneracies(mask);
            Ns = Vblocks(mask);
            
            W = t.codomain.new(dims, false);
            N = t.eye(W, V.domain);
            N.var = fill_matrix_data(N.var, Ns, dims.charges);
%             
%             [~, S, V] = tpermute(t, [p1 p2], [length(p1) length(p2)]);
%             
%             dims = struct;
%             [mblocks, dims.charges] = matrixblocks(t);
%             Ns = cell(size(mblocks));
%             dims.degeneracies = zeros(size(mblocks));
%             
%             for i = 1:length(mblocks)
%                 Ns{i} = rightnull(mblocks{i}, alg, atol);
%                 dims.degeneracies(i) = size(Ns{i}, 1);
%             end
%             
%             mask = dims.degeneracies > 0;
%             dims.charges = dims.charges(mask);
%             dims.degeneracies = dims.degeneracies(mask);
%             Ns = Ns(mask);
%             
%             N = Tensor.eye(t.domain.new(dims, false), t.domain);
%             N.var = fill_matrix_data(N.var, Ns, dims.charges);
        end
        
        function [U, S, V, eta] = tsvd(t, p1, p2, trunc)
            % Compute the singular value decomposition of a tensor. This computes left and
            % right isometries U and V, and a non-negative diagonal tensor S such that
            % :code:`norm(tpermute(t, [p1 p2], [length(p1) length(p2)]) - U * S * V) = 0`
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
            %   truncate such that the dim of S is not larger than this value for any given
            %   charge.
            %
            % TruncTotalDim : int
            %   truncate such that the total dim of S is not larger than this value.
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
            %   :code:`U * S * V = tpermute(t, [p1 p2], [length(p1) length(p2)])`.
            %
            % eta : numeric
            %   truncation error.
            
            arguments
                t
                p1 = 1:t.rank(1)
                p2 = t.rank(1) + (1:t.rank(2))
                trunc.TruncDim
                trunc.TruncTotalDim
                trunc.TruncBelow
                trunc.TruncSpace
                trunc.doPlot = false
            end
            
            t = tpermute(t, [p1 p2], [length(p1) length(p2)]);
            
            dims = struct;
            [mblocks, dims.charges] = matrixblocks(t);
            Us = cell(size(mblocks));
            Ss = cell(size(mblocks));
            Vs = cell(size(mblocks));
            dims.degeneracies = zeros(size(mblocks));
            
            doTrunc = ~isempty(fieldnames(trunc));
            eta = 0;
            for i = 1:length(mblocks)
                if doTrunc
                    [Us{i}, Ss{i}, Vs{i}] = svd(mblocks{i}, 'econ');
                    dims.degeneracies(i) = size(Ss{i}, 1);
                else
                    [Us{i}, Ss{i}, Vs{i}] = svd(mblocks{i});
                end
                Vs{i} = Vs{i}';
            end

            if isfield(trunc, 'doPlot') && trunc.doPlot
                Ss_old = Ss;
                W_old = t.domain.new(dims, false);
                S_old = Tensor.zeros(W_old, W_old);
                S_old.var = fill_matrix_data(S_old.var, Ss_old, dims.charges);
            end
            
            if isfield(trunc, 'TruncBelow')
                for i = 1:length(mblocks)
                    s = diag(Ss{i});
                    dims.degeneracies(i) = sum(s > trunc.TruncBelow);
                    eta = eta + sum(s(dims.degeneracies(i) + 1:end).^2 * qdim(dims.charges(i)));
                    Us{i} = Us{i}(:, 1:dims.degeneracies(i));
                    Ss{i} = diag(s(1:dims.degeneracies(i)));
                    Vs{i} = Vs{i}(1:dims.degeneracies(i), :);
                end
            end
            if isfield(trunc, 'TruncTotalDim')
                qdims = qdim(dims.charges);
                totaldim = sum(dims.degeneracies .* qdims);
                minvals = cellfun(@(x) min(diag(x)), Ss);
                while totaldim > trunc.TruncTotalDim
                    [~, i] = min(minvals);
                    eta = eta + Ss{i}(end, end)^2 * qdims(i);
                    dims.degeneracies(i) = dims.degeneracies(i) - 1;
                    Us{i} = Us{i}(:, 1:end-1);
                    Vs{i} = Vs{i}(1:end-1, :);
                    Ss{i} = Ss{i}(1:end-1, 1:end-1);
                    if dims.degeneracies(i) > 0
                        minvals(i) = Ss{i}(end,end);
                    else
                        minvals(i) = Inf;
                    end
                    totaldim = totaldim - qdims(i);
                end
            end
            if isfield(trunc, 'TruncDim')
                for i = 1:length(mblocks)
                    dims.degeneracies(i) = min(dims.degeneracies(i), trunc.TruncDim);
                    s = diag(Ss{i});
                    eta = eta + sum(s(dims.degeneracies(i) + 1:end).^2 * qdim(dims.charges(i)));
                    Us{i} = Us{i}(:, 1:dims.degeneracies(i));
                    Ss{i} = diag(s(1:dims.degeneracies(i)));
                    Vs{i} = Vs{i}(1:dims.degeneracies(i), :);
                end
            end
            if isfield(trunc, 'TruncSpace')
                truncspace = trunc.TruncSpace;
                [b, ind] = ismember(truncspace.dimensions.charges,dims.charges);
                assert(all(b),"Truncation space contains charges that S does not.")
                assert(all(dims.degeneracies(ind)>=truncspace.dimensions.degeneracies), "Truncation space has degeneracies larger than S.")
                dims.degeneracies(ind) = truncspace.dimensions.degeneracies;
                dims.degeneracies(setxor(ind,1:numel(dims.degeneracies))) = 0;
                for i = 1:length(mblocks)
                    s = diag(Ss{i});
                    eta = eta + sum(s(dims.degeneracies(i) + 1:end).^2 * qdim(dims.charges(i)));
                    Us{i} = Us{i}(:, 1:dims.degeneracies(i));
                    Ss{i} = diag(s(1:dims.degeneracies(i)));
                    Vs{i} = Vs{i}(1:dims.degeneracies(i), :);
                end
            end
            
            if isfield(trunc, 'doPlot') && trunc.doPlot
                ax = gca;
                ax = plot_schmidtspectrum(S_old, ax);
                hold on
                for i = 1:length(mblocks)
                    s = diag(Ss_old{i});
                    s(dims.degeneracies(i)+1:end) = 0;
                    Ss_old{i} = diag(s);
                end
                S_old.var = fill_matrix_data(S_old.var, Ss_old, dims.charges);
                ax = plot_schmidtspectrum(S_old, ax, 'Marker', 'o');
                hold off
            end

            if ~doTrunc
                W1 = prod(t.codomain);
                W2 = prod(t.domain);
                
                U = Tensor.eye(t.codomain, W1);
                S = Tensor.zeros(W1, W2);
                V = Tensor.eye(W2, t.domain);
            else
                mask = dims.degeneracies ~= 0;
                dims.charges = dims.charges(mask);
                dims.degeneracies = dims.degeneracies(mask);
                W = t.domain.new(dims, false);
                
                U = Tensor.eye(t.codomain, W);
                S = Tensor.zeros(W, W);
                V = Tensor.eye(W, t.domain);
                
                Us = Us(mask);
                Vs = Vs(mask);
                Ss = Ss(mask);
            end
            
            U.var = fill_matrix_data(U.var, Us, dims.charges);
            S.var = fill_matrix_data(S.var, Ss, dims.charges);
            V.var = fill_matrix_data(V.var, Vs, dims.charges);

            if nargout <= 1
                U = S;
            end
            eta = sqrt(eta);
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
            
            assert(isisometric(t.codomain, t.domain), 'tensors:ArgumentError', ...
                'Input should be square.');
            
            mblocks = matrixblocks(t);
            for i = 1:length(mblocks)
                mblocks{i} = inv(mblocks{i});
            end

            if isequal(t.codomain, t.domain)
                t.var = fill_matrix_data(t.var, mblocks);
            else
                t = t';
                t.var = fill_matrix_data(t.var, mblocks);
            end
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
            
            if numel(t) > 1
                bool = arrayfun(@(x) isisometry(x, side, ...
                    'RelTol', tol.RelTol, 'AbsTol', tol.AbsTol), t);
                return
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
            bool = true;
        end
        
        function bool = iszero(t)
            bool = isempty(t.var) || iszero(t.var);
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
        function v = vectorize(t, type)
            % Collect all parameters in a vector, weighted to reproduce the correct
            % inner product.
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
            
            v = vectorize([t.var], type);
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
            
            if numel(t) > 1
                X = devectorize(v, [t.var], type);
                for i = 1:numel(t)
                    t(i).var = X(i);
                end
                return
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
            
            if isa(t.var, 'TrivialBlock')
                blocks = tensorblocks(t);
                a = blocks{1};
                return
            end
            
            [blocks, trees] = tensorblocks(t);
            
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
        
        function a = sparse(t)
            a = SparseTensor(t);
        end
        
        function tdst = desymmetrize(tsrc, mode)
            arguments
                tsrc
                mode {mustBeMember(mode, {'Cartesian', 'Complex'})} = 'Complex'
            end
            
            tdst = repartition(Tensor(double(tsrc)), rank(tsrc));
            if strcmp(mode, 'Complex')
                if ~isempty(tsrc.domain), tdst.domain = ComplexSpace(tsrc.domain); end
                if ~isempty(tsrc.codomain), tdst.codomain = ComplexSpace(tsrc.codomain); end
            end
        end
        
        function mps = FiniteMps(t, varargin)
            A = MpsTensor.decompose_local_state(t, varargin{:});
            mps = FiniteMps(A);
        end
    end
    
    
    %% Utility
    methods
        function [I, J, V] = find(t, k, which)
            arguments
                t
                k = []
                which = 'first'
            end
            assert(strcmp(which, 'first'), 'not implemented yet')
            if ~isempty(k)
                assert(k <= numel(t));
            else
                k = numel(t);
            end
            
            if isempty(t)
                I = [];
                J = [];
                V = [];
                return
            end
            
            if ~isempty(k)
                if strcmp(which, 'first')
                    I = 1:k;
                else
                    I = numel(t):-1:numel(t)+1-k;
                end
            end
            I = reshape(I, [], 1);
            if nargout < 2, return; end
            
            sz = size(t);
            subs = ind2sub_([sz(1) prod(sz(2:end))], I);
            I = subs(:, 1);
            J = subs(:, 2);
            V = reshape(t, [], 1);
        end
        
        function s = GetMD5_helper(t)
            s = {t.codomain t.domain};
        end
        
        function type = underlyingType(t)
            type = underlyingType(t(1).var);
        end
        
        function [codomain, domain] = deduce_spaces(t)
            spaces = cell(1, nspaces(t));
            subs = repmat({1}, 1, nspaces(t));
            for i = 1:length(spaces)
                for j = flip(1:size(t, i))
                    subs{i} = j;
                    spaces{i}(j) = space(t(subs{:}), i);
                end
                subs{i} = 1;
            end
            Nout = indout(t);
            if Nout > 0
                codomain = SumSpace(spaces{1:Nout});
            else
                codomain = SumSpace([]);
            end
            if Nout == length(spaces)
                domain = SumSpace([]);
            else
                domain = SumSpace(spaces{(Nout+1):end})';
            end
        end
        
        function disp(t, details)
            if nargin == 1 || isempty(details), details = false; end
            if isscalar(t)
                r = t.rank;
                fprintf('Rank (%d, %d) %s:\n', r(1), r(2), class(t));
                s = space(t);
                for i = 1:length(s)
                    fprintf('\t%d.\t', i);
                    disp(s(i));
                    fprintf('\b');
                end
                fprintf('\n');
                if details
                    [blocks, charges] = matrixblocks(t);
                    for i = 1:length(blocks)
                        if ~isempty(blocks)
                            fprintf('charge %s:\n', string(charges(i)));
                        end
                        disp(blocks{i});
                    end
                end
            else
                fprintf('%s of size %s:\n', class(t), ...
                    regexprep(mat2str(size(t)), {'\[', '\]', '\s+'}, {'', '', 'x'}));
                subs = ind2sub_(size(t), 1:numel(t));
                spc = floor(log10(max(double(subs), [], 1))) + 1;
                if numel(spc) == 1
                    fmt = strcat("\t(%", num2str(spc(1)), "u)");
                else
                    fmt = strcat("\t(%", num2str(spc(1)), "u,");
                    for i = 2:numel(spc) - 1
                        fmt = strcat(fmt, "%", num2str(spc(i)), "u,");
                    end
                    fmt = strcat(fmt, "%", num2str(spc(end)), "u)");
                end
                for i = 1:numel(t)
                    fprintf('%s\t\t', compose(fmt, subs(i, :)));
                    disp(t(i), details);
                end
            end
        end
    end
end
