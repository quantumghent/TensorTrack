classdef SparseArray
    % Class for multi-dimensional sparse arrays.
    %
    %   Limited to arrays with a total number of elements of at most 2^48-1.
    
    %% Properties
    properties (Access = private)
        var = sparse([]) % (:, 1) double
        sz = [] % (1, :)
    end
    
    
    %% Constructors
    methods        
        function a = SparseArray(varargin)
            % Create a sparse array.
            %
            % Usage
            % -----
            % :code:`a = SparseArray(subs, vals, sz)`
            %   uses the rows of :code:`subs` and :code:`vals` to generate a sparse array
            %   :code:`a` of size :code:`sz = [m1 m2 ... mn]`. :code:`subs` is a
            %   :code:`p` x :code:`n` array specifying the subscripts of the nonzero values
            %   to be inserted into :code:`a`. The k-th row of :code:`subs` specifies the
            %   subscripts for the k-th value in :code:`vals`.
            %   The argument :code:`vals` may be scalar, in which case it is expanded to be
            %   the same length as :code:`subs`, i.e., it is equivalent to
            %   :code:`vals * (p, 1)`.
            %   In the case of duplicate subscripts in :code:`subs`, the corresponding
            %   values are added.
            %
            % :code:`a = SparseArray`
            %   Empty constructor.
            %
            % :code:`a = SparseArray(b)`
            %   Copies/converts :code:`b` if it is a :class:`SparseArray`, a dense array or
            %   a sparse matrix.
            %
            % :code:`a = SparseArray(b, sz)`
            %   Copies/converts :code:`b` if it is a :class:`SparseArray`, a dense array or
            %   a sparse matrix, and sets the size of :code:`a` to :code:`sz`
            % 
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> subs = [1 1 1; 1 1 3; 2 2 2; 4 4 4; 1 1 1; 1 1 1];
            %   >> vals = [0.5; 1.5; 2.5; 3.5; 4.5; 5.5];
            %   >> sz = [4 4 4];
            %   >> a = SparseArray(subs, vals, sz) %<-- sparse 4x4x4
            
            if (nargin == 0) || ((nargin == 1) && isempty(varargin{1}))
                % Empty constructor
                return;

            elseif nargin == 1
                % Single argument

                source = varargin{1};

                switch(class(source))

                    % copy constructor
                    case 'SparseArray'
                        a.var = source.var;
                        a.sz = source.sz;
                                                
                    % sparse matrix, dense array
                    case {'numeric', 'logical', 'double'}
                        a.var = reshape(source, [], 1);
                        if ~issparse(source)
                            a.var = sparse(a.var);
                        end
                        a.sz = size(source);
                        
                    otherwise
                        error('sparse:UnsupportedConstructor', 'Unsupported use of SparseArray constructor.');

                end

                
            elseif nargin == 2
                % Two arguments
                
                var = varargin{1};
                sz = varargin{2}(:).';
                if isempty(sz)
                    sz = size(var);
                end
                if numel(sz) < 2
                    sz = [sz, 1];
                end
                if numel(var) ~= prod(sz)
                    error('sparse:ArgumentError', 'Incompatible size vector and input array.');
                end

                switch(class(var))

                    % copy constructor
                    case 'SparseArray'
                        a.var = var.var;
                        a.sz = sz;

                    % sparse matrix, dense array
                    case {'numeric', 'logical', 'double'}
                        a.var = reshape(var, [], 1);
                        if ~issparse(var)
                            a.var = sparse(a.var);
                        end
                        a.sz = sz;
                        
                    otherwise
                        error('sparse:UnsupportedConstructor', 'Unsupported use of SparseArray constructor.');

                end
                
                
            elseif nargin == 3
                % Three arguments
                subs = varargin{1};
                vals = varargin{2};
                sz = varargin{3}(:)';
                if numel(sz) < 2
                    sz = [sz, 1];
                end
                
                ind = sub2ind_(sz, subs);
                a.var = sparse(ind, ones(size(ind)), vals, prod(sz), 1);
                a.sz = sz;

            else
                error('sparse:UnsupportedConstructor', 'Unsupported use of SparseArray constructor.');
            end
            
        end
                
    end
    
    
    %% Methods
    methods
        
        function a = abs(a)
            % Absolute value.
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array.
            a.var = abs(a.var);
        end
        
        function a = chop(a, tol)
            % Set all nonzero values in :class:`SparseArray` who's absolute value is below
            % a given threshold to zero.
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % tol : :class:`float` , optional
            %   threshold tolerance for absolute values of entries, defaults to
            %   :code:`1e-15`.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   sparse array with entries of absolute value below :code:`tol` set to zero.
            arguments
                a
                tol = 1e-15
            end
            
            [r, c, v] = find(a.var);
            drop = abs(v) < tol;
            a.var(r(drop), c(drop)) = 0;
        end
        
        function a = conj(a)
            % Complex conjugate.
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array.
            a.var = conj(a.var);
        end
        
        function a = ctranspose(a)
            % Complex conjugate transpose.
            %
            % Only defined for 2D sparse arrays.
            %
            % Usage
            % -----
            % :code:`b = ctranspose(a)`
            %
            % :code:`b = a'`
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array.
            assert(ismatrix(a), 'sparse:RankError', 'ctranspose is only defined for 2D sparse arrays.');
            a = permute(conj(a), [2, 1]);
        end
                
        function d = diag(a)
            assert(ismatrix(a), 'sparse:RankError', 'diag is only defined for 2D sparse arrays.');
            d = diag(reshape(a.var, a.sz));
        end
        
        function disp(a)
            nz = nnz(a);

            if nz == 0
                fprintf('all-zero %s of size %s\n', class(a), ...
                    dim2str(a.sz));
                return
            end

            fprintf('%s of size %s with %d nonzeros:\n', class(a), ...
                dim2str(a.sz), nz);

            if (nz > 1000)
                r = input('Big array, do you want to print all nonzeros? (y/n) ', 's');
                if ~strcmpi(r, 'y'), return, end
            end

            [subs, ~, vals] = find(a);
            spc = floor(log10(max(double(subs), [], 1))) + 1;
            fmt_subs = sprintf('%%%du,', spc(1:end-1));
            fmt_subs = sprintf('%s%%%du', fmt_subs, spc(end));
            fmt = sprintf('\t(%s)%%s\n', fmt_subs);
            S = evalc('disp(vals)'); % abuse builtin display
            S = splitlines(S);
            if contains(S{1}, '*') % big numbers
                fprintf('%s\n', S{1});
                S = S(2:end);
            end
            for i = 1:nz
                fprintf(fmt, subs(i,:), S{i});
            end
        end
                
        function varargout = eig(a, varargin)
            % Find eigenvalues and eigenvectors of a 2D sparse array.
            %
            % See Also
            % --------
            % `Documentation for builtin Matlab eig <https://mathworks.com/help/matlab/ref/eig.html>`_.
            assert(ismatrix(a), 'sparse:RankError', 'eig is only defined for 2D sparse arrays.');
            [varargout{1:nargout}] = eig(reshape(a.var, a.sz), varargin{:});
            varargout = cellfun(@maybesparse_, varargout, 'UniformOutput', false);
        end
        
        function varargout = eigs(a, varargin)
            % Find a few eigenvalues and eigenvectors of a 2D sparse array.
            %
            % See Also
            % --------
            % `Documentation for builtin Matlab eigs <https://mathworks.com/help/matlab/ref/eigs.html>`_.
            assert(ismatrix(a), 'sparse:RankError', 'eigs is only defined for 2D sparse arrays.');
            [varargout{1:nargout}] = eigs(reshape(a.var, a.sz), varargin{:});
            varargout = cellfun(@maybesparse_, varargout, 'UniformOutput', false);
        end
        
        function [subs, idx, vals] = find(a)
            % Find subscripts of nonzero elements in a sparse array.
            %
            % Usage
            % -----
            % :code:`idx = find(a)`
            %
            % :code:`[subs, idx, vals] = find(a)`
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array
            %
            % Returns
            % -------
            % subs : (:, :) :class:`int`
            %   subscripts of nonzero array entries.
            %
            % idx : (:, 1) :class:`int`
            %   linear indices of nonzero array entries.
            %
            % var : (:, 1) :class:`double`
            %   values of nonzero array entries.
            [idx, ~, vals] = find(a.var);
            if nargout == 1
                subs = idx;
                return
            end
            subs = ind2sub_(a.sz, idx);
        end
        
        function d = full(a)
            % Convert a sparse array to a dense array.
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % d : :class:`double`
            %   dense output array.
            d = reshape(full(a.var), a.sz);
        end
        
        function a = groupind(a, g)
            % Group (block) specified sets of contiguous indices.
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            % g : (1, :) :class:`int`
            %   list of number of contiguous indices to be grouped in each index of the
            %   output tensor.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array with grouped indices.
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> a = SparseArray.random([2, 3, 4, 5, 6], .1);
            %   >> b = groupind(a, [3, 2]); %<-- sparse 24x30
            assert(sum(g) == ndims(a), 'sparse:InvalidGrouping', 'Invalid grouping.')
            r = zeros(1, length(g));
            offset = 0;
            for i = 1:length(g)
                r(i) = prod(a.sz(1+offset:g(i)+offset));
                offset = offset + g(i);
            end
            a = reshape(a, r);
        end
        
        function a = imag(a)
            % Complex imaginary part of sparse array.
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array with real entries corresponding to the imaginary part of the
            %   entries of :code:`a`.
            a.var = imag(a.var);
        end
        
        function bool = ismatrix(a)
            bool = ndims(a) == 2;
        end
        
        function bool = isnumeric(~)
            % Determine whether input is numeric.
            % 
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % bool : :class:`logical`
            %   defaults to :code:`true` for :class:`SparseArray`.
            bool = true;
        end
        
        function bool = isscalar(~)
            % Determine whether input is scalar.
            % 
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % bool : :class:`logical`
            %   defaults to :code:`false` for :class:`SparseArray`.
            bool = false;
        end
        
        function bool = isrow(a)
            % Determine whether input is row vector.
            % 
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % bool : :class:`logical`
            bool = ismatrix(a) && size(a, 1) == 1;
        end
        
        function bool = isstruct(~)
            bool = false;
        end
        
        function bool = istriu(a)
            bool = istriu(spmatrix(a));
        end
        
        function c = ldivide(a, b)
            % Elementwise left division for sparse arrays.
            %
            % :code:`ldivide(a, b)` is called for the syntax :code:`a .\ b` where :code:`a`
            % or :code:`b` is a :class:`SparseArray`. :code:`a`and :code:`b` must have the
            % same size, unless one is a scalar. 
            %
            % Arguments
            % ---------
            % a, b : :class:`SparseArray` or :class:`double`
            %   input arrays to be divided.
            %
            % Returns
            % -------
            % c : :class:`SparseArray`
            %   elementwise left division of :code:`b` by :code:`a`.
            c = rdivide(b, a);
        end
        
        function c = minus(a, b)
            % Elementwise subtraction for sparse arrays. 
            %
            % :code:`minus(a, b)` is called for the syntax :code:`a - b` where :code:`a`
            % or :code:`b` is a :class:`SparseArray`. :code:`a`and :code:`b` must have the
            % same size, unless one of the is scalar. Scalar can be subtracted from a sparse
            % array of any size, resulting in a dense array.
            %
            % Arguments
            % ---------
            % a, b : :class:`SparseArray` or :class:`double`
            %   intput arrays.
            %
            % Returns
            % -------
            % c : :class:`SparseArray` or :class:`double`
            %   output array.
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> a = SparseArray.random([4 3 2], .1);
            %   >> b = SparseArray.random([4 3 2], .1);
            %   >> a - b %<-- sparse
            %   >> a - 5 %<-- dense
            %   >> a - 0 %<-- dense
            %   >> a - full(a) %<-- dense
            c = plus(a, -b);
        end
        
        function c = mrdivide(a, b)
            % Matrix right division for sparse arrays.
            %
            % Usage
            % -----
            % :code:`mrdivide(a, b)`
            %
            % :code:`a / b`
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   intput array.
            % b : :class:`double`
            %   scalar to divide by.
            %
            % Returns
            % -------
            % c : :class:`SparseArray`
            %   output array.
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> a = SparseArray.random([4 3 2], .1);
            %   >> a / 3 %<-- sparse
            assert(isscalar(b), 'sparse:ScalarDivide', 'SparseArray.mrdivide only supports the scalar case.')
            c = a;
            c.var = c.var/b;
            return;
        end
        
        function c = mtimes(a, b)
            % Matrix multiplication for 2D sparse arrays.
            %
            % :code:`mtimes(a, b)` is called for the syntax :code:`a * b` where :code:`a`
            % or :code:`b` is a :class:`SparseArray`.
            %
            % See Also
            % --------
            % `mtimes <https://mathworks.com/help/matlab/ref/mtimes.html>`_
            if isscalar(a)
                c = b;
                c.var = a * c.var;
                return
            elseif isscalar(b)
                c = a;
                c.var = b * c.var;
                return
            end

            sz = [size(a,1) size(b,2)];

            if isa(a, 'SparseArray')
                a = reshape(a.var, a.sz);
            end

            if isa(b, 'SparseArray')
                b = reshape(b.var, b.sz);
            end

            % always return sparse array; TODO: make sure this doesn't break things
            c = SparseArray(a * b, sz);
        end
        
        function n = ndims(t)
            % Number of dimensions of a sparse array.
            %
            % Note that unlike the `builtin Matlab behavior <https://mathworks.com/help/matlab/ref/double.ndims.html>`
            % trailing singleton dimensions are not ignored.
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> a = SparseArray.random([4 3 1], .1);
            %   >> ndims(a) %<-- returns 3
            n = size(t.sz, 2);
        end
        
        function n = nnz(a)
            % Number of nonzero elements in sparse array.
            n = nnz(a.var);
        end
        
        function nrm = norm(a)
            % Frobenius norm of a sparse array.
            nrm = norm(a.var);
        end
        
        function n = numArgumentsFromSubscript(varargin)
            n = 1;
        end

        function n = numel(a)
            % Number of elements in a sparse array.
            n = prod(size(a));
        end
        
        function c = outer(a, b)
            % Outer product of two sparse arrays.
            %
            % Warning
            % -------
            % Not :code:`kron`.
            varC = kron(b.var, a.var); % ordering is counterintuitive here
            c = SparseArray(varC, [size(a), size(b)]);
        end
        
        function a = permute(a, order)
            % Permute sparse array dimensions.
            if length(a.sz) <  length(order)
                a = reshape(a, [a.size, ones(1,length(order)-length(a.size))]);
            end
            if issorted(order) % don't rebuild t if it was a trivial permute
                return
            end
            sz1 = a.sz;
            sz2 = sz1(order);
            [idx, ~, v] = find(a.var);
            idx2 = sub2ind_(sz2, ind2sub_(sz1, idx, order));
            a = SparseArray(sparse(idx2, ones(size(idx2)), v, numel(a.var), 1), sz2);
        end
        
        function c = plus(a, b)
            % Elementwise addition for sparse arrays. 
            %
            % Usage
            % -----
            % 
            % :code:`plus(a, b)`
            %
            % :code:`a + b`
            %
            % :code:`a` and :code:`b` must have the same size,
            % unless one is a scalar. A scalar can be added to a sparse array of any size.   
            %
            % Arguments
            % ---------
            % a, b : :class:`SparseArray` or :class:`double`
            %   input arrays.
            %
            % Returns
            % -------
            % 
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> a = SparseArray.random([4 3 2], .1);
            %   >> a + b %<-- sparse
            %   >> a + 5 %<-- dense
            %   >> a + 0 %<-- dense
            %   >> a + full(a) %<-- dense
            if (isa(a, 'double') && ~issparse(a)) || (isa(b, 'double') && ~issparse(b))
                c = full(a) + full(b);
                return;
            end
            c = SparseArray(a); b = SparseArray(b);
            c.var = c.var + b.var;
        end
        
        function a = power(a, b)
            % Elementwise power for sparse array.
            assert(isscalar(b), 'sparse:NonScalarPower', 'SparseArray only supports power with a scalar.')
            a.var = a.var.^b;
        end
        
        function varargout = qr(a, varargin)
            % Orthogonal-triangular decomposition for a sparse array.
            %
            % See Also
            % --------
            % `Documentation for builtin Matlab qr <https://mathworks.com/help/matlab/ref/qr.html>`_.
            assert(ismatrix(a), 'sparse:RankError', 'qr is only defined for 2D sparse arrays.');
            [varargout{1:nargout}] = qr(reshape(a.var, a.sz), varargin{:});
            varargout = cellfun(@maybesparse_, varargout, 'UniformOutput', false);
        end
        
        function [Q, R] = qrpos(a, varargin)
            % Positive orthogonal-triangular decomposition for a sparse array.
            assert(ismatrix(a), 'sparse:RankError', 'qrpos is only defined for 2D sparse arrays.');
            [Q, R] = qr(a, varargin{:});
            if isrow(Q)
                Q = Q * sign(R(1));
                R = R * sign(R(1));
            else
                D = diag(R);
                D(abs(D) < 1e-10) = 1;
                D = sign(D);
                Q = Q * diag(D); % TODO: implement broadcasting for .*
                R = diag(D) * R;
            end
        end
        
        function c = rdivide(a, b)
            % Elementwise right division for sparse arrays.
            %
            % :code:`rdivide(a, b)` is called for the syntax :code:`a ./ b` where :code:`a`
            % or :code:`b` is a :class:`SparseArray`. :code:`a`and :code:`b` must have the
            % same size, unless one is a scalar. 
            %
            % Arguments
            % ---------
            % a, b : :class:`SparseArray` or :class:`double`
            %   input arrays to be divided.
            %
            % Returns
            % -------
            % c : :class:`SparseArray`
            %   elementwise left division of :code:`a` by :code:`b`.
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> a = SparseArray.random([4 3 2], .1);
            %   >> a ./ 5 %<-- sparse
            %   >> 5 ./ a %<-- dense
            %   >> a ./ full(a) %<-- sparse
            %   >> full(a) ./ a %<-- dense
            if isa(a, 'SparseArray') && isa(b, 'SparseArray')
                c = a;
                c.var = c.var ./ b.var;
            elseif isa(a, 'SparseArray')
                c = a.var ./ b;
            elseif isa(b, 'SparseArray')
                c = a ./ b.var;
            end
        end
        
        function a = real(a)
            % Complex real part of sparse array.
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array with real entries corresponding to the real part of the
            %   entries of :code:`a`.
            a.var = real(a.var);
        end
        
        function a = reshape(a, varargin)
            % Reshape sparse array.
            if nargin == 2
                a.sz = varargin{1};
            else
                sz = ones(1, numel(varargin));
                I = [];
                for i = 1:numel(varargin)
                    if isempty(varargin{i})
                        I = i;
                    else
                        sz(i) = varargin{i};
                    end
                end
                sz(I) = prod(a.sz) / prod(sz);
                a.sz = sz;
            end
        end
        
        function a = sign(a)
            % Signum function.
            a.var = sign(a.var);
        end

        function varargout = size(a, i)
            if nargin == 1
                sz = a.sz;
            else
                sz = ones(1, max(i));
                sz(1:length(a.sz)) = a.sz;
                sz = sz(i);
            end
            
            if nargout <= 1
                varargout = {sz};
            else
                varargout = num2cell(sz);
            end
        end
        
        function b = sparse(a)
            % Convert 2D sparse array to a sparse matrix.
            assert(ismatrix(a), 'sparse:NotMatrix', sprintf('Cannot convert array of rank %d to sparse matrix.', ndims(a)))
            b = reshape(a.var, a.size);
        end
        
        function b = spmatrix(a)
            % Convert 2D sparse array to a sparse matrix.
            %
            % See Also
            % --------
            % :code:`sparse`
            b = sparse(a);
        end
        
        function a = squeeze(a)
            % Remove singleton dimensions from a sparse array.
            %
            % Usage
            % -----
            % :code:`b = squeeze(a)` 
            %   returns a sparse array :code:`b` with the same elements as :code:`a` but
            %   with all the singleton dimensions removed.
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> squeeze(SparseArray.random([2, 1, 3], 0.5)) %<-- returns a 2 x 3 SparseArray
            %   >> squeeze(SparseArray([1, 1, 1], 1, [1, 1, 1])) %<-- returns a scalar
            if sum(a.sz > 1) == 0
                a = full(a.var);
                return
            end
            % always give n x 1 SparseArray in case of only 1 non-singleton dimension,
            % consistent with class constructor
            a.sz = [a.size(a.size>1), ones(1, 2-sum(a.size>1))];
        end
        
        function a = subsasgn(a, s, rhs)
            % Subscripted assignment for sparse array.
            assert(strcmp(s(1).type, '()'), 'sparse:index', 'only () indexing allowed');
            
            if length(s(1).subs) > 1 % non-linear indexing
                assert(length(s(1).subs) == size(a.sz, 2), 'sparse:index', ...
                    'number of indexing indices must match tensor size.');
                assert(all(a.sz >= cellfun(@max, s(1).subs)), 'sparse:bounds', ...
                    'out of bounds assignment disallowed');
                s(1).subs = {sub2ind(a.sz, s(1).subs{:})};
            end
            
            rhs = full(rhs);
            
            [I, ~, V] = find(a.var);
            [lia, locb] = ismember(s(1).subs{1}, I);
            newI = vertcat(I, s(1).subs{1}(~lia));
            newJ = ones(size(newI));
            V(locb(lia)) = rhs(lia);
            newV = vertcat(V, rhs(~lia));
            
            a.var = sparse(newI, newJ, newV, ...
                size(a.var, 1), size(a.var, 2));
        end
        
        function a_sub = subsref(a, s)
            % Subscripted reference for a sparse array.
            %
            % Usage
            % -----
            % :code:`a_sub = a(i)`
            %   linear indexing for SparseArray with only one non-singleton dimension,
            %   returns a scalar.
            %
            % :code:`a_sub = a(i1, i2, ..., iN)`
            %   where each :code:`in` is an integer index, returns a scalar.
            %
            % :code:`a_sub = a(R1, R2, ..., RN)`
            %   where each :code:`Rn` is either a colon, an array of integers representing a
            %   slice in this dimension or a logical array representing the logical indexing
            %   of this dimension. Returns a sparse array.
            %
            % :code:`a_sub = a(S)`
            %   where :code:`S` is a :code:`1` x :code:`n` array of integer subscripts (for
            %   dimensions that are indexed) or zeroes (for dimensions that are not
            %   indexed). Returns a sparse array.
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> a = SparseArray([4, 4, 4; 2, 2, 1; 2, 3, 2], [3; 5; 1], [4, 4, 4]);
            %   >> a(1, 2, 1) %<-- returns zero
            %   >> a(4, 4, 4) %<-- returns 3
            %   >> a(2, :, :) %<-- returns a 1 x 4 x 4 sparse array
            %   >> a([0, 4, 0]) %<-- returns a 4 x 1 x 4 sparse array
            %   >> a([4, 4, 0]) %<-- returns a 1 x 1 x 4 sparse array
            switch s(1).type

                case '()' % regular subscript indexing

                    % linear subscript indexing with only one non-singleton dimension
                    if (numel(s(1).subs) == 1 && numel(s(1).subs{1}) == 1 && sum(size(a)~=1) == 1)
                        subs = a.size;
                        subs(subs~=1) = s(1).subs{1};
                        a_sub = full(a.var(sub2ind_(a.size, subs)));

                    % regular multi-dimensional subscript indexing
                    elseif (numel(s(1).subs) == ndims(a))

                        % xtract the subdimensions to be extracted from t
                        sbs = s(1).subs;

                        % Error check that range is valid
                        okcolon = cellfun(@(x)  ischar(x) && x == ':', sbs);
                        okrange = cellfun(@(x)  isnumeric(x) && isreal(x) && ...
                                                ~any(isnan(x)) && ~any(isinf(x)) && ...
                                                isequal(x,round(x)) && all(x > 0), sbs);
                        oklogical = arrayfun(@(i)   islogical(sbs{i}) && ...
                                                    length(sbs{i}) == a.size(i), 1:ndims(a));
                        if ~all(okcolon | okrange | oklogical)
                            error('Invalid subscript.');
                        end

                        % simple case: all dimensions indexed
                        indexed = cellfun(@isnumeric, sbs) & cellfun(@(x) length(x) == 1, sbs);
                        if all(indexed)
                            a_sub = full(a.var(sub2ind_(a.size, [sbs{:}])));
                            return;
                        end

                        % convert logical indexing to integer ranges for easy handling
                        for i = find(oklogical)
                            sbs{i} = find(sbs{i});
                        end

                        % extract subscripts and values of nonzero elements
                        [subs, ~, nz] = find(a);

                        % some trickery to do actual slicing
                        f = true(size(nz)); % auxiliary array of booleans
                        new_sz = a.sz; % initialize size of sliced tensor
                        for i = fliplr(1:ndims(a))
                            if strcmp(sbs{i}, ':')
                                continue; % do nothing if not sliced
                            end
                            A = sbs{i};
                            if length(A) ~= length(unique(A))
                                error('Repeated index in position %i.', i);
                            end
                            if ~isempty(subs)
                                % fast intersect for this dimension
                                B = subs(:, i);
                                P = false(max(max(A),max(B))+1,1);
                                P(A+1) = true;
                                % throw away anything that has no intersect with the slice in this dimension
                                f = and(f, P(B+1));
                                % relabel subscripts in this dimension to match slice
                                % order
                                [~, ~, temp] = unique([A(:); subs(f, i)], 'stable');
                                subs(f, i) = temp(length(A)+1:end);
                                % adjust size in this dimension
                                new_sz(i) = length(A);
                            end
                        end
                        if isempty(subs)
                            a_sub = SparseArray([], [], new_sz);
                        else
                            a_sub = SparseArray(subs(f, :), nz(f), new_sz);
                        end

                        return;

                    % hacky indexing using an array
                    elseif (numel(s(1).subs) == 1) && (numel(s(1).subs{1}) == ndims(a))

                        % extract the subdimensions to be extracted from t
                        sbs = s(1).subs{1};

                        % error check that range is valid: each subscript must be a nonnegative integer
                        if ~ (isnumeric(sbs) && isreal(sbs) && ~any(isnan(sbs)) && ~any(isinf(sbs)) ...
                                && all(sbs >= 0))
                            error('Invalid subscript.');
                        end

                        indexed = sbs ~= 0;

                        if all(indexed) % if all dimensions are indexed with integer, return scalar
                            a_sub = full(a.var(sub2ind_(a.size, sbs)));
                            return;
                        end

                        new_sz = a.size;
                        new_sz(indexed) = 1; % set indexed dimensions to 1

                        [subs, ~, nz] = find(a);
                        good_ind = all(subs(:, indexed) == sbs(indexed), 2);
                        good_subs = subs(good_ind, :);
                        good_vals = nz(good_ind);
                        good_subs(:, indexed) = 1; % set indexed subscripts to 1 for new subs

                        a_sub = SparseArray(good_subs, good_vals, new_sz); % construct output
                        return;

                    else
                        error('sparse:InvalidIndexing', 'Incorrect indexing into SparseArray.')
                    end

                otherwise
                    % only overload parentheses indexing
                    a_sub = builtin('subsref', a, s);
            end
        end
        
        function s = sum(a)
            % Sum of all elements of a SparseArray.
            s = full(sum(a.var));
        end
        
        function varargout = svd(a, varargin)
            % Singular value decomposition.
            error('Not supported for sparse arrays.')
        end
        
        function varargout = svds(a, varargin)
            % Find a few singular values and vectors.
            % See Also
            % --------
            % `Documentation for builtin Matlab svds <https://mathworks.com/help/matlab/ref/svds.html>`_.
            assert(ismatrix(a), 'sparse:RankError', 'svds is only defined for 2D sparse arrays.');
            [varargout{1:nargout}] = svds(reshape(a.var, a.sz), varargin{:});
            varargout = cellfun(@maybesparse_, varargout, 'UniformOutput', false);
        end
        
        function c = times(a, b)
            % Array multiplication for sparse tensors.
            %
            % :code:`c = times(a, b)` is called for the syntax :code:`a .* b` when :code:`a`
            % or :code:`b` is a sparse array. :code:`a` and :code:`b` must have the same
            % size, unless one is a scalar. A scalar can be be multiplied by a sparse array
            % of any size.
            %
            % Example
            % -------
            % .. code-block:: matlab
            %
            %   >> a = SparseArray.random([4 3 2], .1);
            %   >> a .* b %<-- sparse
            %   >> a .* 5 %<-- sparse
            %   >> a .* 0 %<-- sparse
            %   >> a .* full(a) %<-- sparse
            if isscalar(b) || ~isa(b,'SparseArray')
                c = SparseArray(a.var .* b(:), a.sz);
            elseif isscalar(a) || ~isa(a,'SparseArray')
                c = SparseArray(b.var .* a(:), b.sz);
            else
                c = SparseArray(b.var .* a.var, b.sz);
            end
        end
        
        function a = transpose(a)
            % Transpose.
            %
            % Only defined for 2D sparse arrays.
            %
            % Usage
            % -----
            % :code:`b = transpose(a)`
            %
            % :code:`b = a.'`
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array.
            assert(ismatrix(a), 'sparse:RankError', 'ctranspose is only defined for 2D sparse arrays.');
            a = permute(a, [2, 1]);
        end

        function a = uminus(a)
            % Unary minus.
            %
            % Usage
            % -----
            % :code:`b = uminus(a)`
            %
            % :code:`b = -a`
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array.
            a.var = -a.var;
        end

        function a = uplus(a)
            % Unary plus.
            %
            % Usage
            % -----
            % :code:`b = uplus(a)`
            %
            % :code:`b = +a`
            %
            % Arguments
            % ---------
            % a : :class:`SparseArray`
            %   input array.
            %
            % Returns
            % -------
            % b : :class:`SparseArray`
            %   output array.
            return;
        end
        
    end
    
    methods (Static)
        function a = delta(numinds, inddim)
            % Create delta- (ghz-) array with given number of indices and index dimension.
            %
            % Arguments
            % ---------
            % numinds : :class:`int`
            %   number of indices of delta-array.
            %
            % inddim : :class:`int`
            %   dimension of each index of delta-array.
            %
            % Returns
            % -------
            % a : :class:`SparseArray`
            %   output delta-array.
            a = SparseArray(repmat(1:inddim, numinds, 1)', 1, repmat(inddim, 1, numinds));
        end
        
        function a = zeros(sz)
            a = SparseArray([], [], sz);
        end
        
        function a = random(sz, density)
            % Create a random complex sparse array.
            %
            % Arguments
            % ---------
            % sz : (1, :) :class:`int`
            %   size of the sparse array
            % density : :class:`double`
            %   density of nonzero elements (0 < :code:`density` < 1)
            a = SparseArray([], [], sz);
            a.var = sprandn(numel(a), 1, density);
            idx = find(a.var);
            a.var(idx) = a.var(idx) + 1i * randn(size(a.var(idx)));
        end
    end
    
end
