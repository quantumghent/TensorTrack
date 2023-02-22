classdef AbstractTensor
    
    
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
            
            if issparse(b), b = full(b); end
            
            % Convert input objects to vectors
            b_vec = vectorize(b);
            b_sz = size(b_vec);
            
            if ~isempty(x0) && ~issparse(x0)
                x0_vec = vectorize(x0);
            else
                x0_vec = zeros(b_sz);
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
                
                options.Tol = eps(underlyingType(x0))^(3/4)
                options.Algorithm = 'eigs'
                options.MaxIter = 100
                options.KrylovDim = 20
                options.IsSymmetric logical = false
                options.DeflateDim
                options.ReOrth = 2
                options.NoBuild
                options.Verbosity = 0
            end
            
            assert(isnumeric(sigma) || ismember(sigma, {'largestabs', 'smallestabs', ...
                'largestreal', 'smallestreal', 'bothendsreal', ...
                'largestimag', 'smallestimag', 'bothendsimag'}), ...
                'tensors:ArgumentError', 'Invalid choice of eigenvalue selector.');
            nargoutchk(0, 3);
            
            
            x0_vec = vectorize(x0);
            sz = size(x0_vec);
            
            if isa(A, 'function_handle')
                A_fun = @(x) vectorize(A(devectorize(x, x0)));
            else
                A_fun = @(x) vectorize(A * devectorize(x, x0));
            end
            
            options.KrylovDim = min(sz(1), options.KrylovDim);
            
            [V, D, flag] = eigs(A_fun, sz(1), howmany, sigma, ...
                'Tolerance', options.Tol, 'MaxIterations', options.MaxIter, ...
                'SubspaceDimension', options.KrylovDim, 'IsFunctionSymmetric', ...
                options.IsSymmetric, 'StartVector', x0_vec, ...
                'Display', options.Verbosity >= 3);
            
            if nargout <= 1
                varargout = {D};
            elseif nargout == 2
                for i = howmany:-1:1
                    varargout{1}(:, i) = devectorize(V(:, i), x0);
                end
                varargout{2} = D;
            else
                varargout{2} = D;
                for i = howmany:-1:1
                    varargout{1}(:, i) = devectorize(V(:, i), x0);
                end
                varargout{3} = flag;
            end
            
            if flag && options.Verbosity > 0
                warning('eigsolve did not converge.');
            elseif ~flag && options.Verbosity > 1
                fprintf('eigsolve converged.\n');
            end
        end
        
        function C = contract(tensors, indices, kwargs)
            arguments (Repeating)
                tensors {mustBeA(tensors, 'AbstractTensor')}
                indices (1, :) {mustBeInteger}
            end
            
            arguments
                kwargs.Conj (1, :) logical = false(size(tensors))
                kwargs.Rank = []
                kwargs.Debug = false
                kwargs.CheckOptimal = false
            end
            
            assert(length(kwargs.Conj) == length(tensors));
            
            if kwargs.CheckOptimal
                legcosts = zeros(2, 0);
                for i = 1:length(indices)
                    legcosts = [legcosts [indices{i}; dims(tensors{i})]];
                end
                legcosts = unique(legcosts.', 'rows');
                
                currentcost = contractcost(indices, legcosts);
                [sequence, cost] = netcon(indices, 1, 1, currentcost, 1, legcosts);
                
                if cost < currentcost
                    warning('suboptimal contraction order.\n optimal: %s', ...
                        num2str(sequence));
                end
            end
            
            for i = 1:length(tensors)
                [i1, i2] = traceinds(indices{i});
                tensors{i} = tensortrace(tensors{i}, i1, i2);
                indices{i}([i1 i2]) = [];
            end
            
            debug = kwargs.Debug;
            
            % Special case for single input tensor
            if nargin == 2
                C = tensors{1};
                
                if isempty(indices{1})
                    assert(isnumeric(C));
                    if kwargs.Conj
                        C = C';
                    end
                    return
                end
                
                [~, order] = sort(indices{1}, 'descend');
                
                if kwargs.Conj
                    C = tpermute(C', order(length(order):-1:1), kwargs.Rank);
                else
                    C = tpermute(C, order, kwargs.Rank);
                end
                return
            end
            
            % Generate trees
            contractindices = cellfun(@(x) x(x > 0), indices, 'UniformOutput', false);
            partialtrees = num2cell(1:length(tensors));
            tree = generatetree(partialtrees, contractindices);
            
            % contract all subtrees
            [A, ia, ca] = contracttree(tensors, indices, kwargs.Conj, tree{1}, debug);
            [B, ib, cb] = contracttree(tensors, indices, kwargs.Conj, tree{2}, debug);
            
            % contract last pair
            [dimA, dimB] = contractinds(ia, ib);
            
            if debug, contractcheck(A, ia, ca, B, ib, cb); end
            
            C = tensorprod(A, B, dimA, dimB, ca, cb, 'NumDimensionsA', length(ia));
            ia(dimA) = [];  ib(dimB) = [];
            ic = [ia ib];
            
            % permute last tensor
            if ~isempty(ic) && length(ic) > 1
                [~, order] = sort(ic, 'descend');
                if isempty(kwargs.Rank)
                    kwargs.Rank = [length(order) 0];
                end
                C = tpermute(C, order, kwargs.Rank);
            end
        end
        
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
            
            n = max(ndims(A), ndims(B));
            assert(isequal(size(A, 1:n), size(B, 1:n)) || isscalar(A) || isscalar(B), ...
                'tensors:SizeError', 'Incompatible sizes for vectorized function.');
            
            % make everything a vector
            A = repartition(A);
            B = repartition(B);
            
            d = norm(A - B);
        end
        
        function local_operators = decompose_local_operator(H, kwargs)
            % convert a tensor into a product of local operators.
            %
            % Usage
            % -----
            % :code:`local_operators = decompose_local_operator(H, kwargs)`.
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
            %   :method:`Tensor.tsvd`
            arguments
                H
                kwargs.Trunc = {'TruncBelow', 1e-14}
            end
            
            assert(mod(nspaces(H), 2) == 0, ...
                'InfJMpo:Argerror', 'local operator must have an even amount of legs.');
            H = repartition(H, nspaces(H) ./ [2 2]);
            assert(isequal(H.domain, H.codomain), ...
                'InfJMpo:ArgError', 'local operator must be square.');
            
            N = indin(H);
            local_operators = cell(1, N);
            if N == 1
                local_operators{1} = insert_onespace(insert_onespace(H, 1), 3);
            else
                [u, s, v] = tsvd(H, [1 N+1], [2:N N+2:2*N], kwargs.Trunc{:});
                local_operators{1} = insert_onespace(tpermute(u * s, [1 3 2], [1 2]), 1);
                
                for i = 2:N-1
                    [u, s, v] = tsvd(v, [1 2 N-i+3], [3:(N-i+2) (N-i+4):nspaces(v)], ...
                        kwargs.Trunc{:});
                    local_operators{i} = tpermute(u * s, [1 2 4 3], [2 2]);
                end
                
                local_operators{N} = insert_onespace(repartition(v, [2 1]), 3);
            end
        end
        
        function B = tensortrace(A, i1, i2)
            if isempty(i1) && isempty(i2), B = A; return; end
            assert(length(i1) == length(i2), 'invalid indices');
            
            E = A.eye(conj(space(A, i1)), space(A, i2));
            iA = [i1 i2];
            iE = [1:length(i1) length(i1) + (length(i2):-1:1)];
            B = tensorprod(A, E, iA, iE);
        end
        
        function sz = dims(t, inds)
            sz = dims([t.codomain, t.domain']);
            if nargin > 1
                sz = sz(inds);
            end
        end
    end
end

