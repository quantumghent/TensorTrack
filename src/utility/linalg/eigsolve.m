function varargout = eigsolve(A, v, howmany, sigma, options)
% Find a few eigenvalues and eigenvectors of an operator.
%
% Usage
% -----
% :code:`[V, D, flag] = eigsolve(A, v, howmany, sigma, kwargs)`
%
% :code:`D = eigsolve(A, v, ...)`
%
% Arguments
% ---------
% A : :class:`matrix` or :class:`function_handle`
%   A square matrix.
%   A function handle which implements one of the following, depending on sigma:
%
%   - :code:`A \ x`, if `sigma` is 0 or 'smallestabs'
%   - :code:`(A - sigma * I) \ x`, if sigma is a nonzero scalar
%   - :code:`A * x`, for all other cases
%
% v : :class:`vector`
%   initial guess for the eigenvector. If A is a :class:`.Tensor`, this defaults
%   to a random complex :class:`.Tensor`, for function handles this is a required
%   argument.
%
% howmany : :class:`int`
%   amount of eigenvalues and eigenvectors that should be computed. By default
%   this is 1, and this should not be larger than the total dimension of A.
%
% sigma : :class:`char` or :class:`numeric`
%   selector for the eigenvalues, should be either one of the following:
%
%   - 'largestabs', 'lm': default, eigenvalues of largest magnitude
%   - 'largestreal', 'lr': eigenvalues with largest real part
%   - 'largestimag', 'li': eigenvalues with largest imaginary part.
%   - 'smallestabs', 'sm': default, eigenvalues of smallest magnitude
%   - 'smallestreal', 'sr': eigenvalues with smallest real part
%   - 'smallestimag', 'si': eigenvalues with smallest imaginary part.
%   - numeric : eigenvalues closest to sigma.
%
% Keyword Arguments
% -----------------
% Tol : :class:`double`
%   tolerance of the algorithm.
%
% Algorithm : :class:`char`
%   choice of eigensolver algorithm. Currently there is a choice between the use
%   of Matlab's buitin `eigs` specified by the identifiers 'eigs' or
%   'KrylovSchur', or the use of a custom Arnolid algorithm specified by
%   the identifier 'Arnoldi'.
%
% MaxIter : :class:`int`
%   maximum number of iterations, 100 by default.
%
% KrylovDim : :class:`int`
%   number of vectors kept in the Krylov subspace.
%
% IsSymmetric : :class:`logical`
%   flag to speed up the algorithm if the operator is symmetric, false by
%   default.
%
% Verbosity : :class:`.Verbosity`
%   Level of output information, by default nothing is printed if `flag` is
%   returned, otherwise only warnings are given, defaults to :code:`Verbosity.warn`.
%
% Returns
% -------
% V : (1, howmany) :class:`vector`
%   vector of eigenvectors.
%
% D : :class:`numeric`
%   vector of eigenvalues if only a single output argument is asked, diagonal
%   matrix of eigenvalues otherwise.
%
% flag : :class:`int`
%   if flag = 0 then all eigenvalues are converged, otherwise not.

arguments
    A
    v
    howmany = 1
    sigma = 'lm'

    options.Algorithm {mustBeMember(options.Algorithm, ...
        {'eigs', 'KrylovSchur', 'Arnoldi'})} = 'Arnoldi'

    options.Tol = eps(underlyingType(v))^(3/4)
    options.MaxIter = 100
    options.KrylovDim = 20
    options.DeflateDim
    options.ReOrth = 2
    options.NoBuild
    options.Verbosity = Verbosity.warn
    options.IsSymmetric logical = false
end

switch options.Algorithm

    case {'Arnoldi'}
        alg_opts = rmfield(options, {'Algorithm', 'IsSymmetric'});
        kwargs = namedargs2cell(alg_opts);
        alg = Arnoldi(kwargs{:});
        [varargout{1:nargout}] = eigsolve(alg, A, v, howmany, sigma);

    case {'eigs', 'KrylovSchur'}
        alg_opts = rmfield(options, ...
            {'Algorithm', 'DeflateDim', 'ReOrth', 'NoBuild', 'IsSymmetric'});
        kwargs = namedargs2cell(alg_opts);
        alg = KrylovSchur(kwargs{:});
        [varargout{1:nargout}] = eigsolve(alg, A, v, howmany, sigma, ...
            'IsSymmetric', options.IsSymmetric);

end

end
