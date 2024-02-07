function N = leftnull(A, alg, atol)
% Compute the left nullspace of a matrix, such that :code:`N' * A == 0`.
%
% Arguments
% ---------
% A : :class:`numeric`
%   input matrix.
%
% alg : :class:`char`, 'qr' or 'svd'
%   choice of algorithm, default 'svd'.
%
% atol : :class:`double`
%   absolute tolerance for the null space, defaults to
%   :code:`max(size(A)) * eps(max(svd(A)))`.
%
% Returns
% -------
% N : :class:`numeric`
%   orthonormal basis for the orthogonal complement of the support of the row space of A.

arguments
    A
    alg = 'svd'
    atol = []
end

[m, n] = size(A);

switch alg
    case 'svd'
        [U, S] = svd(A);

        if n == 1
            s = S(1);
        else
            s = diag(S);
        end

        if isempty(atol)
            atol = max(m, n) * eps(max(s));
        end

        r = sum(s > atol);
        N = U(:, r+1:m);
        
    case 'qr'
        [Q, R] = qr(A);

        r = size(R, 2);
        N = Q(:, r+1:m);
        
    otherwise
        error('Invalid algorithm');
end

end
