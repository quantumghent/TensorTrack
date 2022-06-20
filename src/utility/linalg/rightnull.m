function N = rightnull(A, alg, atol)
% Compute the right nullspace of a matrix, such that A * N' = 0.
%
% Arguments
% ---------
% A : numeric
%   input matrix.
%
% alg : 'lq' or 'svd'
%   choice of algorithm, default 'svd'.
%
% atol : double
%   absolute tolerance for the null space, defaults to max(size(A)) * eps(max(svd(A))).
%
% Returns
% -------
% N : numeric
%   orthonormal basis for the orthogonal complement of the support of the column space of A.

arguments
    A
    alg = 'svd'
    atol = []
end

[m, n] = size(A);

switch alg
    case 'svd'
        [~, S, V] = svd(A, 0);
        
        if m == 1
            s = S(1);
        else
            s = diag(S);
        end
        
        if isempty(atol)
            atol = max(m, n) * eps(max(s));
        end
        
        r = sum(s > atol);
        N = V(:, r+1:n)';
        
    case 'lq'
        [Q, R] = qr(A.');
        
        r = size(R, 2);
        N = Q(:, r+1:n).';

    otherwise
        error('Invalid algorithm');
end

end
