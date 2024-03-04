function [Q, R] = leftorth(A, alg)
% Factorize a matrix into an orthonormal basis `Q` and remainder `R`, such that
% :code:`A = Q * R`.
%
% Usage
% -----
% :code:`[Q, R] = leftorth(A, alg)`
%
% Arguments
% ---------
% A : :class:`numeric`
%   input matrix to factorize.
%
% alg : :class:`char` or :class:`string`
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
% Q : :class:`numeric`
%   orthonormal basis matrix
%
% R : :class:`numeric`
%   remainder matrix, depends on selected algorithm.

% TODO have a look at https://github.com/iwoodsawyer/factor

switch alg
    case 'qr'
        [Q, R] = qr(A, 0);
        
    case 'qrpos'
        [Q, R] = qrpos(A, 0);
        
    case 'ql'
        [Q, R] = qr(flip(A, 2), 0);
        Q = flip(Q, 2);
        R = flip(flip(R, 1), 2);
        
    case 'qlpos'
        [Q, R] = qr(flip(A, 2), 0);
        if isrow(Q)
            Q = Q * sign(R(1));
            R = R * sign(R(1));
        else
            D = diag(R);
            D(abs(D) < 1e-12) = 1;
            D = sign(D);
            Q = Q .* D';
            R = D .* R;
        end
        Q = flip(Q, 2);
        R = flip(flip(R, 1), 2);
        
    case 'polar'
        [U, S, V] = svd(A, 0);
        [m, n] = size(A);
        if m < n                % m > n handled by economy svd
            S = S(:, 1:m);
            V = V(:, 1:m);
        end
        Q = U * V';
        R = V * S * V';
        R = (R + R') / 2;   % force hermitian
        
    case 'svd'
        % TODO add tol?
        [Q, S, V] = svd(A, 0);
        R = S * V';
        
    otherwise
        error('Invalid alg (%s)', alg);
        
end

end
