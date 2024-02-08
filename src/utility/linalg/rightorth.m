function [R, Q] = rightorth(A, alg)
% Factorize a matrix into an orthonormal basis `Q` and remainder `L`, such that
% :code:`A = L * Q`.
%
% Usage
% -----
% :code:`[R, Q] = rightorth(A, alg)`
%
% Arguments
% ---------
% A : :class:`numeric`
%   input matrix to factorize.
%
% alg : :class:`char` or :class:`string`
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
% R : :class:`numeric`
%   remainder matrix, depends on selected algorithm.
%
% Q : :class:`numeric`
%   orthonormal basis matrix.

% TODO have a look at https://github.com/iwoodsawyer/factor

switch alg
    case 'svd'
        [U, S, Q] = svd(A, 0);
        R = U * S;
        Q = Q';
        return;
        
    case 'polar'
        newalg = 'polar';
    case 'lq'
        newalg = 'qr';
    case 'lqpos'
        newalg = 'qrpos';
    case 'rq'
        newalg = 'ql';
    case 'rqpos'
        newalg = 'qlpos';
end

[Q, R] = leftorth(A.', newalg);
Q = Q.';
R = R.';

end
