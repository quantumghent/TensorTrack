function [Q, R] = leftorth(A, alg)



% TODO have a look at https://github.com/iwoodsawyer/factor


switch alg
    case 'qr'
        [Q, R] = qr(A, 0);
        
    case 'qrpos'
        [Q, R] = qr(A, 0);
        D = diag(R);
        D(abs(D) < 1e-12) = 1;
        D = sign(D);
        Q = Q .* D';
        R = D .* R;
        
    case 'ql'
        [Q, R] = qr(flip(A, 2), 0);
        Q = flip(Q, 2);
        R = flip(flip(R, 1), 2);
        
    case 'qlpos'
        [Q, R] = qr(flip(A, 2), 0);
        D = diag(R);
        D(abs(D) < 1e-12) = 1;
        D = sign(D);
        Q = Q .* D';
        R = D .* R;
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
