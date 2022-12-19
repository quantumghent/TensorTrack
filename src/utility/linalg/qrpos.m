function [Q, R] = qrpos(A, varargin)
% Positive orthogonal-triangular decomposition.

[Q, R] = qr(A, varargin{:});
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

end
