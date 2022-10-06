function I = ind2sub_(sz, ind, perm)
% Faster implementation of builtin ind2sub

if isempty(ind)
    I = [];
    return;
elseif size(sz, 2) == 1
    I = ind;
    return;
end

if nargin < 3
    perm = 1:numel(sz);
else
    perm(perm) = 1:numel(sz);
end

nout = numel(sz);
I = zeros(numel(ind), numel(sz));

if nout > 2
    k = cumprod(sz);
    for i = nout:-1:3
        I(:, perm(i)) = floor((ind-1) / k(i-1)) + 1;
        ind = rem(ind-1, k(i-1)) + 1;
    end
end

if nout >= 2
    I(:, perm(2)) = floor((ind-1)/sz(1)) + 1;
    I(:, perm(1)) = rem(ind-1, sz(1)) + 1;
else 
    I(:, perm(1)) = ind;
end

end