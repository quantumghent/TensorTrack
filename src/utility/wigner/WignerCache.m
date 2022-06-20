function [key, sign] = WignerCache(j1, j2, j3, m1, m2, m3)

R = [-j1 + j2 + j3, j1 - j2 + j3, j1 + j2 - j3
    j1 - m1, j2 - m2, j3 - m3
    j1 + m1, j2 + m2, j3 + m3];
sign = 1;

% check if magic square
rowsums = sum(R, 1);
colsums = sum(R, 2);
if any(rowsums(1) ~= rowsums) || any(rowsums(1) ~= colsums)
    key = 1;
    return
end

sign = 1;

% bring S to the first row and column:
[S, i] = min(R, [], 'all', 'linear');
[row, col] = ind2sub([3 3], i);
if row ~= 1
    R([1 row], :) = R([row 1], :);
    sign = -sign;
end
if col ~= 1
    R(:, [1 col]) = R(:, [col, 1]);
    sign = -sign;
end

% bring L to the first row and second column:
[L1, i] = max(R(1,:), [], 'all', 'linear');
[L2, j] = max(R(:,1), [], 'all', 'linear');

if L2 > L1
    % make row 2 and then transpose
    if j ~= 2
        R([3 2], :) = R([2 3], :);
        sign = -sign;
    end
    R = R.';
    L = L2;
else
    if i ~= 2
        R(:, [3 2]) = R(:, [2 3]);
        sign = -sign;
    end
    L = L1;
end

if R(2,2) > R(3,2) || (R(2,2) == R(3,2) && R(2,3) > R(3,3))
    R([2 3], :) = R([3 2], :);
    sign = -sign;
end

X = R(2,1);
B = R(2,2);
T = R(3,3);

key = 1 / 120 * (L * (24 + L * (50 + L * (35 + L * (10 + L))))) + ...
    1 / 24 * (X * (6 + X * (11 + X * (6 + X)))) + ...
    1 / 6 * (T * (2 + T * (3 + T))) + ...
    1 / 2 * (B * (B + 1)) + ...
    S + 2;

end

