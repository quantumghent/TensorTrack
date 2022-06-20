function s = perm2swap(p)
% Convert a permutation into a sequence of swaps on neighbouring indices.
%
% Arguments
% ---------
% p : int
%   permutation vector.
%
% Returns
% -------
% s : int
%   list of swaps that compose into the permutation vector, where `i` indicates a swap
%   between indices `i` and `i+1`.

N = length(p);
s = [];
for k = 1:N - 1
    s = [s flip(k:p(k) - 1)]; %#ok<AGROW>
    p2 = p(k + 1:N);
    p2(p2 < p(k)) = p2(p2 < p(k)) + 1;
    p(k + 1:N) = p2;
    p(k) = k;
end

end
