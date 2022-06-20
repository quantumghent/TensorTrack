function inds = unique2(inds)
% unique2 - Find the elements that appear more than once.
%   inds = unique2(inds)
%       deletes all elements that appear only once.

[~, ia1] = unique(inds, 'first');
[~, ia2] = unique(inds, 'last');
inds(ia1 == ia2) = [];

end
