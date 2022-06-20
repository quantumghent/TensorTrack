function inds = unique1(inds)
% unique1 - Find the elements that appear exactly once.
%   inds = unique1(inds)
%       deletes all elements that appear more than once.

inds = inds(sum(inds(:) == inds) == 1);

end
