function inds = unique1(inds)
% Find the elements that appear exactly once.
%
% Usage
% -----
% :code:`inds = unique1(inds)`
% deletes all elements that appear more than once.

inds = inds(sum(inds(:) == inds) == 1);

end
