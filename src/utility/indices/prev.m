function j = prev(i, total)
%PREV Give the previous index in a cyclic loop.
%   j = prev(i, total)
%       gives the i - 1'th index, but loops back to total when j < 1.
%
%   See also next

j = mod(i - 2, total) + 1;

end
