function j = next(i, total)
%NEXT Give the next index in a cyclic loop.
%   j = next(i, total)
%       gives the i + 1' th index, but loops back to 1 when j > total.
%
%   See also prev

j = mod(i, total) + 1;

end

