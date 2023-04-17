function j = prev(i, total)
% Give the previous index in a cyclic loop.
%
% Usage
% -----
% :code:`j = prev(i, total)`
% 	gives the :math:`i - 1`'th index, but loops back to total when :math:`j < 1`.
%
% See Also
% --------
% :func:`next`

j = mod(i - 2, total) + 1;

end
