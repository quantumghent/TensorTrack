function j = next(i, total)
% Give the next index in a cyclic loop.
%
% Usage
% -----
% :code:`j = next(i, total)`
%	gives the :math:`i + 1`'th index, but loops back to 1 when :math:`j > \text{total}`.
%
% See Also
% --------
% :func:`prev`

j = mod(i, total) + 1;

end

