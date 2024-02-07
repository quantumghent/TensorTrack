function m = mod1(x, y)
% Modulus after division, counting from :code:`1:y`.
%
% Arguments
% ---------
% x : :class:`int`
%   numerator.
%
% y : :class:`int`
%   divisor.
%
% Returns
% -------
% m : :class:`int`
%   remainder after division, where a value of 0 is replaced with y.

m = mod(x-1, y) + 1;

end

