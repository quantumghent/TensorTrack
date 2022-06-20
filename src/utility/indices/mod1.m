function m = mod1(x, y)
% Modulus after division, counting from 1:y.
%
% Arguments
% ---------
% x : int
%   numerator.
%
% y : int
%   divisor.
%
% Returns
% -------
% m : int
%   remainder after division, where a value of 0 is replaced with y.

m = mod(x-1, y) + 1;

end

