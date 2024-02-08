function R = randc(varargin)
% Generate an array containing pseudorandom complex values, both real and imaginary part 
% drawn from the standard uniform distribution on the open interval (0, 1).
%
% Usage
% -----
% :code:`R = randc([m n])`
%
% :code:`R = randc(m, n, p, ...)`
%
% :code:`R = randc(..., classname)`
%
% :code:`R = randc(..., 'like', Y)`
%
% Arguments
% ---------
% m, n, p, ... : :class:`int`
%   integers defining the size of the output array.
%
% classname : :class:`char`
%   datatype of the array, default :code:`'double'`.
%
% Y : :class:`numeric`
%   create an array of the same class as Y.
%
% Returns
% -------
% R : :class:`numeric`
%   complex pseudorandom values, real and imaginary part distributed from the uniform
%   distribution.

R = complex(rand(varargin{:}), rand(varargin{:}));

end

