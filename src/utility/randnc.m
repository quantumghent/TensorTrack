function R = randnc(varargin)
% Generate an array containing pseudorandom complex values, both real and imaginary part 
% drawn from the standard normal distribution.
%
% Usage
% -----
% :code:`R = randnc([m n])`
%
% :code:`R = randnc(m, n, p, ...)`
%
% :code:`R = randnc(..., classname)`
%
% :code:`R = randnc(..., 'like', Y)`
%
% Arguments
% ---------
% m, n, p, ... : int
%   integers defining the size of the output array.
%
% classname : :class:`char`
%   datatype of the array, default :class:`double`.
%
% Y : numeric
%   create an array of the same class as Y.
%
% Returns
% -------
% R : numeric
%   complex pseudorandom values, real and imaginary part distributed from the normal
%   distribution.

R = complex(randn(varargin{:}), randn(varargin{:}));

end
