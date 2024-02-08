function bool = isapprox(A, B, tol)
% Verify whether two arrays are approximately equal, based on their Euclidean distance.
%
% Arguments
% ---------
% A, B : :class:`numeric`
%   input arrays of the same size.
%
% Keyword Arguments
% -----------------
% AbsTol : :class:`double`
%   absolute tolerance.
%
% RelTol : :code:`double`
%   relative tolerance.
%
% Returns
% -------
% bool : :code:`logical`
%   true if :code:`A` and :code:`B` are approximately equal.

arguments
    A
    B
    tol.RelTol = max(sqrt(eps(underlyingType(A))), ...
        sqrt(eps(underlyingType(B))))
    tol.AbsTol = 0
end

% assert(all(size(A) == size(B)), 'Incompatible sizes'); 
bool = distance(A, B) <= max(tol.AbsTol, tol.RelTol * ...
    max(norm(reshape(A, 1, [])), norm(reshape(B, 1, []))));
end
