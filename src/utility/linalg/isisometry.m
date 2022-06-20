function bool = isisometry(A, side, tol)
% Verify whether a matrix is an isometry.
%
% Arguments
% ---------
% A : numeric
%   input matrix.
%
% side : 'left', 'right' or 'both'
%   check if A' * A == I, A * A' == I, or both by default.
%
% Keyword Arguments
% -----------------
% AbsTol : double
%   absolute tolerance
%
% RelTol : double
%   relative tolerance
%
% Returns
% -------
% bool : logical
%   true if A is an isometry.

arguments
    A
    side {mustBeMember(side, {'left', 'right', 'both'})} = 'both'
    tol.RelTol = sqrt(eps(underlyingType(A)))
    tol.AbsTol = 0
end

[m, n] = size(A);

switch side
    case 'left'
        if m < n, bool = false; return; end
        
        bool = isapprox(A' * A, eye(n), ...
            'RelTol', tol.RelTol, 'AbsTol', tol.AbsTol);
        
    case 'right'
        if n < m, bool = false; return; end
        
        bool = isapprox(A * A', eye(m), 'RelTol', tol.RelTol, 'AbsTol', tol.AbsTol);
        
    case 'both'
        if m ~= n, bool = false; return; end
        I = eye(m);
        bool = isapprox(A * A', I, 'RelTol', tol.RelTol, 'AbsTol', tol.AbsTol) && ...
            isapprox(A' * A, I, 'RelTol', tol.RelTol, 'AbsTol', tol.AbsTol);
end

end
