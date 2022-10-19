function bool = isapprox(A, B, tol)

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
