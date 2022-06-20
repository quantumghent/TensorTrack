function bool = isapprox(A, B, tol)

arguments
    A
    B
    tol.RelTol = max(sqrt(eps(underlyingType(A))), ...
        sqrt(eps(underlyingType(B))))
    tol.AbsTol = 0
end

assert(all(size(A) == size(B)), 'Incompatible sizes'); 
bool = distance(A, B) <= max(tol.AbsTol, tol.RelTol * max(norm(A(:)), norm(B(:))));
end
