function contractcheck(A, ia, ca, B, ib, cb)

Aspaces = space(A);
if ca, Aspaces = conj(Aspaces); end
Bspaces = space(B);
if cb, Bspaces = conj(Bspaces); end

[dimA, dimB] = contractinds(ia, ib);

for i = 1:length(dimA)
    assert(Aspaces(dimA(i)) == conj(Bspaces(dimB(i))), 'tensors:SpaceMismatch', ...
        'Invalid index %d:\n\t%s\n\tis incompatible with\n\t%s', ...
        ia(dimA(i)), string(Aspaces(dimA(i))), string(Bspaces(dimB(i))));
end

end