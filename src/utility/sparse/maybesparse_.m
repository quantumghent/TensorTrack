function a = maybesparse_(a)
if issparse(a)
    a = SparseArray(a);
end
end

