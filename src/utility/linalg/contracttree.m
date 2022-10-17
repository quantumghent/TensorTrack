function [C, ic, cc] = contracttree(tensors, indices, conjlist, tree, debug)

if isnumeric(tree)
    C = tensors{tree};
    ic = indices{tree};
    cc = conjlist(tree);
    return
end

[A, ia, ca] = contracttree(tensors, indices, conjlist, tree{1}, debug);
[B, ib, cb] = contracttree(tensors, indices, conjlist, tree{2}, debug);
[dimA, dimB] = contractinds(ia, ib);

if debug
    contractcheck(A, ia, ca, B, ib, cb);
end

C = tensorprod(A, B, dimA, dimB, ca, cb, 'NumDimensionsA', length(ia));

ia(dimA) = [];
ib(dimB) = [];
ic = [ia ib];
cc = false;

end
