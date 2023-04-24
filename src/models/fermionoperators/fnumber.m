function n = fnumber()

pspace = fZ2Space([0 1], [1 1], false);

n = fill_matrix(Tensor.zeros(pspace, pspace), {0, 1});

end
