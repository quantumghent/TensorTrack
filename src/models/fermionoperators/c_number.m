function n = c_number()
% Fermionic number operator.
%
% Returns
% -------
% n : :class:`.Tensor`
%   number operator represented as a 2-leg tensor with :math:`fZ_2` symmetry.

pspace = fZ2Space([0 1], [1 1], false);

n = fill_matrix(Tensor.zeros(pspace, pspace), {0, 1});

end
