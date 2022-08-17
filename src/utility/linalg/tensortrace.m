function C = tensortrace(A, ia)
% tensortrace - Compute the (partial) trace of a tensor.
%   [C, ic] = tensortrace(A, ia)
%       traces over the indices that appear twice in ia.
%
%   [C, ic] = tensortrace(A, ia, ic)
%       optionally specifies the output indices' order.

arguments
    A
    ia
end

[dimA1, dimA2] = traceinds(ia);

szA1 = size(A, dimA1);
szA2 = size(A, dimA2);

assert(all(szA1 == szA2));

E = reshape(eye(prod(szA1)), [szA1 szA2]);

indsA = -(1:ndims(A));
indsA(dimA1) = 1:length(dimA1);
indsA(dimA2) = (1:length(dimA2)) + length(dimA1);
C = contract(A, indsA, E, [1:length(dimA1) (1:length(dimA2)) + length(dimA1)]);

end
