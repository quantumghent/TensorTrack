function C = tensortrace(A, i1, i2)
% Compute the (partial) trace of a tensor.
%
% Usage
% -----
%
% :code:`[C, ic] = tensortrace(A, ia)`
% 	traces over the indices that appear twice in ia.
%
% :code:`[C, ic] = tensortrace(A, ia, ic)`
% 	optionally specifies the output indices' order.

if isempty(i1) && isempty(i2), C = A; return; end
assert(length(i1) == length(i2), 'invalid indices');

szA1 = size(A, i1);
szA2 = size(A, i2);

assert(all(szA1 == szA2));

E = reshape(eye(prod(szA1)), [szA1 szA2]);

indsA = -(1:ndims(A));
indsA(i1) = 1:length(i1);
indsA(i2) = (1:length(i2)) + length(i1);
C = contract(A, indsA, E, [1:length(i1) (1:length(i2)) + length(i1)]);

end
