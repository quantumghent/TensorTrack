function C = tensorprod(A, B, dimA, dimB, cA, cB, options)
% Tensor products between two tensors.
%
% Usage
% -----
% :code:`C = tensorprod(A, B, dimA, dimB)`
% 
% returns the tensor product of tensors :code:`A` and :code:`B`. The arguments :code:`dimA`
% and :code:`dimB` are vectors that specify which dimensions to contract in :code:`A` and
% :code:`B`. The size of the
% output tensor is the size of the uncontracted dimensions of :code:`A` followed by the size
% of the uncontracted dimensions of :code:`B`.
%
% :code:`C = tensorprod(A, B)`
% 
% returns the outer product of tensors :code:`A` and :code:`B`. This is equivalent to the
% previous syntax with :code:`dimA = dimB = []`.
%
% :code:`C = tensorprod(_, NumDimensionsA=ndimsA)`
% 
% optionally specifies the number of dimensions in tensor :code:`A` in addition to combat the
% removal of trailing singleton dimensions.

arguments
    A
    B
    dimA = []
    dimB = []
    cA = false
    cB = false
    options.NumDimensionsA = ndims(A)
end

% persistent version
% if isempty(version), version = ~isMATLABReleaseOlderThan("R2022a", "release", 1); end
% if version
%     C = builtin('tensorprod', dimA, dimB, 'NumDimensionsA', options.NumDimensionsA);
%     return
% end

szA = size(A, 1:options.NumDimensionsA);

szB = size(B);
szB = [szB ones(1, max(0, max(dimB) - length(szB)))];

uncA = 1:length(szA); uncA(dimA) = [];
uncB = 1:length(szB); uncB(dimB) = [];

if isempty(uncA)
    if isempty(uncB)
        szC = [1 1];
    elseif length(uncB) == 1
        szC = [1 szB(uncB)];
    else
        szC = szB(uncB);
    end
elseif isempty(uncB)
    if length(uncA) == 1
        szC = [szA(uncA) 1];
    else
        szC = szA(uncA);
    end
else
    szC = [szA(uncA) szB(uncB)];
end

C = reshape( ...
    reshape(permute(A, [uncA dimA]), prod(szA(uncA)), prod(szA(dimA))) * ...
    reshape(permute(B, [dimB uncB]), prod(szB(dimB)), prod(szB(uncB))), ...
    szC);

end
