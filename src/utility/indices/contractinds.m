function [dimA, dimB] = contractinds(ia, ib)
% Find the contracted dimensions.
%
% :code:`[dimA, dimB] = contractinds(ia, ib)`
% locates the repeated indices in two vectors of integers.

ind = find(ia(:) == ib).' - 1;
dimA = mod(ind, length(ia)) + 1;
dimB = floor(ind / length(ia)) + 1;

% slower than the automatic:
% [dimA, dimB] = ind2sub([length(ia) length(ib)], find(ia(:) == ib));

end
