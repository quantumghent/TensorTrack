function [dimA, dimB] = contractinds(ia, ib)
% contractinds - Find the contracted dimensions.
%   [dimA, dimB] = contractinds(ia, ib)
%       locates the repeated indices.

ind = find(ia(:) == ib).' - 1;
dimA = mod(ind, length(ia)) + 1;
dimB = floor(ind / length(ia)) + 1;

% slower than the automatic:
% [dimA, dimB] = ind2sub([length(ia) length(ib)], find(ia(:) == ib));

end
