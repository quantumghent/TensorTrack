function [dimA1, dimA2] = traceinds(ia)
% Find the traced dimensions.
%
% :code:`[dimA1, dimA2] = traceinds(ia)`
% locates the repeated indices in a vectors of integers.

ind = find(ia(:) == ia & ~tril(true(length(ia)))).' - 1;
dimA1 = mod(ind, length(ia)) + 1;
dimA2 = floor(ind / length(ia)) + 1;

end
