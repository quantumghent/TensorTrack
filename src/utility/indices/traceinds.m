function [dimA1, dimA2] = traceinds(ia)

ind = find(ia(:) == ia & ~tril(true(length(ia)))).' - 1;
dimA1 = mod(ind, length(ia)) + 1;
dimA2 = floor(ind / length(ia)) + 1;

end
