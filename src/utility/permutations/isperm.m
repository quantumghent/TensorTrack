function bool = isperm(p)
% isperm - Check if a vector is a permutation.
%   bool = isperm(p)

bool = all(sort(p) == 1:length(p));

end

