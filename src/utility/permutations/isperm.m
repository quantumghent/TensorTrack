function bool = isperm(p)
% Check if a vector is a permutation.
%
% Usage
% -----
% :code:`bool = isperm(p)`

bool = all(sort(p) == 1:length(p));

end

