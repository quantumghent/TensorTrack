function bool = isposdef(A)
% Verify whether a matrix is positive definite.

[~, flag] = chol(A);
bool = flag == 0;

end
