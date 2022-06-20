function n = treeindsequence(n)
% Compute the number of edge labels needed in a splitting or fusing tree.
%
% Usage
% -----
% ```t = treeindsequence(n)```
% 
% Arguments
% ---------
% n : int
%   number of external edges
%
% Returns
% -------
% t : int
%   total number of edges
%
% Examples
% --------
% ```treeindsequence(0:4)```
%   ans = [0 1 2 4 6]

n = max(n, 2 * n - 2);

end
