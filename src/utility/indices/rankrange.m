function r = rankrange(rank)
% Convert tensor rank into a contiguous range of dimensions.

r = [1:rank(1) rank(1) + (rank(2):-1:1)];

end
