function s = dim2str(sz)
% Convert size array to string output.

s = regexprep(mat2str(sz), {'\[', '\]', '\s+'}, {'', '', 'x'});

end

