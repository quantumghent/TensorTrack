function s = dim2str(sz)

s = regexprep(mat2str(sz), {'\[', '\]', '\s+'}, {'', '', 'x'});

end

