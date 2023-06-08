function safesave(filename, variable)
% Safe wrapper around save.
% Clears path and makes backup of existing files, and creates directories if not existing.

[filepath, name, ext] = fileparts(filename);
if ~isdir(filepath)
    mkdir(filepath);
else
    clear_path(filename);
end

save(filename, '-struct', 'variable');

end
