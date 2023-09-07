function clear_path(filename)
    

if isfile(filename)
    [filepath, name, ext] = fileparts(filename);
    
    pat = "_backup" + digitsPattern;
    if endsWith(name, pat)
        number = extractAfter(name, "_backup");
        newname = replace(name, pat, "_backup" + (double(number) + 1));
    else
        newname = name + "_backup1";
    end
    newfilename = fullfile(filepath, newname + ext);
    clear_path(newfilename);
    movefile(filename, newfilename);
end


end
