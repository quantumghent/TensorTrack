function obj = mat2jl(obj)
% Recursively convert an object into something that can be read by Julia (GlueKit)

if strcmp(class(obj), "Tensor")
    % Reshuffle matrixblocks to compensate for differences between MatLab
    % and Julia order inside blocks.
    obj = juliasort(obj);
end

if isobject(obj) && ~strcmp(class(obj), "string")
    name = class(obj);
    orig = warning('off', 'MATLAB:structOnObject');
    obj = arrayfun(@struct, obj);
    warning(orig);
    for i = 1:numel(obj)
        obj(i).classname = name;
    end
end

if isstruct(obj)
    fns = fieldnames(obj);
    for i = 1:numel(obj)
        for fn = fns'
            obj(i).(fn{1}) = mat2jl(obj(i).(fn{1}));
        end
    end
elseif iscell(obj)
    obj = cellfun(@mat2jl, obj, 'UniformOutput', false);
end

end

