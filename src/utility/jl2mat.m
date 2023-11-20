function obj = jl2mat(obj)
% Recursively convert a Julia object back into something read by TensorTrack.

if isstruct(obj)
    fns = fieldnames(obj);
    for i = 1:numel(obj)
        for fn = fns'
            obj(i).(fn{1}) = jl2mat(obj(i).(fn{1}));
        end
    end
    if isfield(obj, 'classname')
        obj = str2obj(obj);
    end
elseif iscell(obj)
    obj = cellfun(@jl2mat, obj, 'UniformOutput', false);
end

end

function o = str2obj(s)

if numel(s) > 1
    for i = numel(s):-1:1
        o(i) = str2obj(s(i));
    end
    o = reshape(o, size(s));
    return
end

if ismember('Data', fieldnames(s))
    o = feval(s.classname, s.Data);
else
    o = feval(s.classname);
end

for name = fieldnames(s)'
    if strcmp(name{1}, 'classname') || strcmp(name{1}, 'Data')
        continue;
    end
    o.(name{1}) = s.(name{1});
end

end