function [bytes, unit] = memsize(in, unit)
% Check memory size of object in prefered unit ('GB', 'MB', 'KB' or 'B').
%
% Usage
% -----
% :code:`[bytes, unit] = memsize(in, unit)`

if isa(in, 'containers.Map')
    warning('off', 'MATLAB:structOnObject');
    in = struct(in);
    warning('on', 'MATLAB:structOnObject');
end

if ~isobject(in)
    s = whos('in');
    bytes = s.bytes;
else
    props = properties(in);
    bytes = 0;
    for ii = 1:length(props)
        bytes = bytes + memsize({in.(props{ii})}, 'B');
    end
end

if nargin == 1 || isempty(unit)
    if bytes > 1024^3
%         bytes = bytes / 1024^3;
        unit = 'GB';
    elseif bytes > 1024^2
%         bytes = bytes / 1024^2;
        unit = 'MB';
    elseif bytes > 1024
%         bytes = bytes / 1024;
        unit = 'KB';
    else
        unit = 'B';
    end
end

if strcmpi(unit, 'B')
    return;
end

if strcmpi(unit, 'KB')
    bytes = bytes / 1024;
    return;
end

if strcmpi(unit, 'MB')
    bytes = bytes / 1024^2;
    return;
end

if strcmpi(unit, 'GB')
    bytes = bytes / 1024^3;
    return
end

warning('Unrecognised unit for bytes, using B instead.');
unit = 'B';

end
