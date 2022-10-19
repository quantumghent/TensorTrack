function timeStr = time2str(time, unit, precision)
%TIMESTRING Returns string version of a duration in seconds.

if nargin < 3, precision = 1; end


if nargin == 1 || isempty(unit)
    if time >= 3600
        unit = 'h';
    elseif time >= 60
        unit = 'm';
    else
        unit = 's';
    end
end


switch unit
    case 'h'
        time = time / 3600;
        
    case 'm'
        time = time / 60;
        
end


timeStr = sprintf('%0.*f%c', precision, time, unit);


end


