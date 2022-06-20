classdef Options
    % Various package settings.

    
    methods (Static)
        function bool = CacheEnabled(bool)
            persistent isenabled
            if nargin > 0
                isenabled = bool;
                return
            end
            
            if isempty(isenabled), isenabled = true; end
            bool = isenabled;
        end
    end
end

