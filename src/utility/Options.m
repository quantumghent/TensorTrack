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
        
        function bool = Debug(bool)
            persistent dodebug
            if nargin > 0
                dodebug = bool;
                return
            end
            
            if isempty(dodebug), dodebug = false; end
            bool = dodebug;
        end
    end
end

