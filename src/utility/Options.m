classdef Options
    % Various package settings.
    
    methods (Static)
        function bool = CacheEnabled(bool)
            % Enable cache.
            %
            % Usage
            % -----
            % :code:`bool = Options.CacheEnabled(bool)` sets the cache enabling to :code:`bool`.
            %
            % :code:`bool = Options.CacheEnabled(bool)` returns the current cache enabling.
            persistent isenabled
            if nargin > 0
                isenabled = bool;
                return
            end
            
            if isempty(isenabled), isenabled = true; end
            bool = isenabled;
        end
        
        function bool = Debug(bool)
            % Enable cache.
            %
            % Usage
            % -----
            % :code:`bool = Options.Debug(bool)` sets the debug mode to :code:`bool`.
            %
            % :code:`bool = Options.Debug(bool)` returns the current debug mode.
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

