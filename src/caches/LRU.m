classdef LRU < handle
    % A least-recently-used cache. Stores data up to a preset memory limit, then removes
    % the least-recently-used elements to free up space for additional data.
    %
    % Properties
    % ----------
    % sentinel : :class:`.DLL`
    %   sentinel of DLL, +sentinel is MRU, -sentinel is LRU
    %
    % map : :class:`containers.Map`
    %   map of key --> dll
    %
    % itemlimit : :class:`int`
    %   maximum size of cache in number of items
    %
    % memlimit : :class:`double`
    %   maximum size of cache in bytes
    %
    % mem : :class:`double`
    %   current memory usage in bytes.
    
    
    properties (Access = private)
        sentinel
        map 
        itemlimit = Inf
        memlimit = 20 * 2^30
        mem = 0;
    end
    methods
        function cache = LRU(itemlimit, memlimit)
            % Construct a new LRU cache.
            %
            % Arguments
            % ---------
            % itemlimit : :class:`int`
            %   maximum size of cache in number of items.
            %
            % memlimit : :class:`double`
            %   maximum size of cache in number of bytes.
            %
            % Returns
            % -------
            % cache : :class:`.LRU`
            %   empty LRU cache.
            
            % Initialize data
            cache.sentinel = DLL([]);
            cache.map      = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            % Set config
            if nargin > 0 && ~isempty(itemlimit)
                cache.itemlimit = itemlimit; 
            end
            if nargin > 1 && ~isempty(memlimit)
                cache.memlimit = memlimit;
            end
        end
        
        function val = get(cache, key)
            % Retrieve a value from the cache.
            %
            % Arguments
            % ---------
            % cache : :class:`.LRU`
            %   data cache.
            %
            % key : :class:`.uint8`
            %   data key.
            %
            % Returns
            % -------
            % val : :class:`any`
            %   value that is stored with a key, or empty if key not in cache.
            
            if isKey(cache.map, key)
                dll = cache.map(key);
                val = dll.val{2};
                
                % Re-insert dll to the front of the list
                pop(dll);
                append(cache.sentinel, dll);
            else
                val = [];
            end
        end
        
        function cache = set(cache, key, val)
            % Add a value to the cache.
            %
            % Arguments
            % ---------
            % cache : :class:`.LRU`
            %   data cache.
            %
            % key : :class:`uint8`
            %   data key.
            %
            % val : :class:`any`
            %   data value.
            %
            % Returns
            % -------
            % cache : :class:`.LRU`
            %   updated cache.
            
            % remove previously stored data
            if isKey(cache.map, key)
                dll_old = cache.map(key);
                pop(dll_old);
                cache.mem = cache.mem - memsize(dll_old.val{2}, 'B');
                cache.map.remove(key);
            end
            
            % add new data
            dll_new = DLL({key, val});
            append(cache.sentinel, dll_new);
            cache.map(key) = dll_new;
            
            % update memory
            cache.mem = cache.mem + memsize(val, 'B');
            while cache.mem > cache.memlimit || length(cache.map) > cache.itemlimit
                dll_oldest = -cache.sentinel;
                cache.mem = cache.mem - memsize(dll_oldest.val{2}, 'B');
                cache.map.remove(dll_oldest.val{1});
                pop(dll_oldest);
            end
        end
        
        function disp(cache)
            fprintf('LRU cache:\n');
            fprintf('\tbytes = %.2f / %.2f %s\n', ...
                cache.mem / 1024^2, cache.memlimit / 1024^2, 'MB');
            fprintf('\titems = %d / %d\n', cache.map.Count, cache.itemlimit);
        end
    end
end

% function b = memsize(x)
% b = numel(getByteStreamFromArray(x));
% end