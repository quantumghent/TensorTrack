classdef LRU < handle
    % LRU a least-recently-used cache. Stores data up to a preset memory limit, then removes
    % the least-recently-used elements to free up space for additional data.
    
    
    properties (Access = private)
        sentinel            % sentinel of DLL, +sentinel is MRU, -sentinel is LRU
        map                 % map of key --> dll
        itemlimit = Inf     % maximum size of cache in number of items
        memlimit = 6 * 2^30 % maximum size of cache in bytes
        mem = 0;            % current memory usage in bytes.
    end
    
    methods
        function cache = LRU(itemlimit, memlimit)
            % Construct a new LRU cache.
            %
            % Arguments
            % ---------
            % itemlimit : int
            %   maximum size of cache in number of items.
            %
            % memlimit : numeric
            %   maximum size of cache in number of bytes.
            %
            % Returns
            % -------
            % cache : :class:`LRU`
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
            % cache : :class:`LRU`
            %   data cache.
            %
            % key : :class:`uint8`
            %   data key.
            %
            % Returns
            % -------
            % val : any
            %   value that is stored with a key, or empty if key not in cache.
            
            try
                dll = cache.map(key);
                val = dll.val{2};
                
                % Re-insert dll to the front of the list
                pop(dll);
                append(cache.sentinel, dll);
            catch
                val = [];
            end
        end
        
        function cache = set(cache, key, val)
            % Add a value to the cache.
            %
            % Arguments
            % ---------
            % cache : :class:`LRU`
            %   data cache.
            %
            % key : :class:`uint8`
            %   data key.
            %
            % val : any
            %   data value.
            %
            % Returns
            % -------
            % cache : :class:`LRU`
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
                pop(dll_oldest);
                cache.mem = cache.mem - memsize(dll_oldest.val{2}, 'B');
                cache.map.remove(key);
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