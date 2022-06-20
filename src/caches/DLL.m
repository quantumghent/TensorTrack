classdef DLL < handle
    % An element of a doubly-linked list. This is a list where each element consists of a
    % value, a reference to the previous element and a reference to the next element.
    %
    % Based on the work by Richard Lange (2022). Least-Recently Used (LRU) Cache 
    % (https://www.mathworks.com/matlabcentral/fileexchange/68836-least-recently-used-lru-cache),
    % MATLAB Central File Exchange. Retrieved June 18, 2022. 
    
    properties
        val  % data stored in this element
    end
    
    properties (Access = private)
        next % reference to next DLL object
        prev % reference to previous DLL object
    end
    
    methods
        function obj = DLL(val)
            % Construct an element of a doubly-linked list. By default, the detached links
            % point to the object itself.
            %
            % Arguments
            % ---------
            % val : any
            %   data stored in this element.
            %
            % Returns
            % -------
            % dll : :class:`DLL`
            %   data wrapped in a doubly-linked list format.
            
            obj.val = val;
            obj.next = obj;
            obj.prev = obj;
        end
        
        function obj = pop(obj)
            % Remove an object from a doubly-linked list, updating the links on the next and
            % previous element.
            %
            % Arguments
            % ---------
            % obj : :class:`DLL`
            %   object to remove from the list.
            %
            % Returns
            % -------
            % obj : :class:`DLL`
            %   removed object, with detached links.
            
            obj.prev.next = obj.next;
            obj.next.prev = obj.prev;
            obj.next = obj;
            obj.prev = obj;
        end
        
        function obj = append(obj, other)
            % Append an object to a doubly-linked list, updating the links.
            %
            % Arguments
            % ---------
            % obj : :class:`DLL`
            %   list to append to.
            %
            % other : :class:`DLL`
            %   object to append.
            %
            % Returns
            % -------
            % obj : :class:`DLL`
            %   updated list.
            
            other.next = obj.next;
            other.prev = obj;
            obj.next.prev = other;
            obj.next = other;
        end
        
        function other = uplus(obj)
            % Get the next element in the list.
            %
            % Usage
            % -----
            % :code:`other = +obj;`
            %
            % :code:`other = uplus(obj);`
            %
            % Arguments
            % ---------
            % obj : :class`DLL`
            %   current element in the list.
            %
            % Returns
            % -------
            % other : :class`DLL`
            %   next element in the list.
            
            other = obj.next;
        end
        
        function other = uminus(obj)
            % Get the previous element in the list.
            %
            % Usage
            % -----
            % :code:`other = -obj;`
            %
            % :code:`other = uminus(obj);`
            %
            % Arguments
            % ---------
            % obj : :class`DLL`
            %   current element in the list.
            %
            % Returns
            % -------
            % other : :class`DLL`
            %   previous element in the list.
            
            other = obj.prev;
        end
    end
end