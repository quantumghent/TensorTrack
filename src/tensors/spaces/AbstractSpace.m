classdef (Abstract) AbstractSpace
    % Abstract structure of a tensor index.
    
    properties %(Access = protected)
        dimensions                      % Specification of the internal dimensions
        dual (1,1) logical = false    % Flag to indicate if the space is a dual space
    end
    
    %% Constructors
    methods
        function spaces = AbstractSpace(dimensions, dual)
            % Construct an array of vector spaces.
            %
            % Repeating Arguments
            % -------------------
            % dimensions : int or struct
            %   a variable which represents the internal dimension of the space.
            %
            % isdual : logical
            %   a variable which indicates if a space is dual.
            %
            % Returns
            % -------
            % spaces : :class:`AbstractSpace`
            %   array of spaces.
            
            arguments (Repeating)
                dimensions
                dual
            end
            
            if nargin == 0, return; end
            
            for i = length(dimensions):-1:1
                spaces(i).dimensions = dimensions{i};
                spaces(i).dual = dual{i};
            end
        end
    end
    
    methods (Static)
        function spaces = new(dimensions, dual)
            % Construct a vector space. This is a utility method to be able to access the
            % constructor of a subclass.
            %
            % Repeating Arguments
            % -------------------
            % dimensions : int or struct
            %   a variable which represents the internal dimension of the space.
            %
            % dual : logical
            %   a variable which indicates if a space is dual.
            %
            % Returns
            % -------
            % spaces : :class:`AbstractSpace`
            %   array of spaces.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
    end
    
    
    %% Structure
    methods
        function f = isdual(space)
            f = [space.dual];
        end
        
        function c = charges(spaces)
            % Compute all charge combinations of a space. If no internal structure is
            % present, this yields an empty result.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % c : (:, :) :class:`AbstractCharge`
            %   list of charge combinations, where each row is a combination.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
        
        function d = degeneracies(spaces)
            % Compute all degeneracy combinations of a product space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % d : (:, :) int
            %   list of degeneracy combinations, where each row is an entry.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
        
        function d = dims(spaces)
            % Compute the dimension of the spaces.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % d : (1, :) numeric
            %   total dimension of each of the input spaces.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
        
        function s = subspaces(s)
        end
        
        function trees = fusiontrees(codomain, domain)
            % Compute all allowed fusiontrees that connect domain and codomain. If the
            % spaces have no internal structure than this returns `[]`.
            %
            % Arguments
            % ---------
            % codomain, domain : :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % trees : :class:`FusionTree`
            %   list of fusiontrees that are allowed.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
        
        function style = braidingstyle(codomain, domain)
            % Determine the braiding style of the internal structure of a space.
            %
            % Arguments
            % ---------
            % codomain, domain : (1, :) :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % style : :class:`BraidingStyle`
            %   braiding style of the internal structure.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
        
        function style = fusionstyle(codomain, domain)
            % Determine the fusion style of the internal structure of a space.
            %
            % Arguments
            % ---------
            % codomain, domain : (1, :) :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % style : :class:`FusionStyle`
            %   fusion style of the internal structure.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
    end
    
    
    %% Manipulations
    methods
        function spaces = conj(spaces)
            % Compute the element-wise dual space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % spaces : (1, :) :class:`AbstractSpace`
            %   dual spaces.
            
            for i = 1:length(spaces)
                spaces(i).dual = ~spaces(i).dual;
            end
        end
        
        function spaces = ctranspose(spaces)
            % Compute the dual space of the product space, found by taking the reverse order
            % of the dual spaces.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`AbstractSpace`
            %   input product space.
            %
            % Returns
            % -------
            % spaces : (1, :) :class:`AbstractSpace`
            %   dual product space.
            
            spaces = conj(spaces(length(spaces):-1:1));
        end
        
        function space = mtimes(space1, space2)
            % Fuse two spaces to a single space.
            %
            % Arguments
            % ---------
            % space1, space2 : (1,1) :class:`AbstractSpace`
            %   input spaces.
            % 
            % Returns
            % -------
            % space : (1,1) :class:`AbstractSpace`
            %   fused space.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
        
        function space = plus(space1, space2)
            % Find the direct sum of two spaces.
            %
            % Arguments
            % ---------
            % space1, space2 : (1,1) :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % space : (1,1) :class:`AbstractSpace`
            %   direct sum space.
            
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end 
        
        function space = prod(spaces)
            % Fuse a product space to a single space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`AbstractSpace`
            %   Array of input spaces.
            %
            % Returns
            % -------
            % space : (1, 1) :class:`AbstractSpace`
            %   Fused space which is isomorphic to the input product space.
            
            % TODO this is probably faster by fusing from both ends to the middle.
            space = spaces(1);
            for i = 2:length(spaces)
                space = space * spaces(i);
            end
        end
        
        function space = sum(spaces)
            % Fuse a list of spaces to the direct sum space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`AbstractSpace`
            %   Array of input spaces.
            %
            % Returns
            % -------
            % space : (1, 1) :class:`AbstractSpace`
            %   direct sum space.
            space = spaces(1);
            for i = 2:length(spaces)
                space = space + spaces(i);
            end
        end
        
        function spaces = insertone(spaces, i, dual)
            % insert a trivial space at position i.
            arguments
                spaces
                i = length(spaces) + 1
                dual = false
            end
            
            trivialspace = one(spaces);
            if dual, trivialspace = conj(trivialspace); end
            spaces = [spaces(1:i-1) trivialspace spaces(i:end)];
        end
        
        
        %% Utility
        function bools = eq(spaces1, spaces2)
            % Verify if spaces are element-wise equal.
            %
            % Arguments
            % ---------
            % spaces1, spaces2 : (1, :) :class:`AbstractSpace`
            %   input spaces, either of equal size or scalar.
            %
            % Returns
            % -------
            % bools : (1, :) logical
            %   flags that indicate if the element spaces are equal.
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
        
        function bool = lt(space1, space2)
            bool = ~ge(space1, space2);
        end
        
        function bool = gt(space1, space2)
            bool = ~le(space1, space2);
        end
        
        function bool = ge(space1, space2)
            bool = le(space2, space1);
        end
        
        function bool = isequal(spaces)
            % Check whether all input spaces are equal. Spaces are considered equal if they
            % are of same size, and are equal element-wise. For convenience, empty spaces
            % are considered equal to any empty object.
            %
            % Usage
            % -----
            % bool = isequal(spaces{:})
            %
            % Repeating Arguments
            % -------------------
            % spaces : (1,:) :class:`AbstractSpace`
            %   input spaces to compare.
            %
            % Returns
            % -------
            % bool : (1,1) logical
            %   true if all inputs are equal.
            
            arguments (Repeating)
                spaces
            end
            
            if nargin > 2
                bool = isequal(spaces{1:2}) && isequal(spaces{2:end});
                return
            end
            
            bool = (isempty(spaces{1}) && isempty(spaces{2})) || ...
                (isequal(size(spaces{1}), size(spaces{2})) && ...
                all(spaces{1} == spaces{2}));
        end
        
        function hashable = GetMD5_helper(spaces)
            % Helper function for hash algorithm. This converts the space object to a data
            % structure which can be processed by :func:`GetMD5`.
            %
            % Arguments
            % ---------
            % spaces : :class:`AbstractSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % hashable : cell
            %   data which can be accepted by :func:`GetMD5`.
            
            hashable = {spaces.dimensions, spaces.dual};
        end
        
        function disp(spaces)
            % Custom display of spaces.
            
            s = settings;
            shortformat = strcmp('short', s.matlab.commandwindow.NumericFormat.ActiveValue);
                    
            if isscalar(spaces)
                fprintf('\t%s\n\n', string(spaces, 'IncludeDetails', ~shortformat));
                return
            end
            
            sz = size(spaces);
            assert(length(sz) == 2);
            dim_str = sprintf('%dx%d', sz(1), sz(2));
            fprintf('\t%s Product %s:\n', dim_str, name(spaces));
            for i = 1:length(spaces)
                subspacestring = string(spaces(i), ...
                    'IncludeType', false, 'IncludeDetails', ~shortformat);
                fprintf('\t\t%d.\t%s\n', i, subspacestring);
            end
            fprintf('\n');
        end
    end
end
