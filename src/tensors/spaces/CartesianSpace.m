classdef CartesianSpace < AbstractSpace
    % Tensor index without any internal structure, which is self dual.
    
    %% Constructors
    methods
        function spaces = CartesianSpace(dimensions, ~)
            % Construct an array of vector spaces.
            %
            % Repeating Arguments
            % -------------------
            % dimensions : (1, 1) int
            %   dimension of the vector space.
            %
            % ~ : any
            %   ignored argument for syntax congruency with other spacetypes.
            %
            % Returns
            % -------
            % spaces : :class:`CartesianSpace`
            %   array of cartesian spaces.
            
            arguments (Repeating)
                dimensions (1,1) {mustBePositive}
                ~
            end
            
            if nargin == 0
                args = {};
            else
                isdual = num2cell(false(size(dimensions)));
                args = [dimensions; isdual];
            end
            
            spaces = spaces@AbstractSpace(args{:});
        end
    end
    
    methods (Static)
        function spaces = new(varargin)
            % Construct a vector space. This is a utility method to be able to access the
            % constructor of a subclass.
            %
            % Usage
            % -----
            % spaces = CartesianSpace.new(dimensions, ~, ...)
            %
            % spaces = CartesianSpace.new(dims)
            %
            % Repeating Arguments
            % -------------------
            % dimensions : struct
            %   a variable which represents the internal dimension of the space.
            %
            % dims : (1, :) int
            %   vector of dimensions of all spaces.
            %
            % Returns
            % -------
            % spaces : :class:`CartesianSpace`
            %   array of spaces.
            
            if isstruct(varargin{1})
                assert(mod(nargin, 2) == 0)
                args = cell(2, nargin / 2);
                for i = 1:nargin / 2
                    args{1, i} = varargin{2 * i - 1}.degeneracies;
                end
                spaces = CartesianSpace(args{:});
                return
            end
            
            if nargin == 1 && isnumeric(varargin{1})
                args = cell(2, length(varargin{1}));
                args(1, :) = num2cell(varargin{1});
                spaces = CartesianSpace(args{:});
                return
            end
            
            error('Undefined syntax.');
        end
    end
    
    
    %% Structure
    methods
        function c = charges(~)
            % Compute all charge combinations of a space. No internal structure is present, 
            % this yields an empty result.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`CartesianSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % c : []
            %   empty result.
            
            c = Z1;
        end
        
        function d = degeneracies(spaces)
            % Compute all degeneracy combinations of a product space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`CartesianSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % d : (1, :) int
            %   list of degeneracy combinations, with 1 element.
            
            d = [spaces.dimensions];
        end
        
        function d = dims(spaces)
            % Compute the dimension of the spaces.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`CartesianSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % d : (1, :) int
            %   total dimension of each of the input spaces.
            
            d = [spaces.dimensions];
        end
        
        function trees = fusiontrees(~, ~)
            % Compute all allowed fusiontrees that connect domain and codomain. Only the
            % trivial fusion tree is allowed, so this returns empty.
            %
            % Arguments
            % ---------
            % codomain, domain : :class:`CartesianSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % trees : :class:`FusionTree`
            %   only the trivial tree is allowed.
            
            trees = FusionTree();
        end
        
        function style = braidingstyle(~, ~)
            % Determine the braiding style of the internal structure of a space.
            %
            % Arguments
            % ---------
            % codomain, domain : (1, :) :class:`CartesianSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % style : :class:`BraidingStyle`
            %   Trivial braiding style
            
            style = BraidingStyle.Abelian;
        end
        
        function style = fusionstyle(~, ~)
            % Determine the fusion style of the internal structure of a space.
            %
            % Arguments
            % ---------
            % codomain, domain : (1, :) :class:`CartesianSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % style : :class:`FusionStyle`
            %   fusion style of the internal structure.
            
            style = FusionStyle.Unique;
        end
        
        function space = one(~)
            space = CartesianSpace(1, []);
        end
        
        function spaces = insertone(spaces, i, ~)
            arguments
                spaces
                i = length(spaces) + 1
                ~
            end
            
            spaces = [spaces(1:i-1) one(spaces) spaces(i:end)];
        end
    end
        
        
    %% Manipulations
    methods
        function spaces = conj(spaces)
            % Compute the element-wise dual space, which is equal to itself.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`CartesianSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % spaces : (1, :) :class:`CartesianSpace`
            %   dual spaces.
        end
        
        function space = mtimes(space, space2)
            % Fuse two spaces to a single space.
            %
            % Arguments
            % ---------
            % space1, space2 : (1, 1) :class:`CartesianSpace`
            %   input spaces.
            % 
            % Returns
            % -------
            % space : (1, 1) :class:`CartesianSpace`
            %   fused space.
            
            space.dimensions = space.dimensions * space2.dimensions;
        end
        
        function space = prod(spaces)
            % Fuse a product space to a single space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`CartesianSpace`
            %   Array of input spaces.
            %
            % Returns
            % -------
            % space : (1, 1) :class:`CartesianSpace`
            %   Fused space which is isomorphic to the input product space.
            
            space = CartesianSpace(prod(dims(spaces)), []);
        end
        
        function space = plus(space, space2)
            space.dimensions = space.dimensions + space2.dimensions;
        end
        
        function space = sum(spaces)
            
            space = CartesianSpace(sum(dims(spaces)), []);
        end
    end
        
        
    %% Utility
    methods
        function bools = eq(spaces1, spaces2)
            % Verify if spaces are element-wise equal.
            %
            % Arguments
            % ---------
            % spaces1, spaces2 : (1, :) :class:`CartesianSpace`
            %   input spaces, either of equal size or scalar.
            %
            % Returns
            % -------
            % bools : (1, :) logical
            %   flags that indicate if the element spaces are equal.
            
            bools = [spaces1.dimensions] == [spaces2.dimensions];
        end
        
        function bool = ge(space1, space2)
            bool = le(space2, space1);
        end
        
        function bool = le(space1, space2)
            assert(isscalar(space1) && isscalar(space2));
            bool = degeneracies(space1) <= degeneracies(space2);
        end
        
        function space = infimum(space1, space2)
            assert(isscalar(space1) && isscalar(space2));
            space = CartesianSpace.new(min(dims(space1), dims(space2)));
        end
        
        function hashable = GetMD5_helper(spaces)
            % Helper function for hash algorithm. This converts the space object to a data
            % structure which can be processed by :func:`GetMD5`.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`CartesianSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % hashable : (1, :) int
            %   data which can be accepted by :func:`GetMD5`.
            
            hashable = [spaces.dimensions];
        end
        
        function disp(spaces)
            if isscalar(spaces)
                fprintf('CartesianSpace of dimension %d\n', spaces.dimensions);
                return
            end
            
            fprintf('%dx%d Product space with elements:\n\n', ...
                size(spaces, 1), size(spaces, 2));
            for i = 1:length(spaces)
                fprintf('%d.\t', i);
                disp(spaces(i));
                fprintf('\n');
            end
        end
        
        function s = string(spaces)
            if numel(spaces) > 1
                s = strings(size(spaces));
                for i = 1:numel(spaces)
                    s(i) = string(spaces(i));
                end
                return
            end
            
            s = sprintf('CartesianSpace of dimension %d', spaces.dimensions);
        end
    end    
end
