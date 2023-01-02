classdef ComplexSpace < AbstractSpace
    % Tensor index without any internal structure, which is not self dual.
    
    %% Constructors
    methods
        function spaces = ComplexSpace(dimensions, isdual)
            % Construct an array of vector spaces.
            %
            % Repeating Arguments
            % -------------------
            % dimensions : (1, 1) int
            %   dimension of the vector space.
            %
            % isdual : (1, 1) logical
            %   flag to denote dual spaces.
            %
            % Returns
            % -------
            % spaces : :class:`ComplexSpace`
            %   array of complex spaces.
            
            arguments (Repeating)
                dimensions (1, 1) {mustBePositive}
                isdual (1, 1) logical
            end
            
            if nargin == 0
                args = {};
            else
                args = [dimensions; isdual];
            end
            
            spaces = spaces@AbstractSpace(args{:});
        end
    end
    
    methods (Static)
        function spaces = new(dimensions, isdual)
            % Construct a vector space. This is a utility method to be able to access the
            % constructor of a subclass.
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
            % spaces : :class:`ComplexSpace`
            %   array of spaces.
            
            arguments (Repeating)
                dimensions (1, 1)
                isdual (1, 1) logical
            end
            
            if nargin == 0
                spaces = ComplexSpace;
            else
                for i = 1:length(dimensions)
                    if isstruct(dimensions{i})
                        dimensions{i} = dimensions{i}.degeneracies;
                    end
                end
                args = [dimensions; isdual];
                spaces = ComplexSpace(args{:});
            end
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
            % spaces : (1, :) :class:`ComplexSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % c : []
            %   empty result.
            
            c = [];
        end
        
        function d = degeneracies(spaces)
            % Compute all degeneracy combinations of a product space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`ComplexSpace`
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
            % spaces : (1, :) :class:`ComplexSpace`
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
            % codomain, domain : :class:`ComplexSpace`
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
            % codomain, domain : (1, :) :class:`ComplexSpace`
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
            % codomain, domain : (1, :) :class:`ComplexSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % style : :class:`FusionStyle`
            %   fusion style of the internal structure.
            
            style = FusionStyle.Unique;
        end
        
        function space = one(~)
            space = ComplexSpace(1, false);
        end
    end
    
    
    %% Manipulations
    methods
        function space = mtimes(space, space2)
            % Fuse two spaces to a single space.
            %
            % Arguments
            % ---------
            % space1, space2 : (1, 1) :class:`ComplexSpace`
            %   input spaces.
            % 
            % Returns
            % -------
            % space : (1, 1) :class:`ComplexSpace`
            %   fused space.
            
            space.dimensions = space.dimensions * space2.dimensions;
            space.dual = false;
        end
        
        function space = prod(spaces)
            % Fuse a product space to a single space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`ComplexSpace`
            %   Array of input spaces.
            %
            % Returns
            % -------
            % space : (1, 1) :class:`ComplexSpace`
            %   Fused space which is isomorphic to the input product space.
            
            space = ComplexSpace(prod(dims(spaces)), false);
        end
        
        function space1 = infimum(space1, space2)
            assert(isscalar(space1) && isscalar(space2));
            assert(isdual(space1) == isdual(space2));
            space1.dimensions = min(dims(space1), dims(space2));
        end
    end
    
    
    %% Utility
    methods
        function bools = eq(spaces1, spaces2)
            % Verify if spaces are element-wise equal.
            %
            % Arguments
            % ---------
            % spaces1, spaces2 : (1, :) :class:`ComplexSpace`
            %   input spaces, either of equal size or scalar.
            %
            % Returns
            % -------
            % bools : (1, :) logical
            %   flags that indicate if the element spaces are equal.
            
            bools = [spaces1.dimensions] == [spaces2.dimensions] & ...
                [spaces1.isdual] == [spaces2.isdual];
        end
        
        function bool = ge(space1, space2)
            bool = le(space2, space1);
        end
        
        function bool = le(space1, space2)
            assert(isscalar(space1) && isscalar(space2));
            bool = degeneracies(space1) <= degeneracies(space2);
        end
        
        function hashable = GetMD5_helper(spaces)
            % Helper function for hash algorithm. This converts the space object to a data
            % structure which can be processed by :func:`GetMD5`.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`ComplexSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % hashable : (1, :) int
            %   data which can be accepted by :func:`GetMD5`.
            
            hashable = [spaces.dimensions spaces.isdual];
        end
        
        function disp(spaces)
            if isscalar(spaces)
                fprintf('%s\n', string(spaces));
                return
            end
            
            sz = size(spaces);
            dim_str = sprintf('%dx%d', sz(1), sz(2));
            fprintf('%s ProductSpace with elements:\n\n', dim_str);
            for i = 1:length(spaces)
                fprintf('%d.\t', i);
                disp(spaces(i));
                fprintf('\n');
            end
        end
        
        function s = string(spaces)
            if numel(spaces) > 1
                s = arrayfun(@string, spaces);
                return
            end
            
            if spaces.dual
                s = sprintf("Space (%d)*", dims(spaces));
            else
                s = sprintf("Space (%d)", dims(spaces));
            end
        end
    end

end