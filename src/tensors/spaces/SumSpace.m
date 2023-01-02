classdef SumSpace < AbstractSpace
    % Direct product structure for a tensor index.
    
    %% Constructors
    methods
        function spaces = SumSpace(subspaces)
            % Construct an array of vector spaces.
            
            arguments (Repeating)
                subspaces
            end
            
            if nargin == 0
                args = {};
            else
                dual = cell(size(subspaces));
                for i = 1:length(subspaces)
                    dual{i} = isdual(subspaces{i}(1));
                    assert(all(isdual(subspaces{i}) == dual{i}), 'Direct product structure should have all dual or all regular spaces.');
                end
                args = [subspaces; dual];
            end
            
            spaces = spaces@AbstractSpace(args{:});
        end
        
        function space = one(spaces)
            space = SumSpace(arrayfun(@one, spaces(1).dimensions));
        end
        
        function space = infimum(space1, space2)
            assert(isscalar(space1) && isscalar(space2));
            assert(isdual(space1) == isdual(space2));
            space = SumSpace(arrayfun(@infimum, space1.dimensions, space2.dimensions));
        end
    end
    
    
    %% Manipulations
    methods
        function spaces = conj(spaces)
            for i = 1:length(spaces)
                spaces(i).dual = ~spaces(i).dual;
                spaces(i).dimensions = conj(subspaces(spaces(i)));
            end
        end
    end
    
    %% Utility
    methods
        function bools = eq(spaces1, spaces2)
            assert(isequal(size(spaces1), size(spaces2)));
            bools = false(size(spaces1));
            for i = 1:numel(bools)
                bools(i) = all(subspaces(spaces1(i)) == subspaces(spaces2(i)));
            end
        end
        
        function s = subspaces(space, i)
            assert(isscalar(space), ...
                'space:argerror', 'method only defined for scalar inputs.');
            s = space.dimensions;
            if nargin > 1
                s = s(i);
            end
        end
    end
end

