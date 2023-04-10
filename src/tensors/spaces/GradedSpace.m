classdef GradedSpace < AbstractSpace
    % Tensor index structure with internal structure, which is not self dual.
    
    %% Constructors
    methods
        function spaces = GradedSpace(dimensions, dual)
            % Construct an array of vector spaces.
            %
            % Repeating Arguments
            % -------------------
            % dimensions : (1, 1) struct
            %   internal structure of the vector space, with fields 'charges' and
            %   'degeneracies'.
            %
            % dual : (1, 1) logical
            %   flag to denote dual spaces.
            %
            % Returns
            % -------
            % spaces : :class:`GradedSpace`
            %   array of cartesian spaces.
            
            arguments (Repeating)
                dimensions (1, 1) struct
                dual (1, 1) logical
            end
            
            if nargin == 0
                args = {};
            else
                for i = 1:length(dimensions)
                    assert(isfield(dimensions{i}, 'charges'));
                    assert(isfield(dimensions{i}, 'degeneracies'));
                    assert(length(dimensions{i}.charges) == ...
                        length(dimensions{i}.degeneracies), ...
                        'tensors:ArgumentError', ...
                        'Charges and degeneracies should have equal size.');
                    assert(length(unique(dimensions{i}.charges)) == ...
                        length(dimensions{i}.charges), ...
                        'tensors:ArgumentError', ...
                        'Charges should be unique.');
                    
                    assert(all(dimensions{i}.degeneracies >= 0), 'tensors:argerror', ...
                        'degeneracies should be positive.');
                    mask = dimensions{i}.degeneracies == 0;
                    if any(mask)
                        warning('degeneracies should be strictly positive, removing 0.');
                        dimensions{i}.charges = dimensions{i}.charges(mask);
                        dimensions{i}.degeneracies = dimensions{i}.degeneracies(mask);
                    end
                    [dimensions{i}.charges, I] = sort(dimensions{i}.charges);
                    dimensions{i}.degeneracies = dimensions{i}.degeneracies(I);
                end
                
                args = [dimensions; dual];
            end
            
            spaces = spaces@AbstractSpace(args{:});
        end
        
        function space = one(spaces)
            space = GradedSpace(...
                struct('charges', one(spaces(1).dimensions.charges), 'degeneracies', 1), ...
                false);
        end
        
        function space = infimum(space1, space2)
            assert(isscalar(space1) && isscalar(space2));
            assert(isdual(space1) == isdual(space2));
            
            [dims.charges, ia, ib] = intersect(charges(space1), charges(space2));
            d1 = degeneracies(space1);
            d2 = degeneracies(space2);
            dims.degeneracies = min(d1(ia), d2(ib));
            
            space = GradedSpace.new(dims, isdual(space1));
        end
    end
    
    methods (Static)
        function spaces = new(varargin)
            % Construct a vector space. This is a utility method to be able to access the
            % constructor of a subclass.
            %
            % Usage
            % -----
            % spaces = GradedSpace.new(charges, degeneracies, dual, ...)
            %
            % spaces = GradedSpace.new(dimensions, dual, ...)
            %
            % Repeating Arguments
            % -------------------
            % dimensions : struct
            %   a variable which represents the internal dimension of the space.
            %
            % charges : :class:`AbstractCharge`
            %   charges for the internal structure of the space.
            %
            % degeneracies : int
            %   degeneracies for the internal structure of the space.
            %
            % dual : logical
            %   a variable which indicates if a space is dual.
            %
            % Returns
            % -------
            % spaces : :class:`GradedSpace`
            %   array of spaces.
            
            if isstruct(varargin{1})    % default 
                assert(mod(nargin, 2) == 0);
                for i = 1:2:nargin
                    if varargin{i+1}, varargin{i}.charges = conj(varargin{i}.charges); end
                end
                spaces = GradedSpace(varargin{:});
                return
            end
            
            if isa(varargin{1}, 'AbstractCharge')
                assert(mod(nargin, 3) == 0, 'Unknown caller syntax');
                args = cell(2, nargin / 3);
                for i = 1:size(args, 2)
                    args{1, i} = struct('charges', varargin{3 * i - 2}, ...
                        'degeneracies', varargin{3 * i - 1});
                    args{2, i} = varargin{3 * i};
                end
                spaces = GradedSpace.new(args{:});
                return
            end
            
            error('Undefined syntax.');
        end
    end
    
    
    %% Structure
    methods
        function c = charges(spaces)
            % Compute all charge combinations of a product space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`GradedSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % c : (:, :) :class:`AbstractCharge`
            %   list of charge combinations, where each row is an entry.
            
            if isdual(spaces(1))
                c = conj(spaces(1).dimensions.charges);
            else
                c = spaces(1).dimensions.charges;
            end
            
            for i = 2:length(spaces)
                if isdual(spaces(i))
                    c = combvec(c, conj(spaces(i).dimensions.charges));
                else
                    c = combvec(c, spaces(i).dimensions.charges);
                end
            end
        end
        
        function d = degeneracies(spaces)
            % Compute all degeneracy combinations of a product space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`GradedSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % d : (:, :) int
            %   list of degeneracy combinations, where each row is an entry.
            
            d = spaces(1).dimensions.degeneracies;
            for i = 2:length(spaces)
                d = combvec(d, spaces(i).dimensions.degeneracies);
            end
        end
        
        function d = dims(spaces)
            % Compute the dimension of the spaces.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`GradedSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % d : (1, :) int
            %   total dimension of each of the input spaces.
            
            d = zeros(size(spaces));
            for i = 1:length(spaces)
                d(i) = sum(qdim(spaces(i).dimensions.charges) .* ...
                    spaces(i).dimensions.degeneracies);
            end
        end
        
        %% Space manipulations
        function space = mtimes(space1, space2)
            % Fuse two spaces to a single space.
            %
            % Arguments
            % ---------
            % space1, space2 : (1, 1) :class:`GradedSpace`
            %   input spaces.
            % 
            % Returns
            % -------
            % space : (1, 1) :class:`GradedSpace`
            %   fused space.
            
            space = prod([space1, space2]);
        end
        
        function space = prod(spaces, isdual)
            % Fuse a product space to a single space.
            %
            % Arguments
            % ---------
            % spaces : (1, :) :class:`GradedSpace`
            %   Array of input spaces.
            %
            % Returns
            % -------
            % space : (1, 1) :class:`GradedSpace`
            %   Fused space which is isomorphic to the input product space.
            arguments
                spaces
                isdual = false
            end
            if fusionstyle(spaces) == FusionStyle.Unique
                c = prod(charges(spaces), 1);
                d = prod(degeneracies(spaces), 1);
            else
                c = charges(spaces).';
                d = prod(degeneracies(spaces), 1);
                c_cell = cell(size(c, 1), 1);
                d_cell = cell(size(c_cell));
                for i = 1:length(c_cell)
                    [c_cell{i}, N] = prod(c(i, :));
                    d_cell{i} = d(i) .* N;
                end
                c = horzcat(c_cell{:});
                d = horzcat(d_cell{:});
            end
            
            newdimensions = struct;
            [newdimensions.charges, ~, ic] = unique(c);
            newdimensions.degeneracies = zeros(size(newdimensions.charges));
            for i = 1:length(newdimensions.degeneracies)
                idx = ic == i;
                newdimensions.degeneracies(i) = sum(d(idx));
            end
            
            space = GradedSpace.new(newdimensions, isdual);
        end
        
        
        %% Utility
        function bools = eq(spaces1, spaces2)
            % Verify if spaces are element-wise equal.
            %
            % Arguments
            % ---------
            % spaces1, spaces2 : (1, :) :class:`GradedSpace`
            %   input spaces, either of equal size or scalar.
            %
            % Returns
            % -------
            % bools : (1, :) logical
            %   flags that indicate if the element spaces are equal.
            
            if isempty(spaces1) && isempty(spaces2)
                bools = [];
                return
            end
            
            bools = [spaces1.dual] == [spaces2.dual];
            
            if isscalar(spaces2)
                for i = 1:length(spaces1)
                    bools(i) = bools(i) && ...
                        isequal(spaces1(i).dimensions, spaces2.dimensions);
                end
            elseif isscalar(spaces1)
                for i = 1:length(spaces2)
                    bools(i) = bools(i) && ...
                        isequal(spaces2(i).dimensions, spaces1.dimensions);
                end
            else
                assert(isequal(size(spaces1), size(spaces2)));
                for i = 1:length(spaces1)
                    bools(i) = bools(i) && ...
                        isequal(spaces1(i).dimensions, spaces2(i).dimensions);
                end
            end
        end
        
        function bools = ne(spaces1, spaces2)
            bools = ~(spaces1 == spaces2);
        end
        
        function bool = ge(space1, space2)
            bool = le(space2, space1);
        end
        
        function bool = le(space1, space2)
            assert(isscalar(space1) && isscalar(space2), 'spaces:scalar', ...
                'method only defined for scalar inputs.');
            [lia, locb] = ismember(charges(space1), charges(space2));
            if ~all(lia)
                bool = false;
                return
            end
            d1 = degeneracies(space1);
            d2 = degeneracies(space2);
            bool = all(d1 <= d2(locb));
        end
        
        function style = braidingstyle(codomain, domain)
            if nargin == 1 || ~isempty(codomain)
                style = braidingstyle(codomain(1).charges);
            else
                style = braidingstyle(domain(1).charges);
            end
        end
        
        function style = fusionstyle(codomain, domain)
            if nargin == 1 || ~isempty(codomain)
                style = fusionstyle(codomain(1).charges);
            else
                style = fusionstyle(domain(1).charges);
            end
        end
        
        function space = plus(space, space2)
            assert(isscalar(space) && isscalar(space2));
            assert(space.dual == space2.dual);
            
            c = charges(space2);
            d = degeneracies(space2);
            [lia, locb] = ismember(c, charges(space));
            for i = find(lia)'
                space.dimensions.degeneracies(locb(i)) = d(i) + ...
                    space.dimensions.degeneracies(locb(i));
            end
            [space.dimensions.charges, p] = sort([space.dimensions.charges, ...
                c(~lia)]);
            space.dimensions.degeneracies = [space.dimensions.degeneracies, ...
                d(~lia)];
            space.dimensions.degeneracies = space.dimensions.degeneracies(p);
        end
    end

    methods
        function s = string(spaces, kwargs)
            arguments
                spaces
                kwargs.IncludeType = true
                kwargs.IncludeDetails = true
            end
            
            if numel(spaces) > 1
                kwargs = namedargs2cell(kwargs);
                s = arrayfun(@(x) string(x, kwargs{:}), spaces);
                return
            end
            
            dimstring = sprintf("%d", dims(spaces));
            if isdual(spaces), dimstring = dimstring + "*"; end
            
            if kwargs.IncludeType
                typestring = name(spaces);
            end
            
            if kwargs.IncludeDetails
                chargestring = "(" + join(compose("%s => %d", ...
                    string(spaces.dimensions.charges).', ...
                    spaces.dimensions.degeneracies.'), ', ') + ")";
            end
            
            if kwargs.IncludeType && kwargs.IncludeDetails
                s = sprintf("%s: %s %s", typestring, dimstring, chargestring);
            elseif kwargs.IncludeType
                s = sprintf("%s: %s", typestring, dimstring);
            elseif kwargs.IncludeDetails
                s = sprintf("%s %s", dimstring, chargestring);
            else
                s = dimstring;
            end
        end
        
        function complexspaces = ComplexSpace(gradedspaces)
            d = num2cell(dims(gradedspaces));
            isdual = num2cell([gradedspaces.dual]);
            args = [d; isdual];
            complexspaces = ComplexSpace(args{:});
        end
        
        function cartesianspaces = CartesianSpace(gradedspaces)
            d = num2cell(dims(gradedspaces));
            cartesianspaces = CartesianSpace(d{:});
        end
        
        function s = name(spaces)
            s = sprintf("%sSpace", name(spaces(1).dimensions.charges));
        end
    end
    
    
    
end
