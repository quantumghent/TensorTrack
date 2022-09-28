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
                spaces = GradedSpace(varargin{:});
                return
            end
            
            if isa(varargin{1}, 'AbstractCharge')
                assert(mod(nargin, 3) == 0);
                args = cell(2, nargin / 3);
                for i = 1:size(args, 2)
                    args{1, i} = struct('charges', varargin{3 * i - 2}, ...
                        'degeneracies', varargin{3 * i - 1});
                    args{2, i} = varargin{3 * i};
                end
                spaces = GradedSpace(args{:});
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
        
        function trees = fusiontrees(codomain, domain)
            % Compute all allowed fusiontrees that connect domain and codomain. Only the
            % trivial fusion tree is allowed, so this returns empty.
            %
            % Arguments
            % ---------
            % codomain, domain : :class:`GradedSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % trees : :class:`FusionTree`
            %   list of all allowed fusion trees.
            
            rank = [length(codomain) length(domain)];
            spaces = [codomain flip(domain)];
            
            args = cell(2, sum(rank));
            for i = 1:size(args, 2)
                args{1, i} = charges(spaces(i));
                args{2, i} = spaces(i).dual;
            end
            
            trees = FusionTree.new(rank, args{:});
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
        
        function space = prod(spaces)
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
            
            space = GradedSpace(newdimensions, false);
        end
        
        function spaces = insertone(spaces, i, dual)
            arguments
                spaces
                i = length(spaces)
                dual = false
            end
            
            trivialspace = one(spaces);
            if dual, trivialspace = conj(trivialspace); end
            spaces = [spaces(1:i) trivialspace spaces(i+1:end)];
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
            assert(isscalar(space1) && isscalar(space2));
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
        
        function bool = hascharge(spaces, charges)
            mustBeEqualLength(spaces, charges);
            bool = false(1, length(space));
            for i = 1:length(space)
                bool(i) = any(space.charges{1} == charge(i));
            end
        end
        
        function space = fuse(spaces)
            space = spaces(1);
            if length(spaces) == 1
                return
            end
            
            allcharges = combcharges(spaces, []).';
            degeneracies = computedegeneracies(spaces, [], allcharges);
            
            switch fusionstyle(allcharges)
                case FusionStyle.Unique
                    fusedcharges = prod(allcharges, 2);
                    fusedbonds = prod(degeneracies, 2);
                    
                    [newcharges, ~, ic] = unique(fusedcharges);
                    newbonds = zeros(size(newcharges));
                    for i = 1:length(newcharges)
                        newbonds(i) = sum(fusedbonds(ic == i));
                    end
                    space = GradedSpace(newcharges, newbonds, false);
                    
                otherwise
                    fusedbonds = prod(degeneracies, 2);
                    fusedcharges_cell = cell(size(allcharges, 1), 1);
                    fusedbonds_cell = cell(size(fusedcharges_cell));
                    for i = 1:size(allcharges, 1)
                        [fusedcharges_cell{i}, N] = prod(allcharges(i, :));
                        fusedbonds_cell{i} = fusedbonds(i) .* N;
                    end
                    
                    fusedbonds = horzcat(fusedbonds_cell{:});
                    [newcharges, ~, ic] = unique(horzcat(fusedcharges_cell{:}));
                    newbonds = zeros(size(newcharges));
                    for i = 1:length(newcharges)
                        newbonds(i) = sum(fusedbonds(ic == i));
                    end
                    space = GradedSpace(newcharges, newbonds, false);
                    
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
        function disp(spaces)
            % Custom display of spaces.
            if isscalar(spaces)
                fprintf('%s', string(spaces));
%                 fprintf('%s GradedSpace of dimension %d:\n', ...
%                     class(spaces.dimensions.charges), dims(spaces));
%                 title_str = strjust(pad(["dual", "charges", "degeneracies"]), 'right');
%                 charge_str = strjust(pad([string(spaces.dimensions.charges)
%                     string(spaces.dimensions.degeneracies)]), 'center');
%                 fprintf('\t%s:\t%s\n', title_str(1), string(spaces.dual));
%                 fprintf('\t%s:\t%s\n', title_str(2), join(charge_str(1, :), char(9)));
%                 fprintf('\t%s:\t%s\n', title_str(3), join(charge_str(2, :), char(9)));
                return
            end
            
            sz = size(spaces);
            assert(length(sz) == 2);
            dim_str = sprintf('%dx%d', sz(1), sz(2));
            fprintf('%s ProductSpace with elements:\n\n', dim_str);
            for i = 1:length(spaces)
                fprintf('%d.\t', i);
                disp(spaces(i));
                fprintf('\n');
            end
        end
        
        function s = string(spaces)
            title_str = strjust(pad(["dual", "charges", "degeneracies"]), 'right');
            charge_str = strjust(pad([string(spaces.dimensions.charges)
                string(spaces.dimensions.degeneracies)]), 'center');
            s = sprintf(...
                '%s GradedSpace of dimension %d:\n\t%s:\t%s\n\t%s:\t%s\n\t%s:\t%s\n', ...
                class(spaces.dimensions.charges), dims(spaces), ...
                title_str(1), string(spaces.dual), ...
                title_str(2), join(charge_str(1, :), char(9)), ...
                title_str(3), join(charge_str(2, :), char(9)));
        end
    end
    
    
    
end
