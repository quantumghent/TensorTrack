classdef (InferiorClasses = {?GradedSpace, ?CartesianSpace, ?ComplexSpace}) SumSpace < AbstractSpace
    % Direct sum structure for a tensor index.
    
    %% Constructors
    methods
        function spaces = SumSpace(subspaces)
            % Construct a direct sum of vector spaces array of vector spaces.
            %
            % Repeating Arguments
            % -------------------
            % subspace : :class:`.AbstractSpace`
            %   individual subspaces of sum space.
            %
            % Returns
            % -------
            % spaces : :class:`.SumSpace`
            %   array of spaces representing a direct sum of vector spaces.
            
            arguments (Repeating)
                subspaces
            end
            
            if nargin == 0 || (nargin == 1 && isempty(subspaces{1}))
                args = {};
            else
                dual = cell(size(subspaces));
                for i = 1:length(subspaces)
                    dual{i} = isdual(subspaces{i}(1));
                    assert(all(isdual(subspaces{i}) == dual{i}), 'Direct sum structure should have all dual or all regular spaces.');
                end
                args = [subspaces; dual];
            end
            
            spaces = spaces@AbstractSpace(args{:});
            if (nargin == 1 && isempty(subspaces{1}))
                spaces = spaces.empty(1, 0);
            end
        end
        
        function space = one(spaces)
            space = SumSpace(arrayfun(@one, spaces(1).dimensions));
        end
        
        function space = infimum(space1, space2)
            arguments
                space1 SumSpace
                space2 SumSpace
            end
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
        
        function space = mtimes(space1, space2)
            % Fuse two spaces to a single space.
            %
            % Arguments
            % ---------
            % space1, space2 : (1, 1) :class:`.SumSpace`
            %   input spaces.
            %
            % Returns
            % -------
            % space : (1, 1) :class:`.SumSpace`
            %   fused space.
            arguments
                space1 SumSpace
                space2 SumSpace
            end
            
            space = SumSpace(sum(subspaces(space1)) * sum(subspaces(space2)));
        end
    end
    
    %% Utility
    methods
        function d = dims(spaces)
            for i = length(spaces):-1:1
                d(i) = sum(dims(subspaces(spaces(i))));
            end
        end
        
        function bools = eq(spaces1, spaces2)
            assert(isequal(size(spaces1), size(spaces2)));
            bools = false(size(spaces1));
            for i = 1:numel(bools)
                bools(i) = all(subspaces(spaces1(i)) == subspaces(spaces2(i)));
            end
        end
        
        function bool = le(space1, space2)
            arguments
                space1 SumSpace
                space2 SumSpace
            end
            
            assert(isscalar(space1) && isscalar(space2));
            
            bool = le(sum(subspaces(space1)), sum(subspaces(space2)));
        end
        
        function s = subspaces(space, i)
            assert(isscalar(space), ...
                'space:argerror', 'method only defined for scalar inputs.');
            s = space.dimensions;
            if nargin > 1
                s = s(i);
            end
        end
        
        function n = nsubspaces(space)
            n = zeros(size(space));
            for i = 1:numel(n)
                n(i) = numel(space(i).dimensions);
            end
        end
        
        function [cod, dom] = slice(sumcod, sumdom, I)
            assert(length(I) == length(sumcod) + length(sumdom));
            if isempty(sumcod)
                cod = [];
            else
                for i = flip(1:length(sumcod))
                    cod(i) = subspaces(sumcod(i), I(i));
                end
            end
            if isempty(sumdom)
                dom = [];
            else
                for i = flip(1:length(sumdom))
                    dom(i) = subspaces(sumdom(i), I(end+1-i));
                end
                
            end
        end
        
        function disp(spaces)
            s = settings;
            shortformat = strcmp('short', s.matlab.commandwindow.NumericFormat.ActiveValue);
            
            if isscalar(spaces)
                subsp = subspaces(spaces);
                dimstr = num2str(sum(dims(subsp)));
                if isdual(spaces), dimstr = dimstr + "*"; end
                sz = size(subsp);
                szstr = sprintf('%dx%d', sz(1), sz(2));
                fprintf('\t%s %s: %s\n', szstr, name(spaces), dimstr);
                if ~shortformat
                    for i = 1:length(subsp)
                        fprintf('\t\t%d.\t%s\n', i, string(subsp(i), 'IncludeType', false));
                    end
                end
                fprintf('\n');
                return
            end
            
            sz = size(spaces);
            assert(length(sz) == 2);
            dim_str = sprintf('%dx%d', sz(1), sz(2));
            fprintf('\t%s Product %s:\n', dim_str, name(spaces));
            for i = 1:length(spaces)
                subsp = subspaces(spaces(i));
                dimstr = num2str(sum(dims(subsp)));
                if isdual(spaces(i)), dimstr = dimstr + "*"; end
                sz = size(subsp);
                szstr = sprintf('%dx%d', sz(1), sz(2));
                fprintf('\t\t%d.\t%s: %s\n', i, szstr, dimstr);
                
                if ~shortformat
                    for j = 1:length(subsp)
                        fprintf('\t\t\t%d.\t%s\n', j, ...
                            string(subsp(j), 'IncludeType', false));
                    end
                end
            end
            fprintf('\n');
        end
        
        function s = name(spaces)
            s = sprintf("Sum%s", name(subspaces(spaces(1))));
        end
    end
end

