classdef AbelianBlock < MatrixBlock
    %ABELIANBLOCK Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods
        function Y = axpby(a, X, b, Y, p, map)
            %% Special cases
            % only overload general case, other cases are unchanged.
            if a == 0 || ...
                    nargin == 4 || ...
                    (isempty(p) && isempty(map)) || ...
                    (all(p == 1:length(p)) && all(X(1).rank == Y(1).rank))
                
                Y = axpby@MatrixBlock(a, X, b, Y);
                return;
            end
            
            
            %% General case: addition with permutation
            % tensor indexing to matrix indexing
            rrx = rankrange(X(1).rank);
            rry = rankrange(Y(1).rank);
            p_eff(rry) = rrx(p);
            
            % extract small tensor blocks
            doA = a ~= 1;
            vars = cell(size(map, 1), 1);
            offset = 0;
            for i = 1:length(X)
                rowsz = X(i).rowsizes;
                colsz = X(i).colsizes;
                
                var_in = X(i).var;
                if doA, var_in = a .* var_in; end
                
                oldtdims = X(i).tdims(:, rrx);
                newmdims = [prod(X(i).tdims(:, p_eff(1:Y(1).rank(1))), 2) ...
                    prod(X(i).tdims(:, p_eff((1:Y(1).rank(2)) + Y(1).rank(1))), 2)];
                
                for k = 1:length(colsz) - 1
                    for j = 1:length(rowsz) - 1
                        ind = j + (k-1) * (length(rowsz) - 1);
                        offset = offset + 1;
                        
                        vars{offset} = reshape(permute(reshape(...
                            var_in(rowsz(j)+1:rowsz(j+1), colsz(k)+1:colsz(k+1)), ...
                            oldtdims(ind, :)), ...
                            p_eff), ...
                            newmdims(ind, :));
                    end
                end
            end
            
            % apply map
            [row, col] = find(map);
            vars(col) = vars(row);
            
            % inject small tensor blocks
            if b == 0
                offset = 0;
                for i = 1:length(Y)
                    rows = length(Y(i).rowsizes) - 1;
                    cols = length(Y(i).colsizes) - 1;
                    if rows < cols
                        m = cell(rows, 1);
                        for n = 1:rows
                            m{n} = cat(2, vars{offset + n + ((1:cols)-1) * rows});
                        end
                        Y(i).var = cat(1, m{:});
                    else
                        m = cell(cols, 1);
                        for n = 1:cols
                            m{n} = cat(1, vars{offset + (n-1) * rows + (1:rows)});
                        end
                        Y(i).var = cat(2, m{:});
                    end
                    offset = offset + rows * cols;
                end
                return
            end
            
            if b == 1
                offset = 0;
                for i = 1:length(Y)
                    rows = length(Y(i).rowsizes) - 1;
                    cols = length(Y(i).colsizes) - 1;
                    if rows < cols
                        m = cell(rows, 1);
                        for n = 1:rows
                            m{n} = cat(2, vars{offset + n + ((1:cols)-1) * rows});
                        end
                        Y(i).var = Y(i).var + cat(1, m{:});
                    else
                        m = cell(cols, 1);
                        for n = 1:cols
                            m{n} = cat(1, vars{offset + (n-1) * rows + (1:rows)});
                        end
                        Y(i).var = Y(i).var + cat(2, m{:});
                    end
                    offset = offset + rows * cols;
                end
                return
            end
            
            offset = 0;
            for i = 1:length(Y)
                rows = length(Y(i).rowsizes) - 1;
                cols = length(Y(i).colsizes) - 1;
                if rows < cols
                    m = cell(rows, 1);
                    for n = 1:rows
                        m{n} = cat(2, vars{offset + n + ((1:cols)-1) * rows});
                    end
                    Y(i).var = b .* Y(i).var + cat(1, m{:});
                else
                    m = cell(cols, 1);
                    for n = 1:cols
                        m{n} = cat(1, vars{offset + (n-1) * rows + (1:rows)});
                    end
                    Y(i).var = b .* Y(i).var + cat(2, m{:});
                end
                offset = offset + rows * cols;
            end
            
%             rrx = rankrange(x(1).rank);
%             rry = rankrange(y(1).rank);
%             p_eff(rry) = rrx(p);
            
%             %% extract blocks
%             doA = a ~= 1;
%             vars = cell(1, size(map, 1));
%             offset = 0;
%             for block = x
%                 rowsz = block.rowsizes;
%                 colsz = block.colsizes;
%                 
%                 if doA
%                     varA = a .* block.var;
%                 else
%                     varA = block.var;
%                 end
%                 
%                 olddims = block.tdims(:, rrx);
%                 newdims = [prod(olddims(:, p_eff(1:y(1).rank(1))), 2) ...
%                     prod(olddims(:, p_eff(y(1).rank(1) + (1:y(1).rank(2)))), 2)];
%                 for k = 1:length(colsz) - 1
%                     for j = 1:length(rowsz) - 1
%                         ind = j + (k-1) * (length(rowsz) - 1);
%                         offset = offset + 1;
%                         %                         vars{offset} = varA(rowsz(j)+1:rowsz(j+1), colsz(k)+1:colsz(k+1));
%                         %                         vars{offset} = reshape(vars{offset}, olddims(ind, :));
%                         %                         vars{offset} = permute(vars{offset}, p_eff);
%                         %                         vars{offset} = reshape(vars{offset}, newdims(ind, :));
%                         
%                         % less readible but faster
%                         vars{offset} = reshape(permute(reshape(...
%                             varA(rowsz(j)+1:rowsz(j+1), colsz(k)+1:colsz(k+1)), ...
%                             olddims(ind, :)), p_eff), ...
%                             newdims(ind, :));
%                     end
%                 end
%             end
%             
%             %% apply map
%             [row, col] = find(map);
%             vars(col) = vars(row);
%             
%             %% inject blocks
%             if b == 0
%                 offset = 0;
%                 for i = 1:length(y)
%                     rows = length(y(i).rowsizes) - 1;
%                     cols = length(y(i).colsizes) - 1;
%                     if rows < cols
%                         m = cell(rows, 1);
%                         for n = 1:rows
%                             m{n} = cat(2, vars{offset + n + ((1:cols)-1) * rows});
%                         end
%                         y(i).var = cat(1, m{:});
%                     else
%                         m = cell(cols, 1);
%                         for n = 1:cols
%                             m{n} = cat(1, vars{offset + (n-1) * rows + (1:rows)});
%                         end
%                         y(i).var = cat(2, m{:});
%                     end
%                     offset = offset + rows * cols;
%                 end
%                 return
%             end
%             
%             if b == 1
%                 offset = 0;
%                 for i = 1:length(y)
%                     rows = length(y(i).rowsizes) - 1;
%                     cols = length(y(i).colsizes) - 1;
%                     if rows < cols
%                         m = cell(rows, 1);
%                         for n = 1:rows
%                             m{n} = cat(2, vars{offset + n + ((1:cols)-1) * rows});
%                         end
%                         y(i).var = y(i).var + cat(1, m{:});
%                     else
%                         m = cell(cols, 1);
%                         for n = 1:cols
%                             m{n} = cat(1, vars{offset + (n-1) * rows + (1:rows)});
%                         end
%                         y(i).var = y(i).var + cat(2, m{:});
%                     end
%                     offset = offset + rows * cols;
%                 end
%                 return
%             end
%             
%             offset = 0;
%             for i = 1:length(y)
%                 rows = length(y(i).rowsizes) - 1;
%                 cols = length(y(i).colsizes) - 1;
%                 if rows < cols
%                     m = cell(rows, 1);
%                     for n = 1:rows
%                         m{n} = cat(2, vars{offset + n + ((1:cols)-1) * rows});
%                     end
%                     y(i).var = b .* y(i).var + cat(1, m{:});
%                 else
%                     m = cell(cols, 1);
%                     for n = 1:cols
%                         m{n} = cat(1, vars{offset + (n-1) * rows + (1:rows)});
%                     end
%                     y(i).var = b .* y(i).var + cat(2, m{:});
%                 end
%                 offset = offset + rows * cols;
%             end
        end

        function typename = underlyingType(b)
            typename = underlyingType(b(1).var);
        end
    end
end
