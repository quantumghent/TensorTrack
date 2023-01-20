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
                    (all(p == 1:length(p)) && all(X.rank == Y.rank))
                
                Y = axpby@MatrixBlock(a, X, b, Y);
                return;
            end
            
            
            %% General case: addition with permutation
            % tensor indexing to matrix indexing
            rx = X.rank;
            ry = Y.rank;
            rrx = rankrange(rx);
            rry = rankrange(ry);
            p_eff(rry) = rrx(p);
            p_eff_1 = p_eff(1:ry(1));
            p_eff_2 = p_eff((1:ry(2)) + ry(1));
            
            % extract small tensor blocks
            doA = a ~= 1;
            vars_in = cell(size(map, 1), 1);
            offset = 0;
            
            Xrowsizes = X.rowsizes;
            Xcolsizes = X.colsizes;
            Xtdims = X.tdims;
            Xvar = X.var;
            for i = 1:length(Xvar)
                rowsz = Xrowsizes{i};
                colsz = Xcolsizes{i};
                
                var_in = Xvar{i};
                if doA, var_in = a .* var_in; end
                
                oldtdims = Xtdims{i}(:, rrx);
                newmdims = [prod(oldtdims(:, p_eff_1), 2) prod(oldtdims(:, p_eff_2), 2)];
                
                for k = 1:length(colsz) - 1
                    for j = 1:length(rowsz) - 1
                        ind = j + (k-1) * (length(rowsz) - 1);
                        offset = offset + 1;
                        vars_in{offset} = reshape(permute(reshape(...
                            var_in(rowsz(j)+1:rowsz(j+1), colsz(k)+1:colsz(k+1)), ...
                            oldtdims(ind, :)), ...
                            p_eff), ...
                            newmdims(ind, :));
                    end
                end
            end
            
            % apply map
            vars_out = cell(size(map, 2), 1);
            [rows, cols, vals] = find(map);
            for i = 1:length(vals)
                vars_out{cols(i)} = vals(i) * vars_in{rows(i)};
            end
            
            % inject small tensor blocks
            if b == 0
                offset = 0;
                for i = 1:length(Y.var)
                    rows = length(Y.rowsizes{i}) - 1;
                    cols = length(Y.colsizes{i}) - 1;
                    if rows < cols
                        m = cell(rows, 1);
                        for n = 1:rows
                            m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                        end
                        Y.var{i} = cat(1, m{:});
                    else
                        m = cell(cols, 1);
                        for n = 1:cols
                            m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                        end
                        Y.var{i} = cat(2, m{:});
                    end
                    offset = offset + rows * cols;
                end
                return
            end
            
            if b == 1
                offset = 0;
                for i = 1:length(Y.var)
                    rows = length(Y(i).rowsizes) - 1;
                    cols = length(Y(i).colsizes) - 1;
                    if rows < cols
                        m = cell(rows, 1);
                        for n = 1:rows
                            m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                        end
                        Y.var{i} = Y.var{i} + cat(1, m{:});
                    else
                        m = cell(cols, 1);
                        for n = 1:cols
                            m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                        end
                        Y.var{i} = Y.var{i} + cat(2, m{:});
                    end
                    offset = offset + rows * cols;
                end
                return
            end
            
            offset = 0;
            for i = 1:length(Y.var)
                rows = length(Y.rowsizes{i}) - 1;
                cols = length(Y.colsizes{i}) - 1;
                if rows < cols
                    m = cell(rows, 1);
                    for n = 1:rows
                        m{n} = cat(2, vars_out{offset + n + ((1:cols)-1) * rows});
                    end
                    Y.var{i} = b .* Y.var{i} + cat(1, m{:});
                else
                    m = cell(cols, 1);
                    for n = 1:cols
                        m{n} = cat(1, vars_out{offset + (n-1) * rows + (1:rows)});
                    end
                    Y.var{i} = b .* Y.var{i} + cat(2, m{:});
                end
                offset = offset + rows * cols;
            end
        end

        function typename = underlyingType(b)
            typename = underlyingType(b.var{1});
        end
    end
end
