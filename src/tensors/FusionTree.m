classdef FusionTree < matlab.mixin.CustomDisplay
    % Splitting and fusion tree pair in canonical form.
    %
    % Properties
    % ----------
    % charges : :class:`AbstractCharge`
    %   labels for the edges of the fusion trees
    %
    % vertices : int
    %   labels for the vertices of the fusion trees
    %
    % isdual : logical
    %   indicator of duality transform on the external edges.
    %
    % rank : int
    %   amount of splitting tree and fusion tree legs.
    
    properties
        charges (:, :)
        vertices (:, :) uint8
        isdual (1, :) logical
        rank (1, 2) int16
    end
    
    properties (Dependent, Hidden)
        coupled
        inner
        uncoupled
    end
    
    %% Constructors
    methods
        function f = FusionTree(charges, vertices, isdual, rank)
            % Usage
            % -----
            % :code:`f = FusionTree(charges, vertices, isdual, rank)`
            %
            % Arguments
            % ---------
            % charges : :class:`AbstractCharge`
            %   array of charges, where each row represents an allowed fusion channel.
            %
            % vertices : int = []
            %   Array of vertex labels. Must be empty or have the same amount of
            %   rows as `charges`. This is optional if the charges have a fusion style
            %   which is multiplicity-free.
            %
            % isdual : logical
            %   mask that shows the presence of duality transformations.
            %
            % rank : int
            %   number of splitting and fusing legs.
            
            if nargin == 0
                return;
            end
            
            %% Trivial tree with no external legs
            if sum(rank) == 0
                if isempty(charges) && isempty(vertices) && isempty(isdual)
                    f.rank = rank;
                    return
                end
                error('tensors:tree:argError', 'Incompatible input sizes.');
            end
            
            %% General tree
            N = length(isdual);
            if N ~= sum(rank) || ...
                    size(charges, 2) ~= sum(treeindsequence(rank)) + 1 || ...
                    (N > 2 && hasmultiplicity(fusionstyle(charges)) && ...
                    size(vertices, 2) ~= max(rank(1) - 1, 0) + max(rank(2) - 1, 0))
                error('tensors:tree:argError', 'Incompatible input sizes.');
            end
            
            f.charges = charges;
            f.vertices = vertices;
            f.isdual = isdual;
            f.rank = rank;

%             assert(all(isallowed(f)));  % comment out for debugging
        end
    end
    
    methods (Static)
        function f = new(rank, charges, isdual)
            % Generate all allowed fusion trees for a set of external charges.
            %
            % Arguments
            % ---------
            % rank : int
            %   number of legs for the splitting and fusion tree.
            %
            % Repeating Arguments
            % -------------------
            % charges : :class:`AbstractCharge`
            %   set of charges for a given leg
            %
            % isdual : logical
            %   presence of a duality transform
            
            arguments
                rank (1, 2) {mustBeNonnegative, mustBeInteger}
            end
            
            arguments (Repeating)
                charges (1, :) {mustBeA(charges, 'AbstractCharge')}
                isdual (1, 1) logical
            end
            
            assert(sum(rank) == length(isdual), 'tensors:trees:ArgError', ...
                'Rank incompatible with the amount of legs specified.');
            
            if sum(rank) == 0
                f = FusionTree([], [], [], rank);
                return
            end
            
            
            %% Splitting tree
            vertices_l = [];
            if rank(1) == 0
                charges_l = one(charges{1});
                vertices_l = uint8.empty(1, 0);
            elseif rank(1) == 1
                charges_l = [charges{1}(:) charges{1}(:)];
                vertices_l = uint8.empty(size(charges_l, 1), 0);
            else
                switch fusionstyle(charges{1})
                    case FusionStyle.Unique
                        unc_l = combvec(charges{1:rank(1)});
                        inn_l = cumprod(unc_l);
                        charges_l = reshape([unc_l(:).'; inn_l(:).'], [2 1] .* size(unc_l));
                        charges_l = charges_l(2:end, :).';
                        vertices_l = uint8.empty(size(charges_l, 1), 0);
                    case FusionStyle.Simple
                        charges_l = [];
                        for unc_l = combvec(charges{1:rank(1)})
                            inn_l = cumprod(unc_l);
                            temp = reshape( ...
                                [repmat(unc_l(:).', 1, size(inn_l, 2)); inn_l(:).'], ...
                                [2 1] .* size(inn_l));
                            charges_l = [charges_l; temp(2:end, :).']; %#ok<AGROW>
                        end
                        vertices_l = uint8.empty(size(charges_l, 1), 0);
                    case FusionStyle.Generic
                        charges_l = [];
                        vertices_l = [];
                        for unc_l = combvec(charges{1:rank(1)})
                            [inn_l, vert_l] = cumprod(unc_l);
                            temp = reshape( ...
                                [repmat(unc_l(:).', 1, size(inn_l, 2)); inn_l(:).'], ...
                                [2 1] .* size(inn_l));
                            charges_l = [charges_l; temp(2:end, :).'];  %#ok<AGROW>
                            vertices_l = [vertices_l; vert_l.'];          %#ok<AGROW>
                        end
                end
            end
            
            
            %% Fusion tree
            if rank(2) == 0
                charges_r = one(charges{1});
                vertices_r = uint8.empty(1, 0);
            elseif rank(2) == 1
                charges_r = [charges{end}(:) charges{end}(:)];
                vertices_r = uint8.empty(size(charges_r, 1), 0);
            else
                switch fusionstyle(charges{end})
                    case FusionStyle.Unique
                        unc_r = combvec(charges{end:-1:rank(1) + 1});
                        inn_r = cumprod(unc_r);
                        charges_r = reshape([unc_r(:).'; inn_r(:).'], [2 1] .* size(unc_r));
                        charges_r = charges_r(2:end, :).';
                        vertices_r = uint8.empty(size(charges_r, 1), 0);
                    case FusionStyle.Simple
                        charges_r = [];
                        for unc_r = combvec(charges{end:-1:rank(1) + 1})
                            inn_r = cumprod(unc_r);
                            temp = reshape( ...
                                [repmat(unc_r(:).', 1, size(inn_r, 2)); inn_r(:).'], ...
                                [2 1] .* size(inn_r));
                            charges_r = [charges_r; temp(2:end, :).']; %#ok<AGROW>
                        end
                        vertices_r = uint8.empty(size(charges_r, 1), 0);
                    case FusionStyle.Generic
                        charges_r = [];
                        vertices_r = [];
                        for unc_r = combvec(charges{end:-1:rank(1) + 1})
                            [inn_r, vert_r] = cumprod(unc_r);
                            temp = reshape( ...
                                [repmat(unc_r(:).', 1, size(inn_r, 2)); inn_r(:).'], ...
                                [2 1] .* size(inn_r));
                            charges_r = [charges_r; temp(2:end, :).'];  %#ok<AGROW>
                            vertices_r = [vertices_r; vert_r.'];          %#ok<AGROW>
                        end
                end
            end
            
            
            %% Combine
            coupled = unique(intersect(charges_l(:, end), charges_r(:, end)));
            if isempty(coupled)
                f = FusionTree;
                return
            end
            
            charge_cell = cell(size(coupled));
            vertices_cell = cell(size(coupled));
            for i = 1:length(coupled)
                ind_l = coupled(i) == charges_l(:, end);
                nnz_l = nnz(ind_l);
                ind_r = coupled(i) == charges_r(:, end);
                nnz_r = nnz(ind_r);
                charge_cell{i} = [reshape(repmat(reshape( ...
                    charges_l(ind_l, 1:end - 1), 1, []), nnz_r, 1), nnz_r * nnz_l, []), ...
                    repmat(coupled(i), nnz_r * nnz_l, 1), ...
                    repmat(charges_r(ind_r, end - 1:-1:1), nnz_l, 1)];
                
                if hasmultiplicity(fusionstyle(charges{1}))
                    vertices_cell{i} = [reshape(repmat(reshape( ...
                        vertices_l(ind_l, :), 1, []), nnz_r, 1), nnz_r * nnz_l, []), ...
                        repmat(vertices_r(ind_r, end:-1:1), nnz_l, 1)];
                end
            end
            
            f = FusionTree(vertcat(charge_cell{:}), vertcat(vertices_cell{:}), ...
                [isdual{:}], rank);
            f = sort(f);
        end
    end
    
    
    %% Tree manipulations
    methods
        function [c, f] = artinbraid(f, i, inv)
            % Compute the coefficients for braiding a set of two neighbouring strands.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %   tree to swap.
            %
            % i : int
            %   `i` and `i+1` are the braided strands.
            %
            % inv : logical = false
            %   flag to indicate whether to perform an overcrossing (false) or an
            %   undercrossing (true).
            %
            % Returns
            % -------
            % c : sparse double
            %   matrix of coefficients that transform input to output trees.
            %   f(i) --> c(i,j) * f(j)
            %
            % f : :class:`FusionTree`
            %   braided trees in canonical form.
            
            if nargin < 3, inv = false; end
            
            assert(f.rank(2) == 0, 'Only defined for splitting trees.');
            
            %% Special case - Abelian
            if braidingstyle(f) == BraidingStyle.Abelian
                if i == 1
                    f.charges(:, 1:2) = f.charges(:, [2 1]);
                    f.isdual(1:2) = f.isdual([2 1]);
                    [f, p] = sort(f);
                    c = sparse(1:length(f), invperm(p), ones(1, length(f)));
                    return
                end
                
                f.charges(:, [2 * i - 2 2 * i]) = f.charges(:, [2 * i 2 * i - 2]);
                f.charges(:, 2 * i - 1) = f.charges(:, 2 * i - 3) * f.charges(:, 2 * i - 2);
                f.isdual(i:i + 1) = f.isdual([i + 1 i]);
                [f, p] = sort(f);
                c = sparse(1:length(f), invperm(p), ones(1, length(f)));
                return
            end
            
            %% Special case - first legs
            %   braiding by only an R-move
            if i == 1
                if ~hasmultiplicity(fusionstyle(f))
                    R = Rsymbol(f.charges(:, 1), f.charges(:, 2), f.charges(:, 3), inv);
                    f.charges(:, 1:2) = f.charges(:, [2 1]);
                    f.isdual(1:2) = f.isdual([2 1]);
                    [f, p] = sort(f);
                    c = sparse(1:length(f), invperm(p), R);
                    return
                end
                
                charges = f.charges;
                [order, diagcharges, diagvertices] = simulsortrows(charges, f.vertices(:, 2:end));
                [~, ia1] = simulunique(diagcharges, diagvertices, 'rows','stable');
                
                
                [abc, ia2, ic2] = unique(diagcharges(ia1, 1:3), 'rows');
                blocks = cell(1, length(ia2));
                for j = 1:length(ia2)
                    blocks{j} = Rsymbol(abc(j, 1), abc(j, 2), abc(j, 3), inv);
                end

                f.charges = f.charges(order, [2 1 3:end]);
                f.vertices = f.vertices(order, :);
                f.isdual(1:2) = f.isdual([2 1]);
                c = sparse(blkdiag(blocks{ic2}));
                [f, p] = sort(f);
                c = c(invperm(order), p);
                return
            end
            
            
            %% General case
            %   braiding by R-move, F-move, R-move
            if ~hasmultiplicity(fusionstyle(f))
                charges = f.charges; %#ok<*PROPLC>
                mask = true(1, size(charges, 2)); mask(2 * i - 1) = false;
                [diagcharges, order] = sortrows(charges(:, mask));
                charges = charges(order, :);
                [~, ia1, ic1] = unique(diagcharges, 'rows');
                charges1 = charges(ia1, :);
                blocks = cell(1, length(ia1));
                newcharges = cell(1, length(ia1));
                
                [~, ia2, ic2] = unique(charges1(:, [2*i-3 2*i-2 2*i 2*i+1]), 'rows');
                
                for j = 1:length(ia2)
                    a = charges1(ia2(j), 2 * i - 3);
                    b = charges1(ia2(j), 2 * i - 2);
                    cs = charges(ic1 == ia2(j), 2 * i - 1);
                    d = charges1(ia2(j), 2 * i);
                    e = charges1(ia2(j), 2 * i + 1);
                    fs = intersect(a * d, e * conj(b)).';
                    blocks{j} = braidingmatrix(a, b, cs, d, e, fs, inv);
                    
                    for k = find(j == ic2).'
                        newcharges{k} = repmat(charges1(k, :), length(fs), 1);
                        newcharges{k}(:, 2 * i - 2:2 * i) = [newcharges{k}(:, 2 * i) fs(:) ...
                            newcharges{k}(:, 2 * i - 2)];
                    end
                end
                
                f.charges = vertcat(newcharges{:});
                f.isdual(i:i + 1) = f.isdual([i + 1 i]);
                c = sparse(blkdiag(blocks{ic2}));
                [f, p] = sort(f);
                c = c(invperm(order), p);
                return
            end
            
            charges = f.charges;
            mask1 = true(1, size(charges, 2)); mask1(2*i-1) = false;
            vertices = f.vertices;
            mask2 = true(1, size(vertices, 2)); mask2(i-1:i) = false;
            [order, diagcharges, diagvertices] = ...
                simulsortrows(charges(:, mask1), vertices(:, mask2));
            charges = charges(order, :); vertices = vertices(order, :);
            [~, ia1] = simulunique(diagcharges, diagvertices, 'rows');
            charges1 = charges(ia1, :);
            vertices1 = vertices(ia1, :);
            
            newcharges = cell(1, length(ia1));
            newvertices = cell(1, length(ia1));
            
            [~, ia2, ic2] = unique(charges1(:, [2*i-3 2*i-2 2*i 2*i+1]), 'rows');
            blocks = cell(1, length(ia2));
            for j = 1:length(ia2)
                a = charges1(ia2(j), 2 * i - 3);
                b = charges1(ia2(j), 2 * i - 2);
                d = charges1(ia2(j), 2 * i);
                e = charges1(ia2(j), 2 * i + 1);
                
                cs = intersect(a * b, conj(conj(e) * d)).';
                fs = intersect(a * d, e * conj(b)).';
                blocks{j} = braidingmatrix(a, b, cs, d, e, fs, inv);
                
                newfs = a.empty(size(blocks{j}, 2), 0);
                newverts = zeros([length(newfs), 2], 'uint8');
                ctr = 0;
                for f1 = fs.'
                    N1 = Nsymbol(a, d, f1);
                    N2 = Nsymbol(e, conj(b), f1);
                    newfs(ctr + (1:(N1*N2))) = f1;
                    newverts(ctr + (1:(N1*N2)), :) = combvec(1:N1, 1:N2).';
                    ctr = ctr + (N1*N2);
                end
                assert(size(blocks{j}, 2) == ctr);
                
                for k = find(j == ic2).'
                    newcharges{k} = repmat(charges1(k, :), length(newfs), 1);
                    newcharges{k}(:, 2 * i - 2:2 * i) = [newcharges{k}(:, 2 * i) newfs(:) ...
                        newcharges{k}(:, 2 * i - 2)];
                    newvertices{k} = repmat(vertices1(k, :), length(newfs), 1);
                    newvertices{k}(:, i-1:i) = newverts;
                end
            end
            
            f.charges = vertcat(newcharges{:});
            f.vertices = vertcat(newvertices{:});
            f.isdual(i:i + 1) = f.isdual([i + 1 i]);
            c = sparse(blkdiag(blocks{ic2}));
            [f, p] = sort(f);
            c = c(invperm(order), p);
        end
        
        function [c, f] = bendleft(f)
            % Compute the coefficients for bending a fusing leg to a splitting leg.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %   tree to bend.
            %
            % Returns
            % -------
            % c : sparse double
            %   matrix of coefficients that transform input to output trees.
            %   f(i) --> c(i,j) * f(j)
            %
            % f : :class:`FusionTree`
            %   bent trees in canonical form.
            
            f = flip(f);
            [c, f] = bendright(f);
            f = flip(f);
            c = conj(c);
            return
        end
        
        function [c, f] = bendright(f)
            % Compute the coefficients for bending a splitting leg to a fusing leg.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %   tree to bend.
            %
            % Returns
            % -------
            % c : sparse double
            %   matrix of coefficients that transform input to output trees.
            %   `f(i) --> c(i,j) * f(j)`
            %
            % f : :class:`FusionTree`
            %   bent trees in canonical form.
            
            if ~hasmultiplicity(fusionstyle(f))
                if f.rank(1) == 0
                    error('tensors:trees:ArgError', 'Invalid rank [%d %d]', f.rank(1), f.rank(2));
                elseif f.rank(1) == 1
                    [bc, ~, ic] = unique(f.charges(:, 1:2), 'rows');
                    b = bc(:, 1); a = repmat(one(b), size(b)); c = bc(:, 2);
                else
                    ind = treeindsequence(f.rank(1)) - 1:treeindsequence(f.rank(1)) + 1;
                    [abc, ~, ic] = unique(f.charges(:, ind), 'rows');
                    a = abc(:, 1); b = abc(:, 2); c = abc(:, 3);
                end
                
                vals = sqrt(qdim(c) ./ qdim(a)) .* Bsymbol(a, b, c);
                if f.isdual(f.rank(1)), vals = vals .* conj(frobeniusschur(conj(b))); end
                
                
                f.isdual(f.rank(1)) = ~f.isdual(f.rank(1));
                
                if f.rank(2) < 2
                    chargesR = [conj(f.charges(:, end - f.rank(2) - 1)) f.charges(:, end - f.rank(2) + 1:end)];
                else
                    chargesR = [conj(f.charges(:, treeindsequence(f.rank(1)))) ...
                        f.charges(:, end - treeindsequence(f.rank(2)):end)];
                end
                
                chargesC = a(ic);
                
                if f.rank(1) == 1
                    chargesL = [];
                elseif f.rank(1) == 2
                    chargesL = f.charges(:, 1);
                else
                    chargesL = f.charges(:, 1:treeindsequence(f.rank(1) - 1));
                end
                
                f.charges = [chargesL chargesC chargesR];
                
                f.rank = f.rank + int16([-1 +1]);
                c = sparse(1:length(f), 1:length(f), vals(ic));
                return
            end
            
            [~, ia1] = simulunique(f.charges, ...
                f.vertices(:, [1:f.rank(1)-2 f.rank(1):end]), 'rows', 'stable');
            if f.rank(1) == 1
                charges1 = [repmat(one(f.charges), size(ia1)) f.charges(ia1, :)];
                [abc, ia2, ic2] = unique(charges1(:, 1:3), 'rows');
            else
                ind = treeindsequence(f.rank(1)) - 1:treeindsequence(f.rank(1)) + 1;
                [abc, ia2, ic2] = unique(f.charges(ia1, ind), 'rows');
            end
            
            blocks = cell(1, length(ia2));
            for j = 1:length(ia2)
                blocks{j} = sqrt(qdim(abc(j, 3)) / qdim(abc(j, 1))) * ...
                    Bsymbol(abc(j, 1), abc(j, 2), abc(j, 3));
                if f.isdual(f.rank(1))
                    blocks{j} = blocks{j} * conj(frobeniusschur(conj(abc(j, 2))));
                end
            end
            c = sparse(blkdiag(blocks{ic2}));
            
            f.isdual(f.rank(1)) = ~f.isdual(f.rank(1));
            
            if f.rank(2) < 2
                chargesR = [conj(f.charges(:, end - f.rank(2) - 1)) f.charges(:, end - f.rank(2) + 1:end)];
            else
                chargesR = [conj(f.charges(:, treeindsequence(f.rank(1)))) ...
                    f.charges(:, end - treeindsequence(f.rank(2)):end)];
            end
            
            if f.rank(1) == 1
                chargesC = repmat(one(charges1), size(f.charges, 1), 1);
            else
                chargesC = f.charges(:, ind(1));
            end
            
            if f.rank(1) == 1
                chargesL = [];
            elseif f.rank(1) == 2
                chargesL = f.charges(:, 1);
            else
                chargesL = f.charges(:, 1:treeindsequence(f.rank(1) - 1));
            end
            
            f.charges = [chargesL chargesC chargesR];
            
            if f.rank(2) == 0 && f.rank(1) > 1
                f.vertices = f.vertices(:, 1:end-1);
            elseif f.rank(1) == 1
                f.vertices = [ones(size(f.charges, 1), 1, 'uint8') f.vertices];
            end
            
            f.rank = f.rank + int16([-1 +1]);
        end
        
        function [c, f] = braid(f, p, lvl, rank)
            % Compute the coefficients that bring a braided tree into canonical form.
            %   This is done by reducing the braid into a composition of elementary swaps
            %   on neighbouring strands.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %   tree to braid.
            %
            % p : int
            %   permutation indices.
            %
            % lvl : int
            %   height of the strands, indicating over- and undercrossings.
            %
            % rank : int
            %   final number of splitting and fusing legs.
            %
            % Returns
            % -------
            % c : sparse double
            %   matrix of coefficients that transform input to output trees.
            %   f(i) --> c(i,j) * f(j)
            %
            % f : :class:`FusionTree`
            %   braided trees in canonical form.
            %
            % Todo
            % ----
            % Currently this is done by first bending all legs to the splitting tree, this
            % step is not necessary when no splitting legs are braided with fusing legs.
            % Additionally this can be sped up significantly for charges with Abelian
            % braiding style.
            
            arguments
                f
                p (1, :) = 1:f.legs
                lvl (1, :) = 1:f.legs
                rank (1, 2) = f.rank
            end
            N = length(p);
            assert(isperm(p), 'tensors:trees:ArgError', 'p is not a valid permutation.');
            assert(N == f.legs, 'tensors:trees:ArgError', ...
                'p has the wrong number of elements.');
            
            if all(p == 1:f.legs)
                [c, f] = repartition(f, rank);
                return
            end
            
            swaps = perm2swap(p);
            [c, f] = repartition(f, [f.legs 0]);
            for s = swaps
                [c_, f] = artinbraid(f, s, lvl(s) > lvl(s + 1));
                c = c * c_;
                lvl(s:s + 1) = lvl([s + 1 s]);
            end
            [c_, f] = repartition(f, rank);
            c = c * c_;
        end
        
        function [c, f] = elementary_trace(f, i)
            N = f.legs;
            assert(N > 1, 'tensors:trees:ArgError', ...
                'Tree needs at least two legs to trace.');
            assert(f.isdual(i) ~= f.isdual(mod1(i+1, f.legs)), 'tensors:trees:ArgError', ...
                'Tracing only defined between space and its dual.');
            
            L1 = length(f);
            
            if i == 1
                mask = find(f.charges(:, 1) == conj(f.charges(:, 2)) & ...
                    f.charges(:, 3) == one(f.charges));
                b = f.charges(mask, 1);
                vals = sqrt(qdim(b));
                if f.isdual(i), vals = frobeniusschur(b) .* vals; end
                if N > 4
                    f.charges = f.charges(mask, 4:end);
                elseif N == 4
                    f.charges = f.charges(mask, [4 6 7]);
                else
                    f.charges = f.charges(mask, 3:end);
                end
                f.isdual = f.isdual(3:end);
                if hasmultiplicity(fusionstyle(f)) && ~isempty(f.vertices)
                    f.vertices = f.vertices(mask, 2:end);
                end
                f.rank(1) = f.rank(1) - 2;
                [f, ~, locb] = unique(f);
                c = sparse(mask, locb, vals, L1, length(f));
                return
            end
            
            if i == N
                assert(all(f.coupled == one(f.coupled)), 'tensors:trees:ArgError', ...
                    'Tracing final and first leg allowed iff coupled charge is trivial.');
                if N == 2
                    [c, f] = elementary_trace(f, 1);
                    c = conj(c);
                    return
                end
                
                [c, f] = cycleanticlockwise(f);
                [c_, f] = elementary_trace(f, 1);
                c = c * c_;
                return
            end
            
            mask = find(f.charges(:, 2 * i - 2) == conj(f.charges(:, 2 * i)) & ...
                f.charges(:, 2 * i - 3) == f.charges(:, 2 * i + 1));
            charges1 = f.charges(mask, (2 * i - 3):(2 * i + 1));
            
            b = charges1(:, 2);
            vals = sqrt(qdim(b));
            if hasmultiplicity(fusionstyle(f))
                u = one(b);
                mu = f.vertices(mask, i-1);
                nu = f.vertices(mask, i);
                for j = 1:length(vals)
                    F = Fsymbol(charges1(j, 1), b(j), charges1(j, 4), ...
                        charges1(j, 1), charges1(j, 3), u);
                    vals(j) = vals(j) * F(mu(j), nu(j), 1, 1);
                end
            else
                vals = vals .* Fsymbol(charges1(:, 1), b, charges1(:, 4), ...
                    charges1(:, 1), charges1(:, 3), repmat(one(charges1), length(mask), 1));
            end
            
            if f.isdual(i), vals = vals .* frobeniusschur(b); end
            
            if f.legs == 3
                f.charges = f.charges(mask, [1 end]);
            else
                f.charges = f.charges(mask, [1:(2 * i - 3) (2 * i + 2):end]);
            end
            f.isdual(i:i+1) = [];
            
            if hasmultiplicity(fusionstyle(f))
                f.vertices = f.vertices(mask, [1:i-2 i+1:end]);
            end
            f.rank(1) = f.rank(1) - 2;
            
            [f, ~, locb] = unique(f);
            
            c = sparse(mask, locb, vals, L1, length(f));
        end
        
        function [c, f] = permute(f, p, r)
            % Compute the coefficients that bring a permuted tree into canonical form.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %   tree to permute.
            %
            % p : int
            %   permutation indices.
            %
            % r : int
            %   final number of splitting and fusing legs.
            %
            % Returns
            % -------
            % c : sparse double
            %   matrix of coefficients that transform input to output trees.
            %   `f(i) --> c(i,j) * f(j)`
            %
            % f : :class:`FusionTree`
            %   permuted trees in canonical form.
            
            arguments
                f
                p (1,:)
                r (1,2) = rank(f)
            end
            
            if isempty(f), c = []; return; end
            
            assert(issymmetric(braidingstyle(f)));
            [c, f] = braid(f, p, 1:f.legs, r);
        end
        
        function [c, f] = repartition(f, newrank)
            % Compute the coefficients that bend legs to a desired rank.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %   tree to repartition.
            %
            % newrank : int
            %   new rank of the fusion tree.
            %
            % Returns
            % -------
            % c : sparse double
            %   matrix of coefficients that transform input to output trees.
            %   `f(i) --> c(i,j) * f(j)`
            %
            % f : :class:`FusionTree`
            %   repartitioned trees in canonical form.
            
            arguments
                f
                newrank (1, 2) int16
            end
            oldrank = f.rank;
            assert(sum(newrank) == sum(oldrank), 'tensors:trees:ArgError', ...
                'New rank is incompatible with old rank.');
            
            if all(newrank == oldrank)
                c = speye(length(f));
                return
            end
            
            if newrank(1) > oldrank(1)
                [c, f] = bendleft(f);
                for i = 2:(newrank(1) - oldrank(1))
                    [c_, f] = bendleft(f);
                    c = c * c_;
                end
            else
                [c, f] = bendright(f);
                for i = 2:(oldrank(1) - newrank(1))
                    [c_, f] = bendright(f);
                    c = c * c_;
                end
            end
            [f, p] = sort(f);
            c = c(:, p);
        end
        
        function [c, f] = twist(f, i, inv)
            % Compute the coefficients that twist legs.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %   tree to repartition.
            %
            % i : int or logical
            %   indices of legs to twist.
            %
            % inv : logical
            %   flag to determine inverse twisting.
            %
            % Returns
            % -------
            % c : sparse double
            %   matrix of coefficients that transform input to output trees.
            %   `f(i) --> c(i,j) * f(j)`
            %
            % f : :class:`FusionTree`
            %   twisted trees in canonical form.
            
            arguments
                f
                i
                inv = false
            end
            
            if isempty(f), c = []; return; end
            
            if istwistless(braidingstyle(f)) || isempty(i) || ~any(i)
                c = speye(length(f));
                return
            end
            
            theta = prod(twist(f.uncoupled(:, i)), 2);
            if inv, theta = conj(theta); end
            c = spdiags(theta, 0, length(f), length(f));
        end
    end
    
    
    %% Utility
    methods
        function style = braidingstyle(f)
            style = f.charges.braidingstyle;
        end
        
        function bool = isallowed(f)
            if ~hasmultiplicity(fusionstyle(f))
                if f.rank(1) == 0
                    bool = f.charges(:, 1) == one(f.charges);
                elseif f.rank(1) == 1
                    bool = f.charges(:, 1) == f.charges(:, 2);
                else
                    bool = true(length(f), 1);
                    for i = 1:f.rank(1)-1
                        bool = bool & ...
                            Nsymbol(f.charges(:, 2*i-1), f.charges(:, 2*i), f.charges(:, 2*i+1));
                    end
                end
                if f.rank(2) == 0
                    bool = bool & f.charges(:, end) == one(f.charges);
                elseif f.rank(2) == 1
                    bool = bool & f.charges(:, end) == f.charges(:, end-1);
                else
                    for i = 1:f.rank(2)-1
                        bool = bool & ...
                            Nsymbol(f.charges(:, end-2*i+2), f.charges(:, end-2*i+1), f.charges(:, end-2*i));
                    end
                end
                return
            end
            if f.rank(1) == 0
                bool = f.charges(:, 1) == one(f.charges);
            elseif f.rank(1) == 1
                bool = f.charges(:, 1) == f.charges(:, 2);
            else
                bool = true(length(f), 1);
                for i = 1:f.rank(1)-1
                    bool = bool & ...
                        f.vertices(:, i) <= ...
                        Nsymbol(f.charges(:, 2*i-1), f.charges(:, 2*i), f.charges(:, 2*i+1));
                end
            end
            if f.rank(2) == 0
                bool = bool & f.charges(:, end) == one(f.charges);
            elseif f.rank(2) == 1
                bool = bool & f.charges(:, end) == f.charges(:, end-1);
            else
                for i = 1:f.rank(2)-1
                    bool = bool & ...
                        f.vertices(:, end-i+1) <= ...
                        Nsymbol(f.charges(:, end-2*i+2), f.charges(:, end-2*i+1), f.charges(:, end-2*i));
                end
            end
        end
        
        function f = flip(f)
            f.charges = fliplr(f.charges);
            f.isdual = fliplr(f.isdual);
            f.rank = fliplr(f.rank);
            if hasmultiplicity(fusionstyle(f))
                f.vertices = fliplr(f.vertices);
            end
        end
        
        function style = fusionstyle(f)
            style = fusionstyle(f.charges);
        end
        
        function [lia, locb] = ismember(f1, f2)
            if hasmultiplicity(fusionstyle(f1))
                error('TBA');
            end
            
            [lia, locb] = ismember(f1.charges, f2.charges, 'rows');
        end
        
        function [f, p] = sort(f)
            % sort - Sort the fusion trees into canonical order.
            %   [f, p] = sort(f)
            if isempty(f),  p = []; return; end
            
            if ~isempty(f.vertices) && hasmultiplicity(fusionstyle(f))
                cols = [treeindsequence(f.rank(1)) + 1 ...                                  % center charge
                    (1:treeindsequence(f.rank(2))) + treeindsequence(f.rank(1)) + 1 ...     % fuse charges
                    (f.rank(1):size(f.vertices, 2)) + size(f.charges, 2) ...                % fuse vertices
                    fliplr(1:treeindsequence(f.rank(1))) ...                                % split charges
                    (1:f.rank(1)-1) + size(f.charges, 2)];                                  % split vertices
                
                [p, f.charges, f.vertices] = simulsortrows(f.charges, f.vertices, ...
                    'Col', cols);
            else
                cols = [treeindsequence(f.rank(1)) + 1 ...                                  % center charge
                (1:treeindsequence(f.rank(2))) + treeindsequence(f.rank(1)) + 1 ...         % fuse charges
                fliplr(1:treeindsequence(f.rank(1)))];                                      % split charges
                if nargout > 1
                    [f.charges, p] = sortrows(f.charges, cols);
                else
                    f.charges = sortrows(f.charges, cols);
                end
            end
        end
        
        function bool = issorted(f)
            if ~isempty(f.vertices) && hasmultiplicity(fusionstyle(f))
                [~, p] = sort(f);
                bool = isequal(p, (1:length(f)).');
            else
                bool = issorted(f.coupled);
                if ~bool, return; end
                cols = [treeindsequence(f.rank(1)) + 1 ...                                  % center charge
                    (1:treeindsequence(f.rank(2))) + treeindsequence(f.rank(1)) + 1 ...         % fuse charges
                    fliplr(1:treeindsequence(f.rank(1)))];
                bool = issortedrows(f.charges(:, cols));
            end
        end
    end
    
    %% Setters and Getters
    methods
        function c = get.coupled(f)
            c = f.charges(:, treeindsequence(f.rank(1)) + 1);
        end
        
        function unc = get.uncoupled(f)
            unc = f.charges(:, ...
                [treeindsequence(1:f.rank(1)) ...
                end + 1 - treeindsequence(f.rank(2):-1:1)]);
        end
        
        function f = split(f)
            f.isdual(f.rank(1)+1:end) = [];
            f.charges(:, treeindsequence(f.rank(1)) + 2:end) = [];
            if f.rank(1) >= 2 && hasmultiplicity(fusionstyle(f))
                f.vertices(:, f.rank(1):end) = [];
            else
                f.vertices = uint8.empty(length(f), 0);
            end
            f.rank(2) = 0;
        end
        
        function f = fuse(f)
            f.isdual(1:f.rank(1)) = [];
            f.charges(:, 1:treeindsequence(f.rank(1))) = [];
            if f.rank(2) >= 2 && hasmultiplicity(fusionstyle(f))
                f.vertices(:, 1:f.rank(1)-1) = [];
            else
                f.vertices = uint8.empty(length(f), 0);
            end
            f.rank(1) = 0;
        end
    end
    
    
    %% Utility
    methods
        function bool = isempty(f)
            bool = size(f, 2) == 0;
        end
        
        function l = legs(f)
            l = sum(f.rank);
        end
        
        function l = length(f)
            l = size(f.charges, 1);
        end
        
        function varargout = size(f, varargin)
            if nargin == 1
                if nargout < 2
                    c = f.charges;
                    varargout{1} = [1, size(c, 1)];
                elseif nargout == 2
                    varargout = {1, size(f.charges, 1)};
                else
                    error('wrong input');
                end
                return
            end
            
            sz = size(f);
            if nargout == 1
                varargout{1} = sz([varargin{:}]);
                return
            end
            
            assert(nargout == length(varargin))
            for i = 1:nargout
                varargout{i} = sz(varargin{i});
            end
        end
        
        function n = numArgumentsFromSubscript(~, ~, ~)
            n = 1;
        end
        
        function varargout = subsref(f, s)
            switch s(1).type
                case '()'
                    if length(s(1).subs) == 1
                        i = s(1).subs{1};
                    elseif length(s(1).subs) == 2 && strcmp(s(1).subs{1}, ':')
                        i = s(1).subs{2};
                    else
                        error('Invalid indexing');
                    end
                    
                    if length(s) == 1
                        % f(i)
                        f.charges = f.charges(i, :);
                        if hasmultiplicity(f.fusionstyle)
                            f.vertices = f.vertices(i, :);
                        end
                        varargout{1:numArgumentsFromSubscript(f)} = f;
                        
                    elseif length(s) == 2 && strcmp(s(2).type, '.')
                        % f(i).property
                        switch s(2).subs
                            case {'charges', 'vertices', 'coupled', ...
                                    'uncoupled', 'inner'}
                                prop = f.(s(2).subs);
                                [varargout{1:nargout}] = prop(i, :);
                                
                            case {'isdual', 'rank'}
                                [varargout{1:nargout}] = subsref(f, s(2));
                                
                            otherwise
                                error('tensors:tree:indexing', ...
                                    'Not a valid indexing expression.');
                        end
                        
                    else
                        [varargout{1:nargout}] = builtin('subsref', f, s);
                    end
                    
                otherwise
                    [varargout{1:nargout}] = builtin('subsref', f, s);
            end
        end
        
        function n = numel(f)
            n = size(f.charges, 1);
        end
        
        function bool = ne(f1, f2)
            bool = any(f1.charges ~= f2.charges, 2);
            if hasmultiplicity(fusionstyle(f1))
                bool = bool | any(f1.vertices ~= f2.vertices, 2);
            end
        end
        
        function bool = eq(f1, f2)
            bool = all(f1.charges == f2.charges, 2);
            if hasmultiplicity(fusionstyle(f1))
                bool = bool & all(f1.vertices == f2.vertices, 2);
            end
        end
    end
    
    
    %% Display
    methods (Access = protected)
        function header = getHeader(f)
            if isscalar(f)
                headerFormat = '  (%d, %d) %s %s:\n';
            else
                headerFormat = '  (%d, %d) %s %s array:\n';
            end
            header = sprintf(headerFormat, f.rank(1), f.rank(2), ...
                class(f.charges), matlab.mixin.CustomDisplay.getClassNameForHeader(f));
        end
        function displayScalarObject(f)
            displayNonScalarObject(f);
        end
        function displayNonScalarObject(f)
            header = getHeader(f);
            disp(header);
            fprintf('    isdual:\n    %s\n', num2str(f.isdual));
            fprintf('    charges:\n');
            chargeStr = strjust(pad(string(f.charges)), 'center');
            chargeFormat = ['    ', ...
                repmat('%s ', 1, max(f.rank(1), 2*f.rank(1)-2)), ...
                '| %s |', ...
                repmat(' %s', 1, max(f.rank(2), 2*f.rank(2)-2)) '\n'];
            fprintf(chargeFormat, chargeStr.');
            if hasmultiplicity(fusionstyle(f))
                fprintf('    vertices:\n');
                disp(f.vertices);
            end
        end
    end
    
    
    %% Constructors and converters
    methods
        function A = double(f)
            % Construct a numeric representation of a fusion tree.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %
            % Returns
            % -------
            % A : double
            %   array representation of a fusion tree.
            
            assert(issymmetric(braidingstyle(f)), ...
                'tensors:tree:invalidCharge', ...
                'Tensor representation is only defined for charges with bosonic braiding.');
            
            sz = zeros(1, f.legs);
            uncs = cell(1, f.legs);
            for i = 1:length(sz)
                uncs{i} = unique(f.uncoupled(:, i));
                sz(i) = length(uncs{i});
            end
            A_cell = cell(sz);
            
            inds = arrayfun(@(x) 1:x, sz, 'UniformOutput', false);
            for ind = num2cell(combvec(inds{:}))
                unc = cellfun(@(x, i) x(i), uncs, ind.');
                locb = all(unc == f.uncoupled, 2);
                if any(locb)
                    A_cell{ind{:}} = fusiontensor(subsref(f, substruct('()', {locb})));
                else
                    A_cell{ind{:}} = zeros(qdim(unc));
                end
            end
            
            A = cell2mat(A_cell);
        end
        
        function C = fusiontensor(f)
            % Construct an array representation of a fusion tree.
            %
            % Arguments
            % ---------
            % f : :class:`FusionTree`
            %   an array of fusion trees to convert.
            %
            % Returns
            % -------
            % C : cell
            %   a cell array containing the array representations of f.
            
            if fusionstyle(f) == FusionStyle.Unique
                C = num2cell(ones(size(f)));
                return;
            end
            
            if length(f) > 1
                C = arrayfun(@fusiontensor, f);
                return
            end
            
            
            multi = hasmultiplicity(fusionstyle(f));
            switch f.rank(1)
                case 0
                    C_split = 1;
                case 1
                    if f.isdual(1)
                        C_split = flipper(f.charges(1));
                    else
                        C_split = squeeze(fusiontensor(one(f.charges), ...
                            f.charges(1), f.charges(2)));
                    end
                otherwise
                    C_split = fusiontensor(f.charges(1), f.charges(2), f.charges(3));
                    if multi
                        C_split = C_split(:, :, :, f.vertices(1));
                    end
                    if f.isdual(1)
                        C_split = contract(C_split, [1 -2 -3], ...
                            flipper(f.charges(1)), [-1 1]);
                    end
                    if f.isdual(2)
                        C_split = contract(C_split, [-1 1 -3], ...
                            flipper(f.charges(2)), [-2 1]);
                    end
                    
                    for i = 3:f.rank(1)
                        C_split_ = fusiontensor(f.charges(2*i-3), f.charges(2*i-2), ...
                            f.charges(2*i-1));
                        if multi
                            C_split_ = C_split_(:, :, :, f.vertices(i-1));
                        end
                        if f.isdual(i)
                            C_split_ = contract(C_split_, [-1 1 -3], ...
                                flipper(f.charges(2*i-2)), [-2 1]);
                        end
                        C_split = contract(C_split, [-(1:i-1) 1], ...
                            C_split_, [1 -i -i-1]);
                    end
            end
            switch f.rank(2)
                case 0
                    C_fuse = 1;
                case 1
                    if f.isdual(end)
                        C_fuse = flipper(f.charges(end));
                    else
                        C_fuse = squeeze(fusiontensor(one(f.charges), ...
                            f.charges(end), f.charges(end-1)));
                    end
                otherwise
                    C_fuse = fusiontensor(f.charges(end), f.charges(end-1), ...
                        f.charges(end-2));
                    if multi
                        C_fuse = C_fuse(:, :, :, f.vertices(end));
                    end
                    if f.isdual(end)
                        C_fuse = contract(C_fuse, [1 -2 -3], flipper(f.charges(end)), [-1 1]);
                    end
                    if f.isdual(end-1)
                        C_fuse = contract(C_fuse, [-1 1 -3], flipper(f.charges(end-1)), [-2 1]);
                    end
                    
                    for i = 3:f.rank(2)
                        C_fuse_ = fusiontensor(f.charges(end-2*i+4), f.charges(end-2*i+3), ...
                            f.charges(end-2*i+2));
                        if multi
                            C_fuse_ = C_fuse_(:, :, :, f.vertices(end-i+2));
                        end
                        if f.isdual(end-i+1)
                            C_fuse_ = contract(C_fuse_, [-1 1 -3], ...
                                flipper(f.charges(end-2*i+3)), [-2 1]);
                        end
                        C_fuse = contract(C_fuse, [-(1:i-1) 1], ...
                            C_fuse_, [1 -i -i-1]);
                    end
            end
            if f.rank(1) == 0
                C = permute(conj(C_fuse), max(f.rank(2), 2):-1:1);
            elseif f.rank(2) == 0
                C = C_split;
            else
                C = contract(C_split, [-(1:f.rank(1)) 1], ...
                    conj(C_fuse), [-(f.rank(2):-1:1)-f.rank(1) 1]);
            end
            C = {C};
        end
    end
end
