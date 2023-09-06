classdef Expand
    % Bond expansion algorithm for uniform matrix product states.
    
    %% Options
    properties
        bondsmethod {mustBeMember(bondsmethod, ...
            {'off', 'factor', 'explicit', 'extrapolate', 'twosite', 'twosite_simple'})} = 'factor'
        
        chargesmethod {mustBeMember(chargesmethod, ...
            {'off', 'fusionproduct', 'twosite', 'twosite_simple'})} = 'off'
        
        % general expansion options
        schmidtcut = 1e-5
        notrunc = false

        noisefactor = 1e-3
        which = 'largestabs' % needed for twosite update; needs to be passed from fixedpoint algorithm!
        
        % bond expansion options
        minbond = 1
        maxbond = 1e9
        tolbond = 0.2
        bondfactor = 1.2
        cutfactor = 1
        explicitbonds = []
        
        % charge expansion options
        mincharges = 2
        
        % optional finalization
        finalize = []
    end
        
    
    %%
    methods
        function v = Expand(kwargs)
            arguments
                kwargs.?Expand
            end
            
            fields = fieldnames(kwargs);
            if ~isempty(fields)
                for field = fields.'
                    v.(field{1}) = kwargs.(field{1});
                end
            end
        end
        
        function [mps2, flag] = changebonds(alg, mpo, mps1)
            % Change sectors and bond dimensions of mps virtual spaces.
            
            % canonicalize before starting
            for d = 1:depth(mps1)
                mps1(d) = canonicalize(mps1(d));
            end

            [new_spaces, flag] = determine_new_spaces(alg, mpo, mps1);
            new_spaces = fix_injectivity(alg, mps1, new_spaces);
            mps2 = expand_mps(alg, mps1, new_spaces);
            
            % finalize
            if ~isempty(alg.finalize)
                mps2 = approximate(alg.finalize, mpo, mps1, mps2);
            end
        end
        
        function [new_spaces, flag] = determine_new_spaces(alg, mpo, mps)
            % Determine which charges and bond dimensions to add/remvove in mps virtual
            % spaces.
            
            for d = depth(mps):-1:1
                for w = period(mps):-1:1
                    if startsWith(alg.bondsmethod, 'twosite') || startsWith(alg.chargesmethod, 'twosite')
                        % twosite expansion takes care of both bonds and charges at the same time
                        [addbond, addcharge] = expand_twosite(alg, mpo, mps, d, w);
                    else
                        % heuristic charge and bond expansion based on current spectrum
                        [S, C] = schmidt_values(mps(d), w);
                        addbond = expand_bonds(alg, S);
                        addcharge = expand_charges(alg, mps, d, w, S, C);
                    end
                    if alg.notrunc
                        addbond = max(0, addbond);
                    end
                    flag = addcharge.flag ~= 0 || nnz(addbond) > 0;
                    
                    old_space = rightvspace(mps(d).AR{w});
                    
                    new_space = struct;
                    new_space.charges = charges(old_space);
                    new_space.degeneracies = degeneracies(old_space) + addbond;
                    
                    if addcharge.flag && isdual(old_space)
                        addcharge.charges = conj(addcharge.charges);
                    end
                    
                    if addcharge.flag == +1          % add charges
                        new_space.charges = [new_space.charges, addcharge.charges];
                        new_space.degeneracies = [new_space.degeneracies, addcharge.bonds];
                    elseif addcharge.flag == -1        % remove charges
                        keep = ~ismember(new_space.charges, addcharge.charges);
                        new_space.charges = new_space.charges(keep);
                        new_space.degeneracies = new_space.degeneracies(keep);
                    end
                    
                    new_spaces(d, w) = old_space.new(new_space, isdual(old_space));
                end
            end
        end

        function spaces = fix_injectivity(alg, mps, spaces)
            % attempt to ensure injectivity
            for d = 1:depth(mps)
                flag = false(1, period(mps));
                while ~all(flag)
                    for w = period(mps):-1:1
                        ww = next(w, period(mps));
                        www = prev(w, period(mps));
                        rhs1 = prod([spaces(d, www) pspace(mps, w)'], ...
                            isdual(spaces(d, w)));
                        rhs2 = prod([spaces(d, ww) pspace(mps, ww)], ...
                            isdual(spaces(d, w)));
                        tmp = spaces(d, w);
                        spaces(d, w) = infimum(spaces(d, w), rhs1);
                        spaces(d, w) = infimum(spaces(d, w), rhs2);
                        assert(dims(spaces(d, w)) > 0, 'no fusion channels left');
                        flag(w) = tmp == spaces(d, w);
                    end
                end
            end
        end
        
        function [addbond, addcharge] = expand_twosite(alg, mpo, mps, d, w)
            % Does a twosite update and adds/removes charges that are above/under the cut
            % in the new spectrum.
            % For larger unit cells, cut is chosen to the right of site w, so
            % couple sites w and w+1.
            
            [svals, charges] = schmidt_values(mps(d), w);
            dd = prev(d, depth(mps));
            ww = next(w, period(mps));
            
            % perform twosite update, take SVD and normalize
            [GL, GR] = environments(Vumps, mpo, mps); % should be able to input this...
            H_AC2 = AC2_hamiltonian(mpo, mps, GL, GR, w);
            AC2 = computeAC2(mps, 1, w);
            [AC2, ~] = eigsolve(H_AC2{1}, AC2.var, 1, alg.which);
            if strcmp(alg.bondsmethod, 'twosite')
                %twosite fixed point
                [AC2, ~] = eigsolve(H_AC2{1}, AC2, 1, alg);
            else
                %single application
                AC2 = apply(H_AC2{1}, AC2);
            end
            [~, C2, ~] = tsvd(AC2.var, ...
                1:mps(dd).AC{w}.plegs+1,  ...
                (1:mps(dd).AR{ww}.plegs+1) + 1 + mps(dd).AC{w}.plegs, ...
                'TruncBelow', alg.schmidtcut, ...
                'TruncDim', alg.maxbond);
            [svals2, charges2] = matrixblocks(C2);
            
            added_charges = setdiff(charges2, charges);
            removed_charges = setdiff(charges, charges2);
            if ~isempty(removed_charges) && ~alg.notrunc
                % need to remove sectors
                addcharge.flag = -1;
                addcharge.charges = removed_charges;
            elseif ~isempty(added_charges)
                % need to add sectors
                addcharge.flag = +1;
                addcharge.charges = added_charges;
                addcharge.bonds = ones(size(added_charges));
                for ii = 1:size(added_charges)
                    addcharge.bonds(ii) = length(svals2{added_charges(ii) == charges2});
                end
            else
                % do nothing
                addcharge.flag = 0;
            end
            
            % determine new bonds
            addbond = zeros(size(svals));
            for ii = 1:length(svals)
                if ismember(charges(ii), charges2)
                    chi_new = length(svals2{charges(ii) == charges2});
                    chi_new = between(alg.minbond, chi_new, alg.maxbond);
                    addbond(ii) = chi_new - length(svals{ii});
                else % charge was removed
                    addbond(ii) = 0;
                end
            end
            
            % enforce alg.notrunc option
            if alg.notrunc
                if addcharge.flag == -1
                    addcharge = struct('flag', 0, 'bonds', [], 'charges', []);
                end
                addbond = max(addbond, 0);
            end

        end

        function addcharge = expand_charges(alg, mps, d, w, svals, charges)
            % Determine whether to add or substract sectors.
            %   addcharge.flag = {+1 [0] -1} for addition/subtraction of sectors.
            %   addcharge.charges = [new_charges]
            %   addcharge.bonds = [new_bonds]
            
            addcharge = struct('flag', false, 'bonds', [], 'charges', []);
            
            % determine cut if alg.maxbond is reached:
            cut = alg.schmidtcut;
            for ii = 1:length(svals)
                if length(svals{ii}) > alg.maxbond
                    cut = max(cut, svals{ii}(alg.maxbond));
                end
            end
            
            highest_schmidt = cellfun(@(x) x(1), svals);
            
            switch alg.chargesmethod
                
                case {'off', []}
                    
                    addcharge.flag = 0;
                    
                case 'fusionproduct'
                    % Standard case, adds all remaining charges in fusion product of virtual
                    % and physical space if none of current charges has a highest schmidt
                    % value above the specified cut.
                    % Removes all but one charge whose highest schmidt value falls below
                    % the specified cut.
                    
                    if all(highest_schmidt > cut)
                        % need to add sectors
                        addcharge.flag = +1;
                        ww = prev(w, period(mps));
                        new_sector = leftvspace(mps(d), ww) * prod(pspace(mps, ww));
                        new_charges = new_sector.dimensions.charges;
                        addcharge.charges = setdiff(new_charges, charges);
                        addcharge.bonds = ones(size(addcharge.charges));
                    elseif length(svals) > alg.mincharges && sum(highest_schmidt < cut) >= 2
                        % need to remove sectors
                        addcharge.flag = -1;
                        below = find(highest_schmidt < cut);
                        [~, keep] = max(highest_schmidt(below));
                        below(keep) = [];
                        addcharge.charges = charges(below);
                    else
                        % do nothing
                        addcharge.flag = 0;
                    end
                                        
            end
            
            % enforce alg.notrunc option
            if addcharge.flag == -1 && alg.notrunc
                addcharge = struct('flag', 0, 'bonds', [], 'charges', []);
            end
            
        end

        function addbond = expand_bonds(alg, svals)
            % Determine whether to add or substract bond dimension. Returns an  1 x n array
            % :code:`addbond` of integers, where :code:`addbond(i)` is the bond dimension to
            % be added to/subtracted from sector :code:`i`.
            
            %new bonds are explicitly provided
            if strcmp(alg.bondsmethod, 'explicit')
                addbond = alg.explicitbonds;
                return
            end

            % recursive
            if iscell(svals)
                addbond = zeros(size(svals));
                for ii = 1:length(svals)
                    addbond(ii) = expand_bonds(alg, svals{ii});
                end
                return
            end
            
            cut = alg.schmidtcut;
            if length(svals) > alg.maxbond
                cut = max(cut, svals(alg.maxbond));
            end
            chi = length(svals);
            chi_tol = ceil(alg.tolbond * chi);
            
            switch alg.bondsmethod
                case 'factor'
                    % new bond is factor * old bond
                    
                    if svals(end) > cut
                        % need to add bond
                        chi_new = ceil(alg.bondfactor * chi);
                    elseif chi > alg.minbond && svals(max(1, chi-chi_tol)) < cut && ~alg.notrunc
                        % need to subtract bond
                        chi_new = min(nnz(svals > cut) + floor(chi_tol / 2), chi);
                    else
                        % do nothing
                        chi_new = chi;
                    end

                case 'extrapolate'
                    % do an extrapolation in order to determine new bond
                    
                    if svals(end) > cut
                        % need to add bond
                        chi_new = extrapolate_schmidt(svals, cut) + floor(chi_tol / 2);
                    elseif chi > alg.minbond && svals(max(1, chi-chi_tol)) < cut && ~alg.notrunc
                        % need to subtract bond
                        chi_new = min(nnz(svals > cut) + 1 + ceil(chi_tol / 2), chi);
                    else
                        % do nothing
                        chi_new = chi;
                    end
                    
                otherwise
                    % no method specified
                    chi_new = chi;
                    
            end
            
            % enforce min/max bonds.
            chi_new = between(alg.minbond, chi_new, alg.maxbond);
            addbond = chi_new - chi;
            
            % enforce alg.notrunc option
            if addbond < 0 && alg.notrunc
                addbond = 0;
            end
            
            % auxiliary function
            function chi = extrapolate_schmidt(svals, cut)
                if length(svals) < 5
                    chi = max(alg.minbond, ceil(alg.bondfactor * length(svals)));
                    return
                end
                y = reshape(log10(svals), [length(svals), 1]);
                x = reshape(1:length(y), [length(svals), 1]);
                range = floor(length(y) * 1/3):length(y);
                ab = [x(range), 1 + 0 * x(range)] \ y(range);
                chi = ceil((log10(cut / alg.cutfactor) -ab(2)) / ab(1));
            end
        end
        
        function mps = expand_mps(alg, mps, new_spaces)
            % Update mps virtual spaces according to given charge and bond expansion.
            
            for d  = 1:depth(mps)
                % expand tensors and convert to UniformMps
                AR = mps(d).AR;
                for w = 1:period(mps)
                    AR{w} = expand(AR{w}, new_spaces(d, prev(w, period(mps))), ...
                        new_spaces(d, w), alg.noisefactor);
                end
                
                mps(d) = UniformMps(AR);
            end
        end
        
    end
end
