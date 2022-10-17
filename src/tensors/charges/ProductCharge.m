classdef ProductCharge < AbstractCharge
% Irreducible representations of a direct product of groups.
    
    properties
        charges (1,:) = {}
    end
    
    methods
        function prodcharge = ProductCharge(charges)
            % Construct a product group charge.
            % 
            % Usage
            % -----
            % charges = ProductCharge(charges1, charges2, ...)
            %   creates an (array of) charges that are representations of the direct product
            %   group group1 x group2 x ...
            %
            % Arguments
            % ---------
            % charges1, charges2, ... : AbstractCharge
            %   charges of the separate groups.
            
            arguments (Repeating)
                charges
            end
            
            if nargin == 0
                return
            end
            
            for i = 2:nargin
                assert(all(size(charges{1}) == size(charges{i})), ...
                    'Charges should have equal size.');
            end
            prodcharge.charges = charges;
        end
        
        function a = cat(dim, a, varargin)
            % Concatenate charges.
            
            for i = 1:length(varargin)
                if isempty(varargin{i}), continue; end
                if isempty(a)
                    a = varargin{i};
                    continue;
                end
                for j = 1:length(a.charges)
                    a.charges{j} = cat(dim, a.charges{j}, varargin{i}.charges{j});
                end
            end
        end
        
        function a = horzcat(a, varargin)
            a = cat(2, a, varargin{:});
        end
        
        function a = vertcat(a, varargin)
            a = cat(1, a, varargin{:});
        end
        
        function s = GetMD5_helper(a)
            s = cellfun(@GetMD5_helper, a.charges, 'UniformOutput', false);
        end
        
        function varargout = size(prodcharge, varargin)
            [varargout{1:nargout}] = size(prodcharge.charges{1}, varargin{:});
        end
        
        function a = reshape(a, varargin)
            for i = 1:numel(a.charges)
                a.charges{i} = reshape(a.charges{i}, varargin{:});
            end
        end
        
        function n = ndims(a)
            n = ndims(a.charges{1});
        end
        
        function n = numel(a)
            n = prod(size(a));
        end
        
        function bool = isempty(a)
            bool = isempty(a.charges) || isempty(a.charges{1});
        end
        
        function l = length(a)
            if isempty(a.charges)
                l = 0;
                return
            end
            l = length(a.charges{1});
        end
        
        function [d, N] = prod(a, dim)
            % Total fusion product of charges.
            
            arguments
                a
                dim = find(size(a) ~= 1, 1);
            end
            
            switch fusionstyle(a)
                case FusionStyle.Unique
                    fusedcharges = cell(size(a.charges));
                    for i = 1:length(fusedcharges)
                        fusedcharges{i} = prod(a.charges{i}, dim);
                    end
                    d = ProductCharge(fusedcharges{:});
                    
                    if nargout > 1
                        N = ones(size(d));
                    end
                    
                case FusionStyle.Simple
                    d = subsref(a, substruct('()', {1}));
                    N = 1;
                    for i = 2:length(a)
                        a_i = subsref(a, substruct('()', {i}));
                        d_ = subsref(d, substruct('()', {1})) * ...
                            a_i;
                        if nargout > 1
                            N_(1:length(d_)) = N(1);
                        end
                        
                        for j = 2:length(d)
                            c = subsref(d, substruct('()', {j})) * ...
                                a_i;
                            d_ = [d_ c];
                            if nargout > 1
                                N_(end + (1:length(c))) = N(j);
                            end
                        end
                        
                        if nargout < 2
                            d = unique(d_);
                        else
                            [d, ~, ic] = unique(d_);
                            N = zeros(size(d));
                            for i = 1:length(N)
                                N(i) = sum(N_(ic == i));
                            end
                        end
                    end
                    
                case FusionStyle.Generic
                    error('TBA');
                    d = a(1);
                    N = 1;
                    for i = 2:length(a)
                        d_ = d(1) * a(i);
                        if nargout > 1
                            N_ = N(1) .* ...
                                Nsymbol(repmat(d(1), 1, length(d_)), ...
                                repmat(a(i), 1, length(d_)), d_);
                        end
                        
                        for j = 2:length(d)
                            c = d(j) * a(i);
                            d_(end + (1:length(c))) = c;
                            if nargout > 1
                                N_(end + (1:length(c))) = N(j) .* ...
                                    Nsymbol(repmat(d(j), 1, length(c)), ...
                                    repmat(a(i), 1, length(c)), c);
                            end
                        end
                        
                        if nargout < 2
                            d = unique(d_);
                        else
                            [d, ~, ic] = unique(d_);
                            N = zeros(size(d));
                            for i = 1:length(N)
                                N(i) = sum(N_(ic == i));
                            end
                        end
                    end
            end
        end

        function varargout = subsref(prodcharge, s)
            % Overload indexing.
            %
            % Usage
            % -----
            % charges_slice = charges(i1, i2, ...)
            %   extracts elements out of the charge array.
            %
            % product_slice = charges{i}
            %   separate out the direct product factors.
            %
            % Arguments
            % ---------
            % charges : ProductCharge
            %   array of charges.
            %
            % s : substruct
            %   structure containing indexing data.
            %
            % Returns
            % -------
            % charges_slice : ProductCharge
            %   sliced array of product charges.
            %
            % product_slice : AbstractCharge
            %   array of factor charges.
            
            switch s(1).type
                case '.'
                    [varargout{1:nargout}] = builtin('subsref', prodcharge, s);
                    
                case '()'
                    for i = 1:numel(prodcharge.charges)
                        prodcharge.charges{i} = prodcharge.charges{i}(s(1).subs{:});
                    end
                    if length(s) == 1
                        varargout = {prodcharge};
                        return
                    end
                    
                    [varargout{1:nargout}] = subsref(prodcharge, s(2:end));
                    
                case '{}'
                    assert(length(s) == 1);
                    assert(length(s(1).subs) == 1);
                    assert(nargout == length(s(1).subs{1}));
                    varargout(1:nargout) = prodcharge.charges(s(1).subs{:});
                    
                otherwise
                    error('Undefined behaviour');
            end
        end
        
        function a = subsasgn(a, s, b)
            % Overload indexed assignment.
            %
            % Usage
            % -----
            % a = subsasgn(a, substruct('()', subs), b)
            % a(subs{:}) = b
            %   assign array slices.
            %
            % a = subsasgn(a, substruct('{}', subs), c)
            % a{i} = c
            %   assign to a factor slice.
            %
            % Arguments
            % ---------
            % a : ProductCharge
            %   array of charges to assign to.
            %
            % s : substruct
            %   structure containing indexing data.
            %
            % b : ProductCharge
            %   slice to assign.
            %
            % c : AbstractCharge
            %   factor to assign.
            %
            % Returns
            % -------
            % a : ProductCharge
            %   assigned array
            
            switch s(1).type
                case '()'
                    assert(length(s) == 1, 'Undefined assignment syntax.');
                    for i = 1:length(a.charges)
                        if isempty(b)
                            a.charges{i}(s(1).subs{:}) = [];
                        else
                            a.charges{i}(s(1).subs{:}) = b.charges{i};
                        end
                    end
                    
                case '{}'
                    if length(s) == 1
                        a.charges{s(1).subs{:}} = b;
                    else
                        a.charges{s(1).subs{:}} = subsasgn(a.charges{s(1).subs{:}}, ...
                            s(2:end), b);
                    end
                    
                case '.'
                    a = builtin('subsasgn', a, s, b);
                otherwise
                    error('Undefined behaviour.');
            end
        end
        
        function ind = end(a, k, n)
            ind = builtin('end', a.charges{1}, k, n);
        end
        
        function a = transpose(a)
            if isempty(a), return; end
            for i = 1:numel(a.charges)
                a.charges{i} = a.charges{i}.';
            end
        end
        
        function a = ctranspose(a)
            if isempty(a), return; end
            for i = 1:numel(a.charges)
                a.charges{i} = a.charges{i}.';
            end
        end
        
        function style = braidingstyle(prodcharge)
            style = braidingstyle(prodcharge.charges{1});
            for i = 2:length(prodcharge.charges)
                style = style & braidingstyle(prodcharge.charges{i});
            end
        end
        
        function a = conj(a)
            for i = 1:length(a.charges)
                a.charges{i} = conj(a.charges{i});
            end
        end
        
        function bool = eq(a, b)
            assert(length(a.charges) == length(b.charges));
            bool = eq(a.charges{1}, b.charges{1});
            for i = 2:length(a.charges)
                bool = bool & eq(a.charges{i}, b.charges{i});
            end
        end
        
        function nu = frobeniusschur(a)
            nu = frobeniusschur(a.charges{1}); 
            for i = 2:length(a.charges)
                nu = nu .* frobeniusschur(a.charges{i});
            end
        end
        
        function style = fusionstyle(a)
            style = fusionstyle(a.charges{1});
            for i = 2:length(a.charges)
                style = style & fusionstyle(a.charges{i});
            end
        end
        
        function C = fusiontensor(a, b, c)
            C = fusiontensor(a.charges{1}, b.charges{1}, c.charges{1});
            for i = 2:length(a.charges)
                sz = size(C, 1:4);
                C_ = fusiontensor(a.charges{i}, b.charges{i}, c.charges{i});
                sz_ = size(C_, 1:4);
                C = reshape(contract(C, -(1:2:8), C_, -(2:2:8)), sz .* sz_);
            end
        end
        
        function c = intersect(a, b)
            c = subsref(a, substruct('()', ...
                {mod(find(reshape(a, [], 1) == reshape(b, 1, [])) - 1, length(a)) + 1}));
        end
        
        function F = Fsymbol(a, b, c, d, e, f)
            if hasmultiplicity(fusionstyle(a))
                error('Not implemented yet.');
            end
            
            F = Fsymbol(a.charges{1}, b.charges{1}, c.charges{1}, ...
                d.charges{1}, e.charges{1}, f.charges{1});
            for i = 2:length(a.charges)
                F = F .* Fsymbol(a.charges{i}, b.charges{i}, c.charges{i}, ...
                    d.charges{i}, e.charges{i}, f.charges{i});
            end
        end
        
        function F = Fmatrix(a, b, c, d, e, f)
            % Compute the full recoupling matrix from ``e`` to ``f``.
            %
            % .. todo::
            %   Add proper definition?
            %
            % Usage
            % -----
            % :code:`F = Fmatrix(a, b, c, d, e, f)` computes the matrix between all allowed
            % channels.
            % 
            % Arguments
            % ---------
            % a, b, c : :class:`.AbstractCharge`
            %   charges being fused
            % d : :class:`.AbstractCharge`
            %   total charges
            % e : :class:`.AbstractCharge` (1, \*)
            %   intermediate charges before recoupling
            % f : :class:`.AbstractCharge` (1, \*)
            %   intermediate charge after recoupling
            %
            % Returns
            % -------
            % F : :class:`double` (\*, \*, \*, \*)
            %   recoupling matrix between all allowed channels
            if a.fusionstyle == FusionStyle.Unique
                if nargin < 5, e = a * b; end
                if nargin < 6, f = b * c; end
                F = Fsymbol(a, b, c, d, e, f);
                return
            end
            
            if nargin < 5, e = intersect(a * b, conj(c * conj(d))); end
            if nargin < 6, f = intersect(b * c, conj(conj(d) * a)); end
            
            if hasmultiplicity(a.fusionstyle)
                Fblocks = cell(length(f), length(e));
                for i = 1:length(e)
                    for j = 1:length(f)
                        Fblocks{j, i} = Fsymbol(a, b, c, d, ...
                            subsref(e, substruct('()', {i})), ...
                            subsref(f, substruct('()', {j})));
                        sz = size(Fblocks{j, i}, 1:4);
                        Fblocks{j, i} = reshape(Fblocks{j, i}, ...
                            sz(1) * sz(2), sz(3) * sz(4)).';
                    end
                end
                F = cell2mat(Fblocks);
                return
            end
            
            F = zeros(length(f), length(e));
            for i = 1:length(e)
                for j = 1:length(f)
                    F(j, i) = Fsymbol(a, b, c, d, ...
                        subsref(e, substruct('()', {i})), ...
                        subsref(f, substruct('()', {j})));
                end
            end
        end
        
        function disp(a)
            for i = 1:length(a.charges)
                disp(a.charges{i});
            end
        end
        
        function N = Nsymbol(a, b, c)
            N = Nsymbol(a.charges{1}, b.charges{1}, c.charges{1});
            for i = 2:length(a.charges)
                N = N .* Nsymbol(a.charges{i}, b.charges{i}, c.charges{i});
            end
        end
        
        function c = mtimes(a, b)
            if fusionstyle(a) == FusionStyle.Unique
                charges = cell(size(a.charges));
                for i = 1:numel(charges)
                    charges{i} = mtimes(a.charges{i}, b.charges{i});
                end
                c = ProductCharge(charges{:});
                return
            end
            
            assert(isscalar(a) && isscalar(b))
            charges = cell(size(a.charges));
            ctr = 1;
            for i = 1:length(charges)
                charges{i} = a.charges{i} * b.charges{i};
                n = numel(charges{i});
                charges{i} = reshape(repmat(reshape(charges{i}, 1, []), ctr, 1), 1, []);
                for j = 1:i-1
                    charges{j} = repmat(charges{j}, 1, n);
                end
                ctr = ctr * n;
            end
            c = ProductCharge(charges{:});
        end
        
        function bool = ne(a, b)
            bool = a.charges{1} ~= b.charges{1};
            for i = 2:length(a.charges)
                bool = bool | a.charges{i} ~= b.charges{i};
            end
        end
        
        function d = qdim(a)
            d = qdim(a.charges{1});
            for i = 2:length(a.charges)
                d = d .* qdim(a.charges{i});
            end
        end
        
        function R = Rsymbol(a, b, c, inv)
            if nargin < 4, inv = []; end
            if hasmultiplicity(fusionstyle(a))
                error('TBA');
            end
            
            R = Rsymbol(a.charges{1}, b.charges{1}, c.charges{1}, inv);
            for i = 2:length(a.charges)
                R = R .* Rsymbol(a.charges{i}, b.charges{i}, c.charges{i}, inv);
            end
        end
        
        function theta = twist(a)
            theta = twist(a.charges{1});
            for i = 2:length(a.charges)
                theta = theta .* twist(a.charges{i});
            end
        end
        
        function bool = issorted(a)
            [a, p] = sort(a);
            bool = isequal(p, (1:length(a)).');
        end
        
        function bool = issortedrows(a)
            [~, p] = sortrows(a);
            bool = isequal(p, (1:size(a, 1)).');
        end
        
        function [a, I] = sort(a, varargin)
            [I, a.charges{:}] = simulsort(a.charges{:}, varargin{:});
        end
        
        function [a, I] = sortrows(a, col, direction)
            arguments
                a
                col (1,:) double = 1:size(a, 2)
                direction = 'ascend'
            end
            
            newcol = reshape(col + (((1:length(a.charges)).' - 1) .* size(a, 2)), 1, []);
            [I, a.charges{:}] = simulsortrows(a.charges{:}, ...
                'Col', newcol, ...
                'Direction', direction);
        end
        
        function s = string(a)
            charge_str = cellfun(@string, a.charges, 'UniformOutput', false);
            s = join(cat(ndims(a) + 1, charge_str{:}), ' x ', ndims(a) + 1);
        end
            
        function a = one(a)
            for i = 1:length(a.charges)
                a.charges{i} = one(a.charges{i});
            end
        end
    end
end

