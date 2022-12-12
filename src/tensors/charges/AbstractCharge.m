classdef (Abstract) AbstractCharge
    % AbstractCharge - Abstract base class for objects in a fusion category.
    
    %% Required categorical data.
    methods
        function style = braidingstyle(a)
            % Trait that describes the braiding style.
            %
            % Arguments
            % --------
            % a : :class:`.AbstractCharge`
            %   input charge
            %
            % Returns
            % -------
            % style : :class:`.BraidingStyle`
            %   braiding style of given charge
            %
            % See Also
            % -------
            % :class:`.BraidingStyle`
            error('AbstractCharge:requiredMethod', ...
                    'Error. \nMethod must be overloaded.')
        end
        
        function abar = conj(a)
            % Compute the dual charge.
            %
            % Arguments
            % --------
            % a : :class:`AbstractCharge`
            %   input charge
            %
            % Returns
            % -------
            % abar : :class:`AbstractCharge`
            %   conjugate charge suche that :code:`one(a)` is an element of :code:`a * abar`
            error('AbstractCharge:requiredMethod', ...
                    'Error. \nMethod must be overloaded.')
        end
        
        function style = fusionstyle(a)
            % Trait that describes the fusion style.
            %
            % Arguments
            % --------
            % a : :class:`.AbstractCharge`
            %   input charge
            %
            % Returns
            % -------
            % style : :class:`.FusionStyle`
            %   fusion style of given charge
            %
            % See Also
            % -------
            % :class:`.FusionStyle`
            error('AbstractCharge:requiredMethod', ...
                    'Error. \nMethod must be overloaded.')
        end
        
        function c = mtimes(a, b)
            % Implement the fusion rules.
            %
            % Arguments
            % ---------
            % a, b : :class:`.AbstractCharge`
            %   charges to be fused
            %
            % Returns
            % -------
            % c : :class:`.AbstractCharge` (1, \*)
            %   the unique elements in the decomposition of the tensor product of ``a`` and
            %   ``b``
            error('AbstractCharge:requiredMethod', ...
                    'Error. \nMethod must be overloaded.')
        end
        
        function F = Fsymbol(a, b, c, d, e, f)
            % Compute the recoupling coefficients.
            %
            % Usage
            % -----
            % :code:`F = Fsymbol(a, b, c, d, e, f)` computes the isomorphism between the
            % following fusion diagrams:
            % 
            % |
            %
            % .. image:: ../img/Fmove.svg
            %    :alt: Fmove
            %    :scale: 6 %
            %    :align: center
            % 
            % |
            %
            % Arguments
            % ---------
            % a, b, c : :class:`.AbstractCharge`
            %   charges being fused
            %
            % d : :class:`.AbstractCharge`
            %   total charges
            %
            % e : :class:`.AbstractCharge`
            %   intermediate charge before recoupling
            %
            % f : :class:`.AbstractCharge`
            %   intermediate charge after recoupling
            %
            % Returns
            % -------
            % F : :class:`double` (\*, \*, \*, \*)
            %   recoupling coefficients
            error('AbstractCharge:requiredMethod', ...
                    'Error. \nMethod must be overloaded.')
        end
        
        function N = Nsymbol(a, b, c)
            % Compute the fusion multiplicities.
            %
            % Arguments
            % ---------
            % a, b : :class:`.AbstractCharge`
            %   charges being fused
            % c : :class:`.AbstractCharge`
            %   resulting charge
            %
            % Returns
            % -------
            % N : :class:`int`
            %   amount of times c appears in the fusion product of a and b
        end
        
        function e = one(a)
            % Compute the trivial charge.
            %
            % Arguments
            % ---------
            % a : :class:`.AbstractCharge`
            %   input charge
            %
            % Returns
            % -------
            % e : :class:`.AbstractCharge`
            %   trivial charge
            error('AbstractCharge:requiredMethod', ...
                    'Error. \nMethod must be overloaded.')
        end
        
    end
    
    %% Optional categorical data
    methods
        function C = fusiontensor(a, b, c)
            % Compute a tensor for the fusion vertex.
            % 
            % Usage
            % -----
            % :code:`C = fusiontensor(a, b, c)` a tensor representation for the following
            % diagram:
            % 
            % |
            %
            % .. image:: ../img/fusiontensor.svg
            %    :alt: fusiontensor
            %    :scale: 6 %
            %    :align: center
            % 
            % |
            %
            % Arguments
            % ---------
            % a, b : :class:`.AbstractCharge`
            %   charges being fused
            % c : :class:`.AbstractCharge`
            %   resulting charge
            %
            % Returns
            % -------
            % C : :class:`double` (\*, \*, \*, \*)
            %   fusion tensor
            error('AbstractCharge:optionalMethod', ...
                    'Error. \nMethod must be overloaded.')
        end
        
        function R = Rsymbol(a, b, c, inv)
            % Compute the braiding coefficients.
            %
            % Usage
            % -----
            % :code:`R = Rsymbol(a, b, c, inv)` computes the isomorphism between the
            % following fusion diagrams:
            % 
            % |
            %
            % .. image:: ../img/Rmove.svg
            %    :alt: Rmove
            %    :scale: 6 %
            %    :align: center
            % 
            % |
            %
            % Arguments
            % ---------
            % a, b : :class:`.AbstractCharge`
            %   charges being fused
            % c : :class:`.AbstractCharge`
            %   resulting charge
            % inv : :class:`logical` , optional
            %   compute inverse braiding coefficient (default false)
            %
            % Returns
            % -------
            % R : :class:`double` (\*, \*)
            %   braiding coefficients
            error('AbstractCharge:optionalMethod', ...
                    'Error. \nMethod must be overloaded.')
        end
    end
    
    %% Categorical data that can be derived.
    methods
        function A = Asymbol(a, b, c)
            % Compute the fusion to splitting coefficient.
            % 
            % .. todo::
            %   Add diagram?
            %
            % Arguments
            % ---------
            % a, b : :class:`.AbstractCharge`
            %   charges being fused
            % c : :class:`.AbstractCharge`
            %   resulting charge
            %
            % Returns
            % -------
            % A : :class:`double`
            %   fusion to splitting coefficient
            A = sqrt(qdim(a) .* qdim(b) ./ qdim(c)) .* ...
                conj(frobeniusschur(a) .* ...
                Fsymbol(conj(a), a, b, b, repmat(one(a), size(a)), c));
        end
        
        function B = Bsymbol(a, b, c)
            % Compute the splitting to fusion coefficient.
            % 
            % .. todo::
            %   Add diagram?
            %
            % Arguments
            % ---------
            % a, b : :class:`.AbstractCharge`
            %   charges being fused
            % c : :class:`.AbstractCharge`
            %   resulting charge
            %
            % Returns
            % -------
            % B : :class:`double`
            %   splitting to fusion coefficient
            B = sqrt(qdim(a) .* qdim(b) ./ qdim(c)) .* ...
                Fsymbol(a, b, conj(b), a, c, repmat(one(a), size(a)));
        end
        
        function B = braidingmatrix(a, b, c, d, e, f, inv)
            % Compute the matrix for general Artin braids.
            % 
            % .. todo::
            %   Complete docstring.
            % 
            if hasmultiplicity(a.fusionstyle)
                R1 = arrayfun(@(x) Rsymbol(d, x, e, inv), c, 'UniformOutput', false);
                R2 = arrayfun(@(x) Rsymbol(d, a, x, ~inv), f, 'UniformOutput', false);
                F = arrayfun(@(x,y) Fsymbol(d, a, b, e, x, y), ...
                    repmat(f(:).', length(c), 1), repmat(c(:), 1, length(f)), ...
                    'UniformOutput', false);
                
                blocks = cell(length(c), length(f));
                for i = 1:length(c)
                    for j = 1:length(f)
                        blocks{i,j} = contract(...
                            R1{i}, [-2 1], ...
                            F{i,j}, [2 -3 -1 1], ...
                            R2{j}, [-4 2]);
                        sz = size(blocks{i,j}, 1:4);
                        blocks{i,j} = reshape(blocks{i,j}, sz(1) * sz(2), sz(3) * sz(4));
                    end
                end
                B = cell2mat(blocks);
                return
            end
            R1 = Rsymbol(c, repmat(d, length(c), 1), repmat(e, length(c), 1), inv);
            R2 = Rsymbol(repmat(d, length(f), 1), repmat(a, length(f), 1), f, ~inv);
            R1 = reshape(R1, [], 1);
            R2 = reshape(R2, 1, []);
            
            F = Fmatrix(d, a, b, e, f, c);
            B = R1 .* F .* R2;
        end
        
        function F = flipper(a)
            % Create a matrix-representation of an arrowflip.
            % 
            % .. todo::
            %   Add diagram or definition?
            %
            % Arguments
            % ---------
            % a : :class:`.AbstractCharge`
            %
            % Returns
            % -------
            % F : :class:`double` (\*, \*)
            %   matrix-representation of an arrowflip
            F = conj(sqrt(qdim(a)) .* fusiontensor(conj(a), a, one(a)));
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
                        Fblocks{j, i} = Fsymbol(a, b, c, d, e(i), f(j));
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
                    F(j, i) = Fsymbol(a, b, c, d, e(i), f(j));
                end
            end
        end
        
        function nu = frobeniusschur(a)
            % Compute the Frobenius-Schur indicator.
            % 
            % Arguments
            % ---------
            % a : :class:`.AbstractCharge`
            %
            % Returns
            % -------
            % nu : :class:`double`
            %   Frobenius-Schur indicator of ``a``
            nu = sign(Fsymbol(a, conj(a), a, a, one(a), one(a)));
        end
        
        function d = qdim(a)
            % Compute the quantum dimension of a charge.
            % 
            % Arguments
            % ---------
            % a : :class:`.AbstractCharge`
            %
            % Returns
            % -------
            % d : :class:`double`
            %   quantum dimension of ``a``
            if fusionstyle(a) == FusionStyle.Unique
                d = ones(size(a));
            else
                d = zeros(size(a));
                for i = 1:numel(d)
                    F = Fsymbol(a(i), conj(a(i)), a(i), a(i), one(a(i)), one(a(i)));
                    d(i) = abs(1 / F(1));
                end
            end
        end
        
        function theta = twist(a)
            % Compute the coefficient obtained by twisting a charge.
            %
            % .. todo::
            %   Add diagram/definition?
            %
            % Arguments
            % ---------
            % a : :class:`.AbstractCharge`
            %
            % Returns
            % -------
            % theta : :class:`double`
            %   twist coefficient of ``a``
            if braidingstyle(a) == BraidingStyle.Bosonic
                theta = ones(size(a));
                return
            end
            if numel(a) > 1, theta = arrayfun(@twist, a); return; end
            
            theta = 0;
            for b = a * a
                theta = theta + qdim(b) / qdim(a) * trace(Rsymbol(a, a, b));
            end
        end
        
        function p = parity(a)
            if istwistless(braidingstyle(a))
                p = false(size(a));
                return
            end
            
            error('Non-bosonic charges should implement a parity.')
        end
    end
    
    %% Utility functions
    methods
        function [d, N] = prod(a, dim)
            % Total fusion product of charges.
            %
            % .. todo::
            %   Complete docstring
            %
            % Arguments
            % ---------
            % a, b, c, ... : :class:`.AbstractCharge`
            %
            % Returns
            % -------
            % c : :class:`.AbstractCharge`
            %   total fusion product determined by subsequently multiplying input charges
            arguments
                a
                dim = find(size(a) ~= 1, 1)
            end
            
            switch fusionstyle(a)
                case FusionStyle.Unique
                    if dim == 1
                        d = a(1, :);
                        for i = 2:size(a, 1)
                            d = d * a(i, :);
                        end
                    else
                        d = a(:, 1);
                        for i = 2:size(a, 2)
                            d = d * a(:, i);
                        end
                    end
                    if nargout > 1
                        N = ones(size(d));
                    end
                    
                case FusionStyle.Simple
                    d = a(1);
                    N = 1;
                    for i = 2:length(a)
                        d_ = d(1) * a(i);
                        if nargout > 1
                            N_(1:length(d_)) = N(1);
                        end
                        
                        for j = 2:length(d)
                            c = d(j) * a(i);
                            d_(end + (1:length(c))) = c;
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
        
        function [charges, vertices] = cumprod(a)
            % Cumulative fusion product of elements.
            %
            % .. todo::
            %   Complete docstring
            %
            % Usage
            % -----
            % :code:`Y = cumprod(X)` computes the cumulative fusion product along the
            % columns of X.
            %
            % For ``FusionStyle.Unique``, ``Y`` has the same size as ``X``, which can be
            % an arbitrarily-sized matrix.
            % For other fusion styles, ``X`` should be a vector, and ``Y`` is a matrix
            % containing all different combinations.
            %
            % Arguments
            % ---------
            % a: :class:`.AbstractCharge` (\*, \*)
            %
            % Returns
            % -------
            % charges : :class:`.AbstractCharge`
            %   explanation
            % vertices : type
            %   explanation
            
            if a.fusionstyle == FusionStyle.Unique
                charges = a;
                for i = 2:size(charges, 1)
%                     charges(i, :) = charges(i-1, :) * charges(i, :);
                    charges = subsasgn(charges, substruct('()', {i, ':'}), ...
                        subsref(charges, substruct('()', {i-1, ':'})) * ...
                        subsref(charges, substruct('()', {i, ':'})));
                end
                vertices = [];
                return
            end
            
            % TODO avoid loop allocations, ideas: vertcat cell + unique?
            % Not so simple for non-unique fusion style.
            if size(a, 2) > 1
                charges = [];
                vertices = [];
                for i = 1:size(a, 2)
                    [chargepart, vertexpart] = cumprod(...
                        subsref(a, substruct('()', {':', i})));
                    charges = horzcat(charges, chargepart);
                    vertices = horzcat(vertices, vertexpart);
                end
                return
            end
            
            if a.fusionstyle == FusionStyle.Simple
                if length(a) == 1
                    charges = a;
                    vertices = [];
                    return
                end
                
                if length(a) == 2
                    f = subsref(a, substruct('()', {1})) * subsref(a, substruct('()', {2}));
                    charges = [repmat(subsref(a, substruct('()', {1})), 1, length(f)); f];
                    vertices = [];
                    return
                end
                
                part = cumprod(subsref(a, substruct('()', {1:length(a)-1})));
                charges = [];
                for i = 1:size(part, 2)
                    f = subsref(part, substruct('()', {size(part, 1), i})) * ...
                        subsref(a, substruct('()', {length(a)}));
                    charges = [charges ...
                        [repmat(subsref(part, substruct('()', {1:size(part, 1), i})), ...
                        1, length(f)); f]];
                end
                vertices = [];
                return
            end
            
            if length(a) == 1
                charges = a;
                vertices = [];
                return
            end
            
            if length(a) == 2
                charges = [];
                vertices = [];
                for f = subsref(a, substruct('()', {1})) * subsref(a, substruct('()', {2}))
                    N = Nsymbol(subsref(a, substruct('()', {1})), ...
                        subsref(a, substruct('()', {2})), f);
                    charges = [charges repmat([subsref(a, substruct('()', {1})); f], 1, N)]; 
                    vertices = [vertices 1:N];
                end
                return
            end
            
            [chargepart, vertexpart] = cumprod(a(1:end-1));
            charges = [];
            vertices = [];
            for i = 1:size(chargepart, 2)
                for f = subsref(chargepart, substruct('()', {size(chargepart, 1), i})) * ...
                        subsref(a, substruct('()', {length(a)}))
                    N = Nsymbol(subsref(chargepart, substruct('()', {size(chargepart, 1), i})), ...
                        subsref(a, substruct('()', {length(a)})), f);
                    charges = [charges, ...
                        repmat([subsref(chargepart, substruct('()', {':', i})); f], 1, N)];
                    vertices = [vertices [repmat(vertexpart(:, i), 1, N); 1:N]];
                end
            end
        end
        
        function Y = combvec(X)
            % Create all combinations of vectors.
            %
            % .. todo::
            %   Complete docstring
            %
            % Usage
            % -----
            % :code:`Y = combvec(X1, X2, ...)` takes any number of inputs ``X``, where
            % each ``Xi`` has ``Ni`` columns, and returns a matrix of ``prod(N)`` column
            % vectors, where the columns consist of all combinations found by
            % combining one column vector from each ``Xi``.
            %
            % Arguments
            % ---------
            % X: :class:`.AbstractCharge`, repeating
            %   input charges
            %
            % Returns
            % -------
            % Y : :class:`.AbstractCharge`
            %   explanation
            %
            % See Also
            % --------
            % :func:`combvec`
            arguments (Repeating)
                X
            end
            
            Y = X{1};
            for i = 2:length(X)
                Y = [repmat(Y, 1, size(X{i}, 2))
                    reshape(repmat(X{i}, size(Y, 2), 1), size(X{i}, 1), [])];
            end
        end
        
        function [lia, locb] = ismember_sorted(a, b)
            if isempty(a) || isempty(b)
                lia = false(size(a));
                locb = zeros(size(a));
                return
            end
            
            if isscalar(a) || isscalar(b)
                if any(size(a) == numel(a)) && any(size(b) == numel(b))
                    lia = a == b;
                else
                    lia = reshape(a == b, [], 1);
                end
                
                if ~any(lia)
                    lia = false(size(a));
                    locb = zeros(size(a));
                    return
                end
                if ~isscalar(b)
                    locb = find(lia);
                    locb = locb(1);
                    lia = any(lia);
                else
                    locb = double(lia);
                end
                return
            end
            
            [sortab, indsortab] = sort([reshape(a, [], 1); reshape(b, [], 1)]);
            d = subsref(sortab, substruct('()', {1:length(sortab)-1})) == ...
                subsref(sortab, substruct('()', {2:length(sortab)}));
            ndx1 = indsortab(d);
            
            if nargout <= 1
                lia = ismember(1:length(a), ndx1);
            else
                szuA = length(a);
                d = find(d);
                [lia, locb] = ismember(1:szuA, ndx1);
                newd = d(locb(lia));
                where = indsortab(newd+1) - szuA;
                locb(lia) = where;
            end
            lia = reshape(lia, size(a));
            if nargout > 1
                locb = reshape(locb, size(a));
            end
        end
    end
end
