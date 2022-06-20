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
                error('Not implemented.');
            end
            if inv
                R1 = conj(Rsymbol(repmat(d, length(c), 1), c, repmat(e, length(c), 1)));
                R2 = Rsymbol(repmat(d, length(f), 1), repmat(a, length(f), 1), f);
            else
                R1 = Rsymbol(c, repmat(d, length(c), 1), repmat(e, length(c), 1));
                R2 = conj(Rsymbol(repmat(a, length(f), 1), repmat(d, length(f), 1), f));
            end
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
                        Fblocks{j, i} = reshape(Fblocks{j, i}, ...
                            size(Fblocks{j, i}, 1) * size(Fblocks{j, i}, 2), ...
                            size(Fblocks{j, i}, 3) * size(Fblocks{j, i}, 4));
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
                F = Fsymbol(a, conj(a), a, a, one(a), one(a));
                d = abs(1 / F(1));
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
                    
                    error('TBA');
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
                    [chargepart, vertexpart] = cumprod(a(:, i));
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
                    f = a(1) * a(2);
                    charges = [repmat(a(1), 1, length(f)); f];
                    vertices = [];
                    return
                end
                
                part = cumprod(a(1:end-1));
                charges = [];
                for i = 1:size(part, 2)
                    f = part(end, i) * a(end);
                    charges = [charges [repmat(part(:, i), 1, length(f)); f]];
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
                for f = a(1) * a(2)
                    N = Nsymbol(a(1), a(2), f);
                    charges = vertcat(charges, ...
                        repmat([a(1) f], N, 1));
                    vertices = vertcat(vertices, (1:N)');
                end
                return
            end
            
            [chargepart, vertexpart] = cumprod(a(1:end-1));
            charges = [];
            vertices = [];
            for i = 1:size(chargepart, 1)
                for f = chargepart(i, end) * a(end)
                    N = Nsymbol(chargepart(i, end), a(end), f);
                    charges = vertcat(charges, ...
                        repmat([chargepart(i, :) f], N, 1));
                    vertices = vertcat(vertices, ...
                        [repmat(vertexpart(i,:), N, 1) (1:N)']);
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
    end
    
    methods
%         bool = issortedrows(A)
%         bools = ne(A, B)
%         [B, I] = sort(A, varargin)
    end
end
