classdef Z2 < AbstractCharge & logical
    % Irreducible representations of :math:`Z_2`.
    %
    % This class implements the trivial or sign representations of :math:`Z_2`, represented
    % using {:code:`false`, :code:`true`} where multiplication is given by
    % :math:`\mathrm{XOR}`, giving the multiplication table:
    %
    % .. list-table::
    %
    %   * - 
    %     - :code:`false`
    %     - :code:`true`
    %   * - :code:`false`
    %     - :code:`false`
    %     - :code:`true`
    %   * - :code:`true`
    %     - :code:`true`
    %     - :code:`false`
    %
    % See Also
    % --------
    % :class:`AbstractCharge`
    
    methods
        function A = Asymbol(a, b, c)
            A = double(Nsymbol(a, b, c));
        end
        
        function B = Bsymbol(a, b, c)
            B = double(Nsymbol(a, b, c));
        end

        function a = conj(a)
        end
        
        function [charges, vertices] = cumprod(a)
            charges = Z2(mod(cumsum(a), 2));
            vertices = [];
        end
        
        function nu = frobeniusschur(a)
            nu = ones(size(a));
        end
        
        function F = Fsymbol(a, b, c, d, e, f)
            F = double(Nsymbol(a, b, e) .* Nsymbol(b, c, f) .* ...
                Nsymbol(e, c, d) .* Nsymbol(a, f, d));
        end
        
        function C = fusiontensor(a, b, c)
            C = double(Nsymbol(a, b, c));
        end
        
        function s = GetMD5_helper(data)
            s = logical(data);
        end
        
        function c = mtimes(a, b)
            c = Z2(xor(a, b));
        end
        
        function N = Nsymbol(a, b, c)
            N = xor(a, b) == c;
        end
        
        function e = one(~)
            e = Z2(false);
        end
        
        function [d, N] = prod(a, dim)
            if nargin < 2 || isempty(dim)
                dim = find(size(a) ~= 1, 1);
            end
            d = Z2(mod(sum(a, dim), 2));
            if nargout > 1
                N = ones(size(d));
            end
        end
        
        function R = Rsymbol(a, b, c, ~)
            R = double(Nsymbol(a, b, c));
        end
        
        function s = string(a)
            s = strings(size(a));
            s(~a) = "+";
            s(a) = "-";
        end
        
        function charge = Z2(varargin)
            if nargin == 0
                labels = [];
            else
                labels = horzcat(varargin{:});
            end
            charge@logical(labels);
        end
        
        function style = braidingstyle(~)
            style = BraidingStyle.Abelian;
        end
        
        function style = fusionstyle(~)
            style = FusionStyle.Unique;
        end
    end
end
