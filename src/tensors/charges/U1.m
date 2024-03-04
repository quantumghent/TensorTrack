classdef U1 < AbstractCharge & int16
    % Irreducible representations of :math:`\mathrm{U}(1)`.
    %
    % This class represents the representations of :math:`\mathrm{U}(1)`, labeled using
    % integers where fusion is given by addition.
    
    methods
        function charge = U1(varargin)
            if nargin == 0
                labels = [];
            else
                labels = horzcat(varargin{:});
            end
            charge@int16(labels);
        end
        
        function A = Asymbol(a, b, c)
            A = double(Nsymbol(a, b, c));
        end
        
        function style = braidingstyle(~)
            style = BraidingStyle.Abelian;
        end
        
        function B = Bsymbol(a, b, c)
            B = double(Nsymbol(a, b, c));
        end
        
        function a = conj(a)
            a = U1(-a);
        end
        
        function [charges, vertices] = cumprod(a)
            charges = U1(cumsum(a));
            vertices = [];
        end
        
        function nu = frobeniusschur(a)
            nu = ones(size(a));
        end
        
        function F = Fsymbol(a, b, c, d, e, f)
            F = double(Nsymbol(a, b, e) .* Nsymbol(b, c, f) .* ...
                Nsymbol(e, c, d) .* Nsymbol(a, f, d));
        end
        
        function style = fusionstyle(~)
            style = FusionStyle.Unique;
        end
        
        function C = fusiontensor(a, b, c)
            C = double(Nsymbol(a, b, c));
        end
        
        function c = mtimes(a, b)
            c = U1(a + b);
        end
        
        function N = Nsymbol(a, b, c)
            N = a + b == c;
        end
        
        function e = one(~)
            e = U1(zeros(1, 1, 'int16'));
        end
        
        function [d, N] = prod(a, dim)
            if nargin < 2 || isempty(dim)
                dim = find(size(a) ~= 1, 1);
            end
            d = U1(sum(a, dim));
            if nargout > 1
                N = ones(size(d));
            end
        end
        
        function R = Rsymbol(a, b, c, ~)
            R = double(Nsymbol(a, b, c));
        end
        
        function s = string(a)
            s = arrayfun(@(x) sprintf("%+d", x), a);
        end
        
        function y = typecast(x, datatype)
            y = typecast(int16(x), datatype);
        end
        
        function s = GetMD5_helper(data)
            s = int16(data);
        end
    end
end