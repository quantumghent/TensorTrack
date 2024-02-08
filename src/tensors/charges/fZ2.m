classdef fZ2 < Z2
    % Fermionic charges, implemented as a graded :math:`Z_2` symmetry.
    
    methods
        function style = braidingstyle(~)
            style = BraidingStyle.Fermionic;
        end
        
        function [charges, vertices] = cumprod(a)
            charges = fZ2(mod(cumsum(a), 2));
            vertices = [];
        end
        
        function c = mtimes(a, b)
            c = fZ2(xor(a, b));
        end
           
        function e = one(~)
            e = fZ2(false);
        end
        
        function p = parity(a)
            p = logical(a);
        end
        
        function [d, N] = prod(a, dim)
            if nargin < 2 || isempty(dim)
                dim = find(size(a) ~= 1, 1);
            end
            d = fZ2(mod(sum(a, dim), 2));
            if nargout > 1
                N = ones(size(d));
            end
        end
        
        function R = Rsymbol(a, b, c, ~)
            R = (-2 .* double(a & b) + 1) .* Nsymbol(a, b, c);
        end
        
        function theta = twist(a)
            theta = -2 * double(a) + 1;
        end
    end
end

