classdef fSU2 < SU2
    % Fermionic :math:`\mathrm{SU}(2)` charges, used to represent fermionspin.
    %
    % This is equivalent to representations of :math:`\mathrm{SU}(2) \otimes fZ_2`, but
    % restricted to only allow for the combinations of integer with trivial and halfinteger
    % with fermion charges.
    
    methods
        function style = braidingstyle(~)
            style = BraidingStyle.Fermionic;
        end
        
        function c = mtimes(a, b)
            c = fSU2(mtimes@SU2(a, b));
        end
        
        function e = one(~)
            e = fSU2(ones(1, 1, 'uint8'));
        end
        
        function p = parity(a)
            p = ~mod(a, 2);
        end
        
        function R = Rsymbol(a, b, c, inv)
            if nargin < 4, inv = false; end
            R = (-2 .* double(parity(a) & parity(b)) + 1) .* Rsymbol@SU2(a, b, c, inv);
        end
        
        function theta = twist(a)
            theta = -2 * double(parity(a)) + 1;
        end
    end
end

