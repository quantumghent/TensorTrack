classdef fU1 < U1
    % Fermionic U1 charges.
    %   This is equivalent to representations of U1 x fZ2, but restricted to only allow
    %   for the combinations of even with trivial and odd with fermion charges.
    
    methods
        function style = braidingstyle(~)
            style = BraidingStyle.Fermionic;
        end
        
        function a = conj(a)
            a = fU1(conj@U1(a));
        end
        
        function varargout = cumprod(a)
            [varargout{1:nargout}] = cumprod@U1(a);
            varargout{1} = fU1(varargout{1});
        end
        
        function c = mtimes(a, b)
            c = fU1(mtimes@U1(a, b));
        end
        
        function e = one(a)
            e = fU1(one@U1(a));
        end
        
        function p = parity(a)
            p = logical(mod(a, 2));
        end
        
        function varargout = prod(varargin)
            [varargout{1:nargout}] = prod@U1(varargin{:});
            varargout{1} = fU1(varargout{1});
        end
        
        function R = Rsymbol(a, b, c, inv)
            if nargin < 4, inv = false; end
            R = (-2 .* double(parity(a) & parity(b)) + 1) .* Rsymbol@U1(a, b, c, inv);
        end
        
        function theta = twist(a)
            theta = -2 * double(parity(a)) + 1;
        end
    end
end

