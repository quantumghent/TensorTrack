classdef A4 < AbstractCharge & uint8
    % Irreducible representations of the alternating group of order 4. This is the group of
    % even permutations on four elements, or of orientation-preserving symmetries of a
    % regular tetrahedron. This is a non-abelian group with multiplicities.
    %
    % The representations are labeled by the integers 1-4, with degrees [1 1 1 3].
    
    methods
        function charge = A4(varargin)
            if nargin == 0
                labels = [];
            else
                labels = horzcat(varargin{:});
            end
            charge@uint8(labels);
        end
        
        function B = Bsymbol(a, b, c)
            B = eye(Nsymbol(a, b, c));
        end
        
        function style = braidingstyle(~)
            style = BraidingStyle.Bosonic;
        end
        
        function b = conj(a)
            b = a;
            b(a == 2) = A4(3);
            b(a == 3) = A4(2);
        end
        
        function varargout = cumprod(a)
            [varargout{1:nargout}] = cumprod@AbstractCharge(a);
        end
        
        function nu = frobeniusschur(a)
            nu = ones(size(a));
        end
        
        function F = Fsymbol(a, b, c, d, e, f)
            persistent Fcache
            if isempty(Fcache)
                load('A4_data.mat', 'Fcache');
            end
            
            if ~(Nsymbol(a, b, e) && Nsymbol(e, c, d) && ...
                    Nsymbol(b, c, f) && Nsymbol(a, f, d))
                assert(isempty(Fcache{a,b,c,d,e,f}));
                F = zeros(Nsymbol(a,b,e), Nsymbol(e,c,d), Nsymbol(b,c,f), Nsymbol(a, f, d));
                return
            end
            
            F = Fcache{a, b, c, d, e, f};
        end
        
        function style = fusionstyle(~)
            style = FusionStyle.Generic;
        end
        
        function C = fusiontensor(a, b, c)
            persistent Ccache
            if isempty(Ccache)
                load('A4_data.mat', 'Ccache');
            end
            
            C = Ccache{a, b, c};
        end
        
        function c = mtimes(a, b)
            if a > b
                c = b * a;
                return
            end
            
            if a == 1
                c = b;
            elseif a == 4
                c = A4(1:4);
            elseif b == 4
                c = b;
            elseif a == 3
                c = A4(2);
            else % a == 2
                if b == 2
                    c = A4(3);
                else % b == 3
                    c = A4(1);
                end
            end
        end
        
        function N = Nsymbol(a, b, c)
            persistent Ncache
            if isempty(Ncache)
                load('A4_data.mat', 'Ncache');
            end
            
            linearinds = double(a(:) + 4 * (b(:)-1) + 16 * (c(:)-1));
            N = reshape(Ncache(linearinds), size(a));
        end
        
        function e = one(~)
            e = A4(1);
        end
        
        function varargout = prod(varargin)
            [varargout{1:nargout}] = prod@AbstractCharge(varargin{:});
            if nargout > 0
                varargout{1} = A4(varargout{1});
            end
        end
        
        function d = qdim(a)
            d = ones(size(a));
            d(a == 4) = 3;
        end
        
        function R = Rsymbol(a, b, c, ~)
            persistent Rcache
            if isempty(Rcache)
                load('A4_data.mat', 'Rcache');
            end
            
            R = Rcache{a, b, c};
        end
        
        function s = string(a)
            s = compose("%d", a);
        end
    end
end

