classdef SU2 < AbstractCharge & uint8
    % SU2 - Irreducible representations of SU2.
    %   This class represents the representations of SU2, represented using uint8, such that
    %   the representative is equal to the quantum dimension. The use of uint8 limits the
    %   maximum to 255. The spin label can be recovered with spin = (j - 1) / 2.
    %
    %   See also AbstractCharge
    
    methods
        function style = braidingstyle(~)
            style = BraidingStyle.Bosonic;
        end
        
        function B = Bsymbol(a, b, c)
            B = sqrt(qdim(a) .* qdim(b)) .* ...
                (-1).^(2 * (spin(a) + spin(b))) .* ...
                Wigner6j(spin(a), spin(b), spin(c), spin(b), spin(a), zeros(size(a)));
        end
        
        function a = conj(a)
        end
        
        function varargout = cumprod(a)
            [varargout{1:nargout}] = cumprod@AbstractCharge(a);
        end
        
        function nu = frobeniusschur(a)
            nu = ones(size(a));
            nu(~mod(a, 2)) = -1;
        end
        
        function F = Fsymbol(a, b, c, d, e, f)

            sa = spin(a);
            sb = spin(b);
            sc = spin(c);
            sd = spin(d);
            se = spin(e);
            sf = spin(f);
            
            F = sqrt(qdim(e)) .* sqrt(qdim(f)) .* ...
                (-1).^(sa + sb + sc + sd) .* ...
                Wigner6j(sa, sb, se, sc, sd, sf);
        end
        
        function F = Fmatrix(a, b, c, d, e, f)
            if nargin < 5
                e = intersect(a * b, conj(c * conj(d)));
            end
            if nargin < 6
                f = intersect(b * c, conj(conj(d) * a)); 
            end
            nf = length(f);
            ne = length(e);
            
            a = repmat(a, nf, ne);
            b = repmat(b, nf, ne);
            c = repmat(c, nf, ne);
            d = repmat(d, nf, ne);
            e = repmat(reshape(e, 1, ne), nf, 1);
            f = repmat(reshape(f, nf, 1), 1, ne);
            
            F =  Fsymbol(a, b, c, d, e, f);
        end
        
        function style = fusionstyle(~)
            style = FusionStyle.Simple;
        end
        
        function C = fusiontensor(a, b, c)
            ma = repmat(reshape(spin(a):-1:-spin(a), qdim(a), 1, 1), [1 qdim(b) qdim(c)]);
            mb = repmat(reshape(spin(b):-1:-spin(b), 1, qdim(b), 1), [qdim(a) 1 qdim(c)]);
            mc = repmat(reshape(spin(c):-1:-spin(c), 1, 1, qdim(c)), [qdim(a) qdim(b) 1]);
            
            ja = repmat(spin(a), [qdim(a) qdim(b) qdim(c)]);
            jb = repmat(spin(b), size(ja));
            jc = repmat(spin(c), size(ja));
            
            C = Wigner3j(ja, jb, jc, ma, mb, mc, 4, true);
        end
        
        function c = intersect(a, b)
            c = a(mod(find(a(:) == b(:).') - 1, length(a)) + 1);
        end
        
        function c = mtimes(a, b)
            c = SU2((double(max(a - b, b - a)):2:double(a + b - 2)) + 1);
        end
        
        function N = Nsymbol(a, b, c)
            N = a <= b + c & b <= a + c & c <= a + b & mod(a + b + c, 2) == 1;
        end
        
        function e = one(~)
            e = SU2(ones(1, 1, 'uint8'));
        end
        
        function varargout = prod(varargin)
            [varargout{1:nargout}] = prod@AbstractCharge(varargin{:});
%             if nargout > 0
%                 varargout{1} = SU2(varargout{1});
%             end
        end
        
        function d = qdim(a)
            d = double(a);
        end
        
        function R = Rsymbol(a, b, c, ~)
            R = (-Nsymbol(a, b, c)).^mod(spin(a) + spin(b) - spin(c), 2);
        end
        
        function s = string(a)
            s = strings(size(a));
            halfints = logical(mod(a-1, 2));
            s(halfints) = arrayfun(@(x) sprintf("%d/2", x-1), a(halfints));
            s(~halfints) = arrayfun(@(x) sprintf("%d", (x-1)/2), a(~halfints));
        end
        
        function s = spin(a)
            s = (double(a) - 1) ./ 2;
        end
        
        function charge = SU2(varargin)
            if nargin == 0
                labels = [];
            else
                labels = horzcat(varargin{:});
            end
            charge@uint8(labels);
        end
    end
end
