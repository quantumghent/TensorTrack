classdef ZN < AbstractCharge & uint8
    % ZN - Irreducible representations of ZN.
    %   This class represents representations of the cyclic group of order N.
    %   
    %   See also AbstractCharge, Z2, Z3, Z4
    
    properties
        N (1,1) uint8 = uint8(1)
    end
    
    methods
        function A = Asymbol(a, b, c)
            A = double(Nsymbol(a, b, c));
        end
        
        function style = braidingstyle(~)
            style = BraidingStyle.Abelian;
        end
        
        function B = Bsymbol(a, b, c)
            B = double(Nsymbol(a, b, c));
        end
        
        function c = cat(dim, varargin)
            if nargin == 2
                c = varargin{1};
                return
            end
            if nargin > 3
                c = cat(dim, cat(dim, varargin{1:floor(end/2)}), ...
                    cat(dim, varargin{floor(end/2)+1:end}));
                return
            end
            if isempty(varargin{1})
                c = varargin{2};
                return
            end
            if isempty(varargin{2})
                c = varargin{1};
                return
            end
            a = varargin{1};
            b = varargin{2};
            assert(a.N == b.N);
            c = ZN(a.N, cat(dim, uint8(a), uint8(b)));
        end
        
        function b = conj(a)
            b = ZN(a.N, mod(a.N - a, a.N));
        end
        
        function [charges, vertices] = cumprod(a)
            charges = ZN(a.N, mod(cumsum(a), a.N));
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
        
        function s = GetMD5_helper(data)
            s.N = data.N;
            s.a = uint8(data);
        end
        
        function c = horzcat(varargin)
            c = cat(2, varargin{:});
        end
        
        function c = mtimes(a, b)
            assert(a.N == b.N);
            c = ZN(a.N, mod(a + b, a.N));
        end
        
        function N = Nsymbol(a, b, c)
            assert(a.N == b.N && b.N == c.N)
            N = mod(a + b, a.N) == c;
        end
        
        function e = one(a)
            e = ZN(a.N, uint8(0));
        end
        
        function [d, N] = prod(a, dim)
            if nargin < 2 || isempty(dim)
                dim = find(size(a) ~= 1, 1);
            end
            d = ZN(a.N, mod(sum(a, dim), a.N));
            if nargout > 1
                N = ones(size(d));
            end
        end
        
        function b = repmat(a, varargin)
            b = ZN(a.N, repmat(uint8(a), varargin{:}));
        end
        
        function b = reshape(a, varargin)
            b = ZN(a.N, reshape(uint8(a), varargin{:}));
        end
        
        function R = Rsymbol(a, b, c, ~)
            R = double(Nsymbol(a, b, c));
        end
        
        function varargout = sort(a, varargin)
            [varargout{1:nargout}] = sort(uint8(a), varargin{:});
            varargout{1} = ZN(a.N, varargout{1});
        end
        
        function s = string(a)
            s = compose("%d", uint8(a));
        end
        
        function varargout = subsref(a, s)
            switch s(1).type
                case '()'
                    b = ZN(a.N, subsref(uint8(a), s(1)));
                    if length(s) == 1
                        varargout = {b};
                    else
                        [varargout{1:nargout}] = subsref(b, s(2:end));
                    end
                case '{}'
                    error('Cannot use {} to index.');
                case '.'
                    [varargout{1:nargout}] = builtin('subsref', a, s);
                otherwise
                    error('Undefined behaviour.');
            end
        end
        
        function a = subsasgn(a, s, b)
            switch s(1).type
                case '()'
                    assert(lenght(s) == 1, 'Undefined assignment syntax.');
                    if isempty(a.N)
                        N = b.N
                    else
                        assert(a.N == b.N)
                        N = b.N;
                    end
                    a = ZN(N, subsasgn(uint8(a), s, b));
                case '{}'
                    error('Cannot use {} to assign.');
                case '.'
                    a = builtin('subsasgn', a, s, b);
                otherwise
                    error('Undefined behaviour.');
            end
        end
        
        function b = transpose(a)
            b = ZN(a.N, transpose(uint8(a)));
        end
        
        function c = vertcat(varargin)
            c = cat(1, varargin{:});
        end
        
        function charge = ZN(N, varargin)
            if nargin < 2
                labels = [];
            else
                labels = horzcat(varargin{:});
            end
            charge@uint8(labels);
            if nargin > 0
                charge.N = N;
            end
        end
    end
end

