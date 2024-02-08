classdef O2 < AbstractCharge
    % Irreducible representations of :math:`\mathrm{O}(2)`.
    %
    %   This class represents the representations of :math:`\mathrm{O}(2)`, also known as the
    %   semi-direct product of :math:`\mathrm{U}(1)` and charge conjugation.
    %
    %   The representations are labeled using :math:`(j, s)`, indicating the behaviour
    %   of :math:`\mathrm{U}(1)` and charge conjugation respectively. This leads to two
    %   1-dimensional representations :math:`(0, 0)` and :math:`(0, 1)`, and for any
    %   non-negative :math:`j` a two-dimensional representation :math:`(j, 2)`.
    %
    % Properties
    % ----------
    % j : :class:`uint8`
    %   :math:`\mathrm{U}(1)` label
    %
    % s : :class:`uint8`
    %   indicator for type of representation
    
    properties
        j uint8
        s uint8
    end
    
    methods
        function style = braidingstyle(~)
            style = BraidingStyle.Bosonic;
        end
        
        function a = conj(a), end
        
        function bool = eq(a, b)
            bool = reshape([a.j], size(a)) == reshape([b.j], size(b)) & ...
                reshape([a.s], size(a)) == reshape([b.s], size(b));
        end
        
        function nu = frobeniusschur(a)
            nu = ones(size(a));
        end
        
        function style = fusionstyle(~)
            style = FusionStyle.Simple;
        end
        
        function F = Fsymbol(a, b, c, d, e, f)
            if numel(a) > 1
                F = reshape(arrayfun(@Fsymbol, a, b, c, d, e, f), size(a));
                return
            end

            if (~Nsymbol(a, b, e) || ~Nsymbol(e, c, d) || ...
                    ~Nsymbol(b, c, f) || ~Nsymbol(a, f, d))
                F = 0;
                return
            end
            
            % Cases where a, b, c = 0+ or 0-
            if (a.j == 0 && a.s == 0) || (b.j == 0 && b.s == 0) || (c.j == 0 && c.s == 0)
                F = 1;
                return
            end
            
            if ((a.j == 0 && a.s == 1) && (b.j == 0 && b.s == 1)) || ...
                    ((a.j == 0 && a.s == 1) && (c.j == 0 && c.s == 1)) || ...
                    ((b.j == 0 && b.s == 1) && (c.j == 0 && c.s == 1))
                F = 1;
                return
            end
            
            if a.j == 0 && a.s == 1
                if d.j == 0
                    F = 1;
                elseif d.j == c.j - b.j
                    F = -1;
                else
                    F = 1;
                end
                return
            end
            
            if b.j == 0 && b.s == 1
                if d.j == max(a.j - c.j, c.j - a.j)
                    F = -1;
                else
                    F = 1;
                end
                return
            end
            
            if c.j == 0 && c.s == 1
                if d.j == a.j - b.j
                    F = -1;
                else
                    F = 1;
                end
                return
            end
            
            % Other cases
            isqrt2 = sqrt(2) / 2;
            if a == b && b == c
                if d == a
                    if e.j == 0
                        if f.j == 0
                            if f.s == 1
                                F = -0.5;
                            else
                                F = 0.5;
                            end
                        else
                            if e.s == 1
                                F = -isqrt2;
                            else
                                F = isqrt2;
                            end
                        end
                    else
                        if f.j == 0
                            F = isqrt2;
                        else
                            F = 0;
                        end
                    end
                else
                    F = 1;
                end
                return
            end
            
            if a == b
                if d == c
                    if f.j == b.j + c.j
                        if e.s == 1
                            F = -isqrt2;
                        else
                            F = isqrt2;
                        end
                    else
                        F = isqrt2;
                    end
                else
                    F = 1;
                end
                return
            end
            
            if b == c
                if d == a
                    if e.j == a.j + b.j
                        F = isqrt2;
                    else
                        if f.s == 1
                            F = -isqrt2;
                        else
                            F = isqrt2;
                        end
                    end
                else
                    F = 1;
                end
                return
            end
            
            if a == c
                if d == b
                    if e.j == f.j
                        F = 0;
                    else
                        F = 1;
                    end
                else
                    if d.s == 1
                        F = -1;
                    else
                        F = 1;
                    end
                end
                return
            end
            
            if d.j == 0 && d.s == 1
                if b.j == a.j + c.j
                    F = -1;
                else
                    F = 1;
                end
                return
            end
            F = 1;
        end
        
        function C = fusiontensor(a, b, c)
            C = zeros(qdim(a), qdim(b), qdim(c), 1);
            if ~Nsymbol(a, b, c), return; end
            if c.j == 0
                if a.j == 0 && b.j == 0
                    C(1, 1, 1, 1) = 1;
                else
                    if c.s == 0
                        C(1, 2, 1, 1) = 1 / sqrt(2);
                        C(2, 1, 1, 1) = 1 / sqrt(2);
                    else
                        C(1, 2, 1, 1) = 1 / sqrt(2);
                        C(2, 1, 1, 1) = -1 / sqrt(2);
                    end
                end
            elseif a.j == 0
                C(1, 1, 1, 1) = 1;
                if a.s == 1
                    C(1, 2, 2, 1) = -1;
                else
                    C(1, 2, 2, 1) = 1;
                end
            elseif b.j == 0
                C(1, 1, 1, 1) = 1;
                if b.s == 1
                    C(2, 1, 2, 1) = -1;
                else
                    C(2, 1, 2, 1) = 1;
                end
            elseif c.j == a.j + b.j
                C(1, 1, 1, 1) = 1;
                C(2, 2, 2, 1) = 1;
            elseif c.j == a.j - b.j
                C(1, 2, 1, 1) = 1;
                C(2, 1, 2, 1) = 1;
            elseif c.j == b.j - a.j
                C(2, 1, 1, 1) = 1;
                C(1, 2, 2, 1) = 1;
            end
        end
        
        function s = GetMD5_helper(data)
            s = [data.j] + [data.s];
        end
        
        function bool = issortedrows(A)
            bool = issortedrows(reshape([A.j] + [A.s], size(A)));
        end
        
        function bool = issorted(A)
            bool = issorted(reshape([A.j] + [A.s], size(A)));
        end
        
        function c = mtimes(a, b)
            assert(isscalar(a) && isscalar(b), ...
                'Simple fusion cannot be vectorised.');
            if a.j == 0 && b.j == 0
                c = O2(0, xor(a.s, b.s));
            elseif a.j == 0
                c = b;
            elseif b.j == 0
                c = a;
            elseif a.j == b.j
                c = O2([0 0 a.j + b.j], [0 1 2]);
            else
                c = O2([a.j + b.j, max(a.j - b.j, b.j - a.j)], [2 2]);
            end
        end
        
        function bool = ne(a, b)
            bool = reshape([a.j], size(a)) ~= reshape([b.j], size(b)) | ...
                reshape([a.s], size(a)) ~= reshape([b.s], size(b));
        end
        
        function N = Nsymbol(a, b, c)
            N = double(reshape(...
                ([c.s] == 0 & [a.j] == [b.j] & [a.s] == [b.s]) | ...
                ([c.s] == 1 & [a.j] == [b.j] & ([a.s] ~= [b.s] | [a.s] == 2)) | ...
                ([c.s] == 2 & ([c.j] == ([a.j] + [b.j]) |  ...
                ([c.j] == max([a.j] - [b.j], [b.j] - [a.j])))), ...
                size(a)));
        end
        
        function charges = O2(j, s)
            if nargin == 0, return; end
            assert(all(size(j) == size(s)));
            if ~true
                for i = numel(j):-1:1
                    charges(i).j = j(i);
                    charges(i).s = s(i);
                end
            else
                j_cell = num2cell(j);
                s_cell = num2cell(s);
                charges = repmat(O2, size(j));
                [charges.j] = deal(j_cell{:});
                [charges.s] = deal(s_cell{:});
            end
            charges = reshape(charges, size(j));
        end
        
        function d = qdim(a)
            d = ones(size(a));
            d([a.j] ~= 0) = 2;
        end
        
        function R = Rsymbol(a, b, c, inv)
            if nargin > 3 && inv
                R = Rsymbol(b, a, c);
                return
            end
            R = Nsymbol(a, b, c);
            R([c.s] == 1 & [a.j] > 0) = -R([c.s] == 1 & [a.j] > 0);
        end
        
        function [B, I] = sort(A, varargin)
            [~, I] = sort(reshape([A.j] + [A.s], size(A)), varargin{:});
            B = A(I);
        end
        
        function str = string(a)
            str = reshape(string([a.j]), size(a));
            str([a.s] == 0) = str([a.s] == 0) + "+";
            str([a.s] == 1) = str([a.s] == 1) + "-";
        end
        
        function e = one(~)
            e = O2(0, 0);
        end
    end
end
