classdef GtPattern
    %GTPATTERN Object that represents a pattern devised by Gelfand and Tsetlin.
    %
    %        /m_{1,N} m_{2,N}  ...  m_{N,N}\
    %       |  m_{1,N-1}  ...  m_{N-1,N-1}  |
    %   M = |             ...               |
    %       |       m_{1,2} m_{2,2}         |
    %        \          m_{1,1}            /
    %
    %   These consist of triangular arrangements of integers M_ij with the
    %   following property:
    %       - for 1 <= k < l <= N:
    %           m_{k,l} >= m_{k,l-1} >= m_{k+1,l}
    %
    %   They will be represented using arrays of size N x N, where the elements
    %   outside of the triangular region are assumed to be zero.
    
    
    properties % (Access = private)
        M double
        N
    end
    
    methods
        function p = GtPattern(N, M)
            if nargin == 0
                return
            end
            if nargin == 1
                M = [];
            else
                assert(2 * size(M, 1) == N * (N + 1), ...
                    'Wrong size for M: %d %d %d', N, size(M, 1), size(M, 2));
            end
            
            p.N = N;
            p.M = M;
        end
        
        function [rows, cols] = size(p)
            sz = [1 size(p.M, 2)];
            if nargout <= 1
                rows = sz;
                return
            end
            rows = sz(1);
            cols = sz(2);
        end
        
        function l = length(p)
            l = size(p.M, 2);
        end
        
        function a = cat(dim, varargin)
            if nargin == 2
                a = varargin{1};
                return
            end
            
            if nargin > 3
                a = cat(dim, cat(dim, varargin{1:floor(end/2)}), ...
                    cat(dim, varargin{floor(end/2)+1:end}));
                return
            end
            
            a = varargin{1};
            b = varargin{2};
            a.M = horzcat(a.M, b.M);
        end
        
        function a = horzcat(varargin)
            a = cat(1, varargin{:});
        end
        
        function a = subsref(a, s)
            assert(length(s) == 1)
            switch s.type
                case '()'
                    
                    if length(s.subs) == 1
                        a.M = a.M(:, s.subs{1});
                    elseif length(s.subs) == 2
                        assert(s.subs{1} == 1 || s.subs{1} == ':');
                        a.M = a.M(:, s.subs{2});
                    else
                        error('Undefined.');
                    end
                    
                case '.'
                    assert(length(s) == 1)
                    switch s.subs
                        case 'M'
                            a = a.M;
                        case 'N'
                            a = a.N;
                        otherwise
                            error('Undefined behaviour.');
                    end
                otherwise
                    error('Undefined behaviour.');
            end
        end
        
        function bool = checkbounds(pat, k, l)
            bool = all(0 < k & k <= l & l <= pat.N);
        end
        
        function m = get(p, k, l)
            %GET Access an element in the pattern.
            %   m = get(pat, k, l)
            %       gets the element specified by m_kl.
            %
            %   m = get(pat, ks, l)
            %       gets a row vector of elements in the l'th row.
            %
            %   m = get(pat, k, ls(:)
            %       gets a col vector of elements in the k'th column.
%             assert(isscalar(p));
            m = p.M(k + bitshift(((l + 1 + p.N) .* (p.N - l)), -1));
        end
        
        function p = set(p, k, l, val)
%             assert(isscalar(p));
            assert(checkbounds(p, k, l));
            p.M(k + bitshift(((l + 1 + p.N) * (p.N - l)), -1)) = val;
        end
        
        function bool = isscalar(a)
            bool = size(a.M, 2) == 1;
        end
        
        function disp(p)
            N = p.N;
            M = p.M;
            for i = 1:size(M, 2)
                %% Special case N == 1
                if N == 1
                    fprintf('(%d)', M(:,i));
                    return
                end
                
                
                %% Other cases
                dig = maxdigits(M(:,i));
                space = repmat(' ', 1, dig);
                out = strings(N, 2*N+2);
                
                % matrix delimiters
                % left parenthesis: char(40)
                % left upper hook: char(9115)
                % left extension:  char(9116)
                % left lower hook: char(9117)
                %
                % right parenthesis: char(41)
                % right upper hook: char(9118)
                % right extension:  char(9119)
                % right lower hook: char(9120)
                
                out(:,1) = vertcat("\t/ ", repmat("\t| ",N-2,1), "\t\\ ");
                out(:,end-1) = vertcat(" \\", repmat(" |", N-2, 1), " /");
                out(:, end) = '\n';
                out(:,2:end-2) = space;
                
                % fill in data
                ctr = 0;
                for ii = 1:N
                    l = N - ii + 1;
                    jj = 1;
                    for k = 1:l
                        out(ii, ii+jj) = sprintf('%*d', dig, M(end-ctr, i));
                        ctr = ctr + 1;
                        jj = jj + 2;
                    end
                end
                assert(ctr == size(M, 1));
                out = join(out, '');
                out = join(out, '');
                fprintf(out);
            end
        end
        
        function bool = lt(a, b)
            assert(isscalar(a) && isscalar(b));
            for l = a.N:-1:1
                for k = 1:l
                    if get(b, k, l) > get(a, k, l)
                        bool = true;
                        return
                    elseif get(b, k, l) < get(a, k, l)
                        bool = false;
                        return
                    end
                end
            end
            bool = false;
        end
        
        function bool = eq(a, b)
            bool = all(a.M == b.M, 1);
            bool = reshape(bool, 1, []);
        end
        
        function sigma = rowsum(p, l)
            %ROWSUM Compute the sum of row `l`.
            %   sigma = rowsum(pat, l)
            %       computes sigma_l = \sum_{k=1:l} m_{k, l}
            
            sigma = sum(get(p, 1:l, l));
        end
        
        function w = pWeight(p)
            %PWEIGHT Compute the p-weight of a pattern.
            %   w = pWeight(pat)
            %       computes the pattern weight W(pat) = (w_1 w_2 ... w_N) where
            %       w_i = sigma_i - sigma_{i-1}.
            %
            %   This is an alternative weight definition to the z-weight.
            %   See also GtPattern.zWeight
            
            M = p.M;
            ctr = 0;
            sigmas = zeros(p.N + 1, length(p));
            for i = 1:p.N
                sigmas(i + 1, :) = sum(M(end - (1:i) - ctr + 1, :), 1);
                ctr = ctr + i;
            end
            w = diff(sigmas).';
%             sigma = [0 arrayfun(@(l) rowsum(p, l), 1:p.N)]
%             w = sigma(2:end) - sigma(1:end-1);
        end
        
        function w = zWeight(p)
            %ZWEIGHT Compute the z-weight of a pattern
            %   w = zWeight(pat)
            %       computes the z-weight L(pat) = (l_1 l_2 ... l_N-1) where
            %       l_i = sigma_l - 1/2(sigma_{l+1} + sigma_{l-1})
            %
            %   This is a generalization of the m-quantum number for angular
            %   momentum.

            sigma = [0 arrayfun(@(l) rowsum(p, l), 1:p.N)];
            w = sigma(2:end-1) - (sigma(1:end-2) + sigma(3:end))/2;
        end
        
        
    end
end


function d = maxdigits(var)
%MAXDIGITS Compute max amount of digits used by positive integers in an array.
%   d = maxdigits(var)

d = 1;
var = var(var ~= 0);
d = max(max(1, 1+floor(log10(double(var)))));

end

function i = lindex(N, k, l)
%LINDEX Give the linear index of element (k,l).
%   i = lindex(N, k, l)
%       computes the linear index of the k,l'th element in a pattern.

i = k + idivide((l + 1 + N) * (N-l), int16(2));

end
