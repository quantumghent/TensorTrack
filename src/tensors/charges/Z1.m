classdef Z1 < AbstractCharge
    % Trivial charges.
    
    methods
        function A = Asymbol(~, ~, ~)
            A = 1;
        end
        
        function B = Bsymbol(~, ~, ~)
            B = 1;
        end
        
        function c = char(a)
            c = '0';
        end
        
        function a = conj(a)
        end
        
        function bool = eq(~, ~)
            bool = true;
        end
        
        function nu = frobeniusschur(~)
            nu = 1;
        end
        
        function F = Fsymbol(~, ~, ~, ~, ~, ~)
            F = 1;
        end
        
        function C = fusiontensor(~, ~, ~)
            C = 1;
        end
        
        function a = intersect(a, ~)
        end
        
        function bool = issortedrows(~)
            bool = true;
        end
        
        function a = mtimes(a, ~), end
        
        function bools = ne(A, B)
            if isscalar(A)
                bools = false(size(B));
            elseif isscalar(B)
                bools = false(size(A));
            else
                assert(all(size(A) == size(B)), 'tensors:charges:SizeError');
                bools = false(size(A));
            end
        end
        
        function N = Nsymbol(~, ~, ~)
            N = 1;
        end
        
        function R = Rsymbol(~, ~, ~, ~)
            R = 1;
        end
        
        function varargout = sort(a, varargin)
            if nargout == 1
                varargout = {a};
            elseif nargout == 2
                varargout = {a, 1:numel(a)};
            end
        end
        
        function s = string(a)
            s = repmat("0", size(a));
        end
        
        function style = braidingstyle(~)
            style = BraidingStyle.Bosonic;
        end
        
        function style = fusionstyle(~)
            style = FusionStyle.Unique;
        end
        
        function e = one(~)
            e = Z1;
        end
    end
end
