classdef FiniteMps
    % Finite Matrix product states
    
    properties
        A MpsTensor
        center
    end
    
    methods
        function mps = FiniteMps(varargin)
            
            if nargin == 0, return; end
            if nargin == 1
                mps.A = varargin{1};
            elseif nargin == 2
                mps.A = varargin{1};
                mps.center = varargin{2};
            end
        end
        
        function p = length(mps)
            p = length(mps.A);
        end
        
        function mps = movecenter(mps, newcenter, alg)
            arguments
                mps
                newcenter {mustBeInteger} = floor(length(mps) / 2)
                alg = 'polar'
            end
            
            assert(newcenter >= 1 && newcenter <= length(mps), ...
                'mps:gauge', ...
                sprintf('new center (%d) must be in 1:%d', newcenter, length(mps)));
            
            for i = max(mps.center, 1):(newcenter - 1)
                [mps.A{i}, L] = leftorth(mps.A{i}, alg);
                mps.A{i + 1} = multiplyleft(mps.A{i + 1}, L);
            end
            for i = min(mps.center, length(mps)):-1:(newcenter + 1)
                [R, mps.A{i}] = rightorth(mps.A{i}, alg);
                mps.A{i - 1} = multiplyright(mps.A{i - 1}, R);
            end
            mps.center = newcenter;
        end
        
        function T = transfermatrix(mps1, mps2, sites)
            arguments
                mps1
                mps2 = mps1
                sites = 1:length(mps1)
            end
            
            assert(all(diff(sites) == 1), 'sites must be neighbouring and increasing.');
            T = transfermatrix(mps1.A(sites), mps2.A(sites));
        end
        
        function o = overlap(mps1, mps2, rholeft, rhoright)
            arguments
                mps1
                mps2
                rholeft = []
                rhoright = []
            end
            
            assert(length(mps1) == length(mps2), 'mps should have equal length');
            assert(length(mps1) >= 2, 'edge case to be added');
            
            n = floor(length(mps1) / 2);
            Tleft = transfermatrix(mps1, mps2, 1:n);
            rholeft = apply(Tleft, rholeft);
            Tright = transfermatrix(mps1, mps2, n+1:length(mps1)).';
            rhoright = apply(Tright, rhoright);
            
            o = overlap(rholeft, rhoright);
        end
        
        function n = norm(mps)
            if isempty(mps.center)
                n = sqrt(abs(overlap(mps, mps)));
                return
            end
            
            n = norm(mps.A(mps.center));
        end
        
        function [mps, n] = normalize(mps)
            if isempty(mps.center), mps = movecenter(mps); end
            n = norm(mps);
            mps.A(mps.center) = mps.A(mps.center) ./ n;
        end
        
        function [svals, charges] = schmidt_values(mps, w)
            arguments
                mps
                w {mustBeInteger} = floor(length(mps) / 2)
            end
            
            assert(0 <= w && w <= length(mps));
            if isempty(mps.center), mps = movecenter(mps, w); end
            
            if w < mps.center
                mps = movecenter(mps, w + 1);
                S = tsvd(mps.A(w + 1), 1, 2:nspaces(mps.A(w + 1)));
            else
                mps = movecenter(mps, w);
                S = tsvd(mps.A(w), 1:nspaces(mps.A(w))-1, nspaces(mps.A(w)));
            end
            [svals, charges] = matrixblocks(S);
            svals = cellfun(@diag, svals, 'UniformOutput', false);
        end
    end
    
    methods (Static)
        function mps = new(fun, vspace1, pspaces, vspaces)
            arguments
                fun
                vspace1
            end
            arguments (Repeating)
                pspaces
                vspaces
            end
            
            if isempty(fun), fun = @randnc; end
            
            L = length(pspaces);
            vspaces = [{vspace1} vspaces];
            
            for w = length(pspaces):-1:1
                rankdeficient = vspaces{w} * pspaces{w} < vspaces{w + 1} || ...
                    vspaces{w} > pspaces{w}' * vspaces{w + 1};
                if rankdeficient
                    error('mps:rank', ...
                        'Cannot create a full rank mps with given spaces.');
                end
                
                A{w} = Tensor.new(fun, [vspaces{w} pspaces{w}], vspaces{w + 1});
            end
            mps = FiniteMps(A);
        end
        
        function mps = randnc(varargin)
            mps = FiniteMps.new(@randnc, varargin{:});
        end
    end
end

