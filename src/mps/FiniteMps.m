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
        
        function mps = movegaugecenter(mps, newcenter, alg)
            arguments
                mps
                newcenter {mustBeInteger} = floor(length(mps) / 2)
                alg = 'polar'
            end
            
            assert(newcenter >= 1 && newcenter <= length(mps), ...
                'mps:gauge', ...
                sprintf('new center (%d) must be in 1:%d', newcenter, length(mps)));
            
            if isempty(mps.center)
                low = 1;
                high = length(mps);
            else
                low = mps.center;
                high = mps.center;
            end
            
            for i = low:(newcenter - 1)
                [mps.A(i), L] = leftorth(mps.A(i), alg);
                mps.A(i + 1) = multiplyleft(mps.A(i + 1), L);
            end
            for i = high:-1:(newcenter + 1)
                [R, mps.A(i)] = rightorth(mps.A(i), alg);
                [mps.A(i - 1)] = multiplyright(mps.A(i - 1), R);
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
        
        function rho = fixedpoint(mps, type, w)
            arguments
                mps
                type {mustBeMember(type, {'l', 'r'})}
                w = floor(length(mps) / 2);
            end
            
            if strcmp(type, 'l')
                if mps.center > w
                    rho = mps.A.eye(leftvspace(mps, w), leftvspace(mps, w));
                else
                    T = transfermatrix(mps, mps, mps.center:w);
                    rho = apply(T, []);
                end
            else
                if mps.center < w
                    rho = mps.A.eye(rightvspace(mps, w)', rightvspace(mps, w)');
                else
                    T = transfermatrix(mps, mps, w:mps.center).';
                    rho = apply(T, []);
                end
            end
                
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
            if isempty(mps.center), mps = movegaugecenter(mps); end
            n = norm(mps);
            mps.A(mps.center) = mps.A(mps.center) ./ n;
        end
        
        function [svals, charges] = schmidt_values(mps, w)
            arguments
                mps
                w {mustBeInteger} = floor(length(mps) / 2)
            end
            
            assert(0 <= w && w <= length(mps));
            if isempty(mps.center), mps = movegaugecenter(mps, w); end
            
            if w < mps.center
                mps = movegaugecenter(mps, w + 1);
                S = tsvd(mps.A(w + 1), 1, 2:nspaces(mps.A(w + 1)));
            else
                mps = movegaugecenter(mps, w);
                S = tsvd(mps.A(w), 1:nspaces(mps.A(w))-1, nspaces(mps.A(w)));
            end
            [svals, charges] = matrixblocks(S);
            svals = cellfun(@diag, svals, 'UniformOutput', false);
        end
        
        function psi = Tensor(mps)
            tensors = num2cell(mps.A);
            indices = cell(size(tensors));
            
            nout = nspaces(mps.A(1)) - 1;
            indices{1} = [-(1:nout), 1];
            
            for i = 2:length(indices)-1
                plegs = mps.A(i).plegs;
                indices{i} = [i-1 -(1:plegs)-nout i];
                nout = nout + plegs;
            end
            
            plegs = mps.A(end).plegs;
            indices{end} = [length(indices)-1 -(1:plegs+1)-nout];
            
            args = [tensors; indices];
            psi = contract(args{:});
        end
    end
    
    methods (Static)
        function mps = new(fun, pspaces, kwargs)
            arguments
                fun
            end
            arguments (Repeating)
                pspaces
            end
            arguments
                kwargs.L
                kwargs.MaxVspace
                kwargs.LeftVspace = one(pspaces{1})
                kwargs.RightVspace = one(pspaces{1})
            end
            
            if isempty(fun), fun = @randnc; end
            
            if isfield(kwargs, 'L')
                pspaces = repmat(pspaces, 1, kwargs.L);
            end
            if isfield(kwargs, 'MaxVspace') && ~iscell(kwargs.MaxVspace)
                kwargs.MaxVspace = {kwargs.MaxVspace};
            end
            
            L = length(pspaces);
            
            vspacesL = cell(1, L + 1);
            vspacesL{1} = kwargs.LeftVspace;
            for i = 1:L
                vspacesL{i+1} = vspacesL{i} * pspaces{i};
                if isfield(kwargs, 'MaxVspace')
                    vspacesL{i + 1} = infimum(vspacesL{i + 1}, ...
                        kwargs.MaxVspace{mod1(i + 1, length(kwargs.MaxVspace))});
                end
            end
            
            vspacesR = cell(1, L + 1);
            vspacesR{end} = kwargs.RightVspace;
            for i = L:-1:1
                vspacesR{i} = pspaces{i}' * vspacesR{i+1};
                if isfield(kwargs, 'MaxVspace')
                    vspacesR{i} = infimum(vspacesR{i}, ...
                        kwargs.MaxVspace{mod1(i, length(kwargs.MaxVspace))});
                end
            end
            
            vspaces = cellfun(@infimum, vspacesL, vspacesR, 'UniformOutput', false);
            
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
        
        
        function mps = new_from_spaces(fun, spaces, kwargs)
            arguments
                fun
            end
            arguments (Repeating)
                spaces
            end
            arguments
                kwargs.MaxVspace
            end
            
            if isfield(kwargs, 'MaxVspace') && ~iscell(kwargs.MaxVspace)
                kwargs.MaxVspace = {kwargs.MaxVspace};
            end
            
            assert(mod(length(spaces) - 1, 2));
            if isempty(fun), fun = @randnc; end
            
            L = (length(spaces) - 1) / 2;
            
            for w = L:-1:1
                Vleft = spaces{2 * w - 1};
                P = spaces{2 * w};
                Vright = spaces{2 * w + 1};
                rankdeficient = Vleft * P < Vright || Vleft > P' * Vright;
                if rankdeficient
                    error('mps:rank', ...
                        'Cannot create a full rank mps with given spaces.');
                end
                A{w} = Tensor.new(fun, [Vleft P], Vright);
            end
            
            mps = FiniteMps(A);
        end
        
        function mps = randnc(varargin)
            mps = FiniteMps.new(@randnc, varargin{:});
        end
    end
end

