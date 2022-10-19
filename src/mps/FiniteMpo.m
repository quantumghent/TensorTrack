classdef FiniteMpo
    % Finite Matrix product operators
    
    properties
        L MpsTensor
        O
        R MpsTensor
    end
    
    methods
        function mpo = FiniteMpo(L, O, R)
            if nargin == 0, return; end
            if ~iscell(O), O = {O}; end
            mpo.L = L;
            mpo.O = O;
            mpo.R = R;
        end
        
        function v = apply(mpo, v)
            N = length(mpo);
            for d = 1:depth(mpo)
                if N == 0
                    v = applytransfer(mpo(d).L, mpo(d).R, v);
                elseif N == 1
                    v = applychannel(mpo(d).O{1}, mpo(d).L, mpo(d).R, v);
                else
                    v = applympo(mpo(d).O{:}, mpo(d).L, mpo(d).R, v);
                end
            end
        end
        
        function [V, D, flag] = eigsolve(mpo, v0, howmany, sigma, options)
            arguments
                mpo
                v0 = []
                howmany = 1
                sigma = 'largestabs'
                
                options.Tol = eps(underlyingType(mpo))^(3/4)
                options.Algorithm = 'eigs'
                options.MaxIter = 100
                options.KrylovDim = 20
                options.IsSymmetric logical = false
                options.DeflateDim
                options.ReOrth = 2
                options.NoBuild
                options.Verbosity = 0
            end
            
            if isempty(v0), v0 = initialize_fixedpoint(mpo(1)); end
            kwargs = namedargs2cell(options);
            [V, D, flag] = eigsolve(@(x) mpo.apply(x), v0, howmany, sigma, kwargs{:});
        end
        
        function v = initialize_fixedpoint(mpo)
            v = Tensor.randnc(domain(mpo), []);
        end
        
        function s = domain(mpo)
            s = conj(...
                [rightvspace(mpo(1).L) cellfun(@(x) pspace(x)', mpo(1).O) leftvspace(mpo(1).R)]);
        end
        
        function s = codomain(mpo)
            s = [leftvspace(mpo(end).L) cellfun(@pspace, mpo(end).O) rightvspace(mpo(end).R)];
        end
        
        function d = depth(mpo)
            d = size(mpo, 1);
        end
        
        function l = length(mpo)
            l = length(mpo(1).O);
        end
        
        function type = underlyingType(mpo)
            type = underlyingType(mpo(1).L);
        end
        
        function mpo = ctranspose(mpo)
            if depth(mpo) > 1
                mpo = flip(mpo, 1);
            end
            
            for d = 1:depth(mpo)
                mpo(d).L = mpo(d).L';
                mpo(d).O = cellfun(@ctranspose, mpo(d).O, ...
                    'UniformOutput', false);
                mpo(d).R = mpo(d).R';
            end
        end
        
        function mpo = transpose(mpo)
            if depth(mpo) > 1
                mpo = flip(mpo, 1);
            end
            
            for d = 1:depth(mpo)
                [mpo(d).L, mpo(d).R] = swapvars(mpo(d).L, mpo(d).R);
                mpo(d).O = cellfun(@transpose, ...
                    fliplr(mpo(d).O), 'UniformOutput', false);
            end
        end
        
        function mpo = slice(mpo, i, j)
            if strcmp(i, ':')
                i = 1:size(mpo(end).O, 2);
            end
            if strcmp(j, ':')
                j = 1:size(mpo(1).O{1}, 4);
            end
            assert(all(1 <= i) && all(i <= size(mpo(end).O{1}, 2)));
            assert(all(1 <= j) && all(j <= size(mpo(1).O{1}, 4)));
            
            assert(depth(mpo) == 1, 'TBA');
            mpo.O{1} = mpo.O{1}(:, i, :, j);
        end
        
        function bool = iszero(mpo)
            if isempty(mpo.O)
                bool = false;
                return
            end
            bool = any(cellfun(@nnz, mpo.O) == 0);
        end
        
        function bool = iseye(mpo)
            if isempty(mpo.O)
                bool = true;
                return
            end
            bool = all(cellfun(@iseye, mpo.O));
        end
%         
%         function v = applyleft(mpo, v)
%             
%         end
%         
%         function mpo1 = plus(mpo1, mpo2)
%             
%         end
%         
%         function v = applyright(mpo, v)
%             arguments
%                 mpo
%                 v MpsTensor
%             end
%             
%             assert(depth(mpo) == 1, 'mps:TBA', 'Not implemented yet.');
%             assert(v.plegs == width(mpo), 'mps:ArgError', 'Incompatible sizes.');
%             
%             w = width(mpo);
%             
%             mpopart = cell(2, w);
%             for i = 1:w
%                 mpopart{1, i} = mpo.O{i};
%                 mpopart{2, i} = [2 * i, 2 * (i + 1) + 1, -(1 + i), 2 * i + 1];
%             end
%             
%             v = MpsTensor(contract(...
%                 v, [1, 2:2:(2 * w), 2 * (w + 1), -(1:v.alegs) - (w + 2)], ...
%                 mpo.L, [-1 3 1], ...
%                 mpo.R, [-(w + 2), 2 * (w + 1) + 1, 2 * (w + 1)], ...
%                 mpopart{:}, 'Rank', rank(v)));
%         end
%         
        function t = Tensor(mpo)
            assert(depth(mpo) == 1, 'not implemented for 1 < depth');
            N = length(mpo);
            inds = arrayfun(@(x) [x -(x+1) (x+1) -(N+x+3)], 1:N, ...
                'UniformOutput', false);
            args = [mpo.O; inds];
            t = contract(mpo.L, [-1 1 -(2*N+4)], args{:}, mpo.R, [-(N+3) N+1 -(N+2)], ...
                    'Rank', [N+2 N+2]);
        end
    end
    
    methods (Static)
        function mpo = randnc(pspaces, vspaces)
            assert(length(pspaces) == length(vspaces) + 1);
            
            L = MpsTensor(Tensor.randnc([pspaces(1) vspaces(1)'], pspaces(1)));
            O = cell(1, length(pspaces)-2);
            for i = 1:length(O)
                O{i} = MpoTensor(Tensor.randnc([vspaces(i) pspaces(i+1)], ...
                    [pspaces(i+1) vspaces(i+1)]));
            end
            R = MpsTensor(Tensor.randnc([pspaces(end)' vspaces(end)], pspaces(end)'));
            mpo = FiniteMpo(L, O, R);
        end
    end
end
