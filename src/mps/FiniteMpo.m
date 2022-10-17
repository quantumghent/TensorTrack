classdef FiniteMpo
    % Finite Matrix product operators
    
    properties
        L
        O
        R
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
            assert(nspaces(v) - 2 == N, 'incompatible vector.');

            inds = arrayfun(@(x) [2*x -(x+1) 2*(x+1) 2*x+1], 1:N, ...
                'UniformOutput', false);
            for d = depth(mpo):-1:1
                args = [mpo(d).O; inds];
                v = contract(v, 1:2:2*N+3, ...
                    mpo(d).L, [-1 2 1], args{:}, mpo(d).R, [2*N+3 2*N+2 -(N+2)], ...
                    'Rank', rank(v));
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
            
            kwargs = [fieldnames(options).'; struct2cell(options)];
            [V, D, flag] = eigsolve(@(x) mpo.apply(x), v0, howmany, sigma, kwargs{:});
        end
        
        function v = initialize_fixedpoint(mpo)
            v = Tensor.randnc(domain(mpo), []);
        end
        
        function s = domain(mpo)
            s = conj(...
                [space(mpo(1).L, 3), cellfun(@(x) space(x, 4), mpo.O), space(mpo.R, 1)]);
        end
        
        function s = codomain(mpo)
            s = [space(mpo(1).L, 1), cellfun(@(x) space(x, 2), mpo.O), space(mpo.R, 3)];
        end
        
        function d = depth(mpo)
            d = size(mpo, 1);
        end
        
        function l = length(mpo)
            l = length(mpo(1).O);
        end
        
        function w = width(mpo)
            w = size(mpo.O, 2);
        end
        
        function v = applyleft(mpo, v)
            
        end
        
        function mpo_d = ctranspose(mpo)
            mpo_d = mpo;
            mpo_d.L = tpermute(conj(mpo.L), [3 2 1], rank(mpo.L));
            mpo_d.R = tpermute(conj(mpo.R), [3 2 1], rank(mpo.R));
            mpo_d.O = cellfun(@(x) tpermute(conj(x), [3 2 1 4], rank(x)), mpo.O, ...
                'UniformOutput', false);
%             mpo_d = FiniteMpo(conj(mpo.L), conj(mpo.O), conj(mpo.R));
        end
        
        function mpo1 = plus(mpo1, mpo2)
            
        end
        
        function v = applyright(mpo, v)
            arguments
                mpo
                v MpsTensor
            end
            
            assert(depth(mpo) == 1, 'mps:TBA', 'Not implemented yet.');
            assert(v.plegs == width(mpo), 'mps:ArgError', 'Incompatible sizes.');
            
            w = width(mpo);
            
            mpopart = cell(2, w);
            for i = 1:w
                mpopart{1, i} = mpo.O{i};
                mpopart{2, i} = [2 * i, 2 * (i + 1) + 1, -(1 + i), 2 * i + 1];
            end
            
            v = MpsTensor(contract(...
                v, [1, 2:2:(2 * w), 2 * (w + 1), -(1:v.alegs) - (w + 2)], ...
                mpo.L, [-1 3 1], ...
                mpo.R, [-(w + 2), 2 * (w + 1) + 1, 2 * (w + 1)], ...
                mpopart{:}, 'Rank', rank(v)));
        end
        
        function t = Tensor(mpo)
            assert(depth(mpo) == 1, 'not implemented for 1 < depth');
            N = length(mpo);
            inds = arrayfun(@(x) [x -(x+1) (x+1) -(N+x+3)], 1:N, ...
                'UniformOutput', false);
            args = [mpo.O; inds];
            t = contract(mpo.L, [-1 1 -(N+3)], args{:}, mpo.R, [-(2*N+4) N+1 -(N+2)], ...
                    'Rank', [N+2 N+2]);
        end
    end
    
    methods (Static)
        
    end
end

