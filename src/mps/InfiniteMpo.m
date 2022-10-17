classdef InfiniteMpo
    % Infinite translation invariant matrix product operator.
    
    properties
        O
    end
    
    methods
        function mpo = InfiniteMpo(varargin)
            if nargin == 0, return; end
            
            if nargin == 1
                O = varargin{1};
                
                if isa(O, 'InfiniteMpo')
                    for i = numel(O):1:1
                        mpo(i).O = O(i).O;
                    end
                    mpo = reshape(mpo, size(O));
                    
                elseif isa(O, 'MpoTensor')
                    mpo.O = {O};
                
                elseif isa(O, 'AbstractTensor')
                    for i = height(O):-1:1
                        mpo(i).O = arrayfun(@MpoTensor, O(i, :), 'UniformOutput', false);
                    end
                    
                elseif iscell(O)
                    mpo.O = O;
                end
                return
            end
        end
    end
    
    methods
        function [GL, lambda] = leftenvironment(mpo, mps1, mps2, GL, eigopts)
            arguments
                mpo
                mps1
                mps2
                GL = []
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            if isempty(GL)
                GL = cell(1, period(mps1));
                for i = size(mpo.O{1}, 1):-1:1
                    GL{1}(1, i, 1) = Tensor.randnc(...
                        [leftvspace(mps2, 1) leftvspace(mpo.O{1}, i)'], leftvspace(mps1, 1));
                end
            end
            
            [GL{1}, lambda] = eigsolve(@(x) transferleft(mpo, mps1, mps2, x), ...
                GL{1}, 1, 'largestabs', 'KrylovDim', eigopts.KrylovDim, 'Tol', eigopts.Tol, ...
                'ReOrth', eigopts.ReOrth, 'MaxIter', eigopts.MaxIter);
            
            n = lambda^(1/period(mps1));
            for w = 2:period(mps1)
                GL{w} = applyleft(mpo.O{w}, mps1.AL(w), conj(mps2.AL(w)), GL{w-1}) / n;
            end
        end
        
        function [GR, lambda] = rightenvironment(mpo, mps1, mps2, GR, eigopts)
            arguments
                mpo
                mps1
                mps2
                GR = []
                eigopts.KrylovDim = 30
                eigopts.MaxIter = 1000
                eigopts.ReOrth = 2
                eigopts.Tol = eps(underlyingType(mps1))^(3/4)
            end
            
            N = period(mps1);
            if isempty(GR)
                GR = cell(1, N);
                for i = size(mpo.O{end}, 3):-1:1
                    GR{1}(1, i, 1) = Tensor.randnc(...
                        [rightvspace(mps1, N)' rightvspace(mpo.O{end}, i)'], ...
                        rightvspace(mps2, N)');
                end
            end
            
            [GR{1}, lambda] = eigsolve(@(x) transferright(mpo, mps1, mps2, x), ...
                GR{1}, 1, 'largestabs', 'KrylovDim', eigopts.KrylovDim, 'Tol', eigopts.Tol, ...
                'ReOrth', eigopts.ReOrth, 'MaxIter', eigopts.MaxIter);
            
            n = lambda^(1/N);
            for w = N:-1:1
                GR{w} = applyright(mpo.O{w}, mps1.AR(w), conj(mps2.AR(w)), ...
                    GR{next(w, N)}) / n;
            end
        end
        
        function x = transferleft(mpo, mps1, mps2, x)
            for w = 1:period(mps1)
                x = applyleft(mpo.O{w}, mps1.AL(w), conj(mps2.AL(w)), x);
            end
        end
        
        function x = transferright(mpo, mps1, mps2, x)
            for w = period(mps1):-1:1
                x = applyright(mpo.O{w}, mps1.AR(w), conj(mps2.AR(w)), x);
            end
        end
    end
    
    methods (Static)
        function mpo = Ising(beta, kwargs)
            arguments
                beta = log(1 + sqrt(2)) / 2;
                kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'Z2'})} = 'Z1'
            end
            
            if strcmp(kwargs.Symmetry, 'Z1')
                t = [exp(beta) exp(-beta); exp(-beta) exp(beta)];
                [v, d] = eig(t);
                t = v * sqrt(d) * v;
                
                o = zeros(2, 2, 2, 2);
                o(1, 1, 1, 1) = 1;
                o(2, 2, 2, 2) = 1;
                
                o = contract(o, 1:4, t, [-1 1], t, [-2 2], t, [-3 3], t, [-4 4]);
                
                O = fill_tensor(Tensor.zeros([2 2 2 2]), o);
                
            else
                s = GradedSpace.new(Z2(0, 1), [1 1], false);
                O = fill_tensor(Tensor([s s], [s s]), ...
                    @(~,f) sqrt(prod(logical(f.uncoupled) .* sinh(beta) + ...
                    ~logical(f.uncoupled) .* cosh(beta))));
            end
            
            mpo = InfiniteMpo(O);
        end
    end
end
