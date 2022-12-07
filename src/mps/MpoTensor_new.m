classdef (InferiorClasses = {?Tensor, ?MpsTensor, ?SparseTensor}) MpoTensor < AbstractTensor
    % Matrix product operator building block
    %   This object represents the MPO tensor at a single site as the sum of rank (2,2)
    %   (sparse) tensors and some scalars, which will be implicitly used as unit tensors.
    %
    %                               4
    %                               ^
    %                               |
    %                        1 ->-- O -->- 3
    %                               |
    %                               ^
    %                               2

    properties
        tensors = []
        scalars = []
    end
    
    methods
        function n = nspaces(~)
            n = 4;
        end
    end
    
    methods
        function t = MpoTensor(varargin)
            if nargin == 0, return; end
            
            if nargin == 1
                if isnumeric(varargin{1})
                    t.scalars = varargin{1};
                    t.tensors = SparseTensor([], [], size(t.scalars));
                elseif isa(varargin{1}, 'MpoTensor')
                    t.scalars = varargin{1}.scalars;
                    t.tensors = varargin{1}.tensors;
                elseif isa(varargin{1}, 'Tensor') || isa(varargin{1}, 'SparseTensor')
                    t.tensors = varargin{1};
                    t.scalars = zeros(size(t.tensors));
                end
                return
            end
            
            if nargin == 2
                t.tensors = varargin{1};
                t.scalars = varargin{2};
                N = max(ndims(t.tensors), ndims(t.scalars));
                assert(isequal(size(t.tensors, 1:N), size(t.scalars, 1:N)), 'mpotensors:dimerror', ...
                    'scalars and tensors should have the same size.');
                if any(ismember(find(t.tensors), find(t.scalars)))
                    warning('mpotensor contains both scalar and tensor at the same site.');
                end
            end
        end
        
        function t = minus(t, t2)
            t.tensors = t.tensors - t2.tensors;
            t.scalars = t.scalars - t2.scalars;
        end
        
        function t = plus(t, t2)
            t.tensors = t.tensors + t2.tensors;
            t.scalars = t.scalars + t2.scalars;
        end
        
        function t = times(t, a)
            if isnumeric(t)
                t = a .* t;
                return
            end
            
            assert(isnumeric(a), 'mpotensor:argerror', 'invalid input types.');
            
            t.tensors = t.tensors .* a;
            t.scalars = t.scalars .* a;
        end
        
        function t = rdivide(t, a)
            assert(isnumeric(a), 'mpotensor:argerror', 'invalid input types.');
            t.tensors = t.tensors ./ a;
            t.scalars = t.scalars ./ a;
        end
        
        function t = mrdivide(t, a)
            assert(isnumeric(a), 'mpotensor:argerror', 'invalid input types.');
            t.tensors = t.tensors / a;
            t.scalars = t.scalars / a;
        end
        
        function v = applychannel(O, L, R, v)
            arguments
                O MpoTensor
                L MpsTensor
                R MpsTensor
                v
            end
            auxlegs_v = nspaces(v) - 3;
            auxlegs_l = nspaces(L) - 3;
            auxlegs_r = nspaces(R) - 3;
            auxlegs = auxlegs_v + auxlegs_l + auxlegs_r;
            
            v = contract(v, [1 3 5 (-(1:auxlegs_v) - 3 - auxlegs_l)], ...
                L, [-1 2 1 (-(1:auxlegs_l) - 3)], ...
                O, [2 -2 4 3], ...
                R, [5 4 -3 (-(1:auxlegs_r) - 3 - auxlegs_l - auxlegs_v)], ...
                'Rank', rank(v) + [0 auxlegs]);
        end
        
        function v = applympo(varargin)
            assert(nargin >= 3)
            v = varargin{end};
            R = varargin{end-1};
            L = varargin{end-2};
            N = nargin - 1;
            
            auxlegs_v = nspaces(v) - N;
            auxlegs_l = L.alegs;
            auxlegs_r = R.alegs;
            auxlegs = auxlegs_v + auxlegs_l + auxlegs_r;
            
            Oinds = arrayfun(@(x) [2*x-2 -x 2*x 2*x-1], 2:N-1, 'UniformOutput', false);
            O = [varargin(1:end-3); Oinds];
            v = contract(v, [1:2:2*N-1 (-(1:auxlegs_v) - N - auxlegs_l)], ...
                L, [-1 2 1 (-(1:auxlegs_l) - N)], ...
                O{:}, ...
                R, [2*N-1 2*N-2 -N (-(1:auxlegs_r) - N - auxlegs_l - auxlegs_v)], ...
                'Rank', rank(v) + [0 auxlegs]);
        end
        
        function O = rot90(O)
            O.tensors = tpermute(O.tensors, [2 3 4 1], [2 2]);
            O.scalars = permute(O.scalars, [2 3 4 1]);
        end
        
        function C = tensorprod(A, B, dimA, dimB, ca, cb, options)
            arguments
                A
                B
                dimA
                dimB
                ca = false
                cb = false
                options.NumDimensionsA = ndims(A)
            end
            
            assert(~isa(A, 'MpoTensor') || ~isa(B, 'MpoTensor'), 'mpotensor:tensorprod', ...
                'cannot contract two mpotensors.');
        
            if isa(A, 'MpoTensor')
                C = tensorprod(A.tensors, B, dimA, dimB, ca, cb);
                
                if nnz(A.scalars) > 0
                    assert(sum(dimA == 1 | dimA == 3, 'all') == 1, ...
                        'Cannot deduce output space unless leg 1 xor leg 3 is connected.');
                    assert(sum(dimA == 2 | dimA == 4, 'all') == 1, ...
                        'Cannot deduce output space unless leg 2 xor leg 4 is connected.');
                    
                    uncA = 1:nspaces(A); uncA(dimA) = [];
                    uncB = 1:nspaces(B); uncB(dimB) = [];
                    
                    A = reshape(permute(A.scalars, [uncA flip(dimA)]), ...
                        [prod(size(A, uncA)) prod(size(A, dimA))]);
                    B = reshape(tpermute(B, [dimB uncB], [length(dimB) length(uncB)]), ...
                        [prod(size(B, dimB)) prod(size(B, uncB))]);
                    C = C + reshape(sparse(A) * B, size(C));
                end
            else
                C = tensorprod(A, B.tensors, dimA, dimB, ca, cb);
                szC = size(C);
                if nnz(B.scalars) > 0
                    assert(sum(dimB == 1 | dimB == 3, 'all') == 1, ...
                        'Cannot deduce output space unless leg 1 xor leg 3 is connected.');
                    assert(sum(dimB == 2 | dimB == 4, 'all') == 1, ...
                        'Cannot deduce output space unless leg 2 xor leg 4 is connected.');
                    
                    uncA = 1:nspaces(A); uncA(dimA) = [];
                    uncB = 1:nspaces(B); uncB(dimB) = [];
                    
                    A = reshape(tpermute(A, [uncA dimA], [length(uncA) length(dimA)]), ...
                        [prod(size(A, uncA)) prod(size(A, dimA))]);
                    B = reshape(permute(B.scalars, [flip(dimB) uncB]), ...
                        [prod(size(B, dimB)) prod(size(B, uncB))]);
                    C = C + reshape(A * sparse(B), szC);
                end
            end
        end
        
        function O = ctranspose(O)
            O.tensors = tpermute(O.tensors', [4 1 2 3], [2 2]);
            O.scalars = conj(permute(O.scalars, [1 4 3 2]));
        end
        
        function O = transpose(O)
            O.tensors = tpermute(O.tensors, [3 4 1 2], [2 2]);
            O.scalars = permute(O.scalars, [3 4 1 2]);
        end
    end
    
    methods
        function s = space(O, i)
            s = space(O.tensors, i);
        end
        
        function s = pspace(O)
            s = space(O.tensors, 2);
        end
        
        function s = leftvspace(O, lvls)
            if nargin == 1
                s = space(O.tensors, 1);
            else
                s = space(O.tensors(lvls, :, :, :), 1);
            end
        end
        
        function s = rightvspace(O, lvls)
            if nargin == 1
                s = space(O.tensors, 3);
            else
                s = space(O.tensors(:, :, lvls, :), 3);
            end
        end
        
        function bool = istriu(O)
            sz = size(O, 1:4);
            sz1 = sz(1) * sz(2);
            sz2 = sz(3) * sz(4);
            bool = istriu(reshape(O.scalars, [sz1, sz2])) && ...
                istriu(reshape(O.tensors, [sz1, sz2]));
        end
        
        function bool = iseye(O)
            bool = nnz(O.tensors) == 0 && isequal(O.scalars, eye(size(O.scalars)));
        end
        
        function n = nnz(O)
            n = nnz(O.tensors) + nnz(O.scalars);
        end
    end
    
    methods 
        function t = subsref(t, s)
            assert(length(s) == 1, 'mpotensor:index', ...
                'only a single level of indexing allowed.');
            switch s.type
                case '.'
                    t = builtin('subsref', t, s);
                case '()'
                    t.scalars = subsref(t.scalars, s);
                    t.tensors = subsref(t.tensors, s);
                otherwise
                    error('mpotensor:index', 'invalid indexing expression');
            end
        end
        
        function t = subsasgn(t, s, v)
            assert(length(s) == 1, 'mpotensor:index', ...
                'only a single level of indexing allowed.');
            assert(strcmp(s.type, '()'), 'mpotensor:index', 'only () indexing allowed.');
            if isempty(t), t = MpoTensor(); end
            if isnumeric(v)
                t.scalars = subsasgn(t.scalars, s, v);
            elseif isa(v, 'MpoTensor')
                t.scalars = subsasgn(t.scalars, s, v.scalars);
                t.tensors = subsasgn(t.tensors, s, v.tensors);
            elseif isa(v, 'Tensor')
                t.tensors = subsasgn(t.tensors, s, v);
            end
        end
        
        function i = end(t, k, n)
            if n == 1
                i = prod(size(t));
                return
            end
            if n > ndims(t)
                i = 1;
            else
                i = size(t, k);
            end
        end
        
        function bools = eq(a, b)
            arguments
                a MpoTensor
                b MpoTensor
            end
            
            bools = a.scalars == b.scalars & a.tensors == b.tensors;
        end
        
        function I = find(O)
            I = union(find(O.tensors), find(O.scalars));
        end
        
        function varargout = size(t, varargin)
            if nargin > 1
                [varargout{1:nargout}] = size(t.scalars, varargin{:});
            else
                [varargout{1:nargout}] = size(t.scalars);
            end
        end
    end
    
    methods (Static)
        function O = zeros(m, n, o, p)
            if nargin == 1
                sz = m;
            elseif nargin == 4
                sz = [m n o p];
            else
                error('mpotensor:argerror', 'invalid amount of inputs.');
            end
            O = MpoTensor(SparseTensor([], [], sz), zeros(sz));
        end
    end
end
