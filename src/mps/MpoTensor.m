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
            auxlegs_l = L.alegs;
            auxlegs_r = R.alegs;
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
            
            Oinds = cellfun(@(x) [2*x-2 -x 2*x 2*x-1], 2:N-1, 'UniformOutput', false);
            O = [varargin(1:end-3); Oinds];
            v = contract(v, [1:2:2*N-1 ((-1:auxlegs_v) - N - auxlegs_l)], ...
                L, [-1 2 1 (-(1:auxlegs_l) - N)], ...
                O{:}, ...
                R, [2*N-1 2*N-2 -N (-(1:auxlegs_r) - N - auxlegs_l - auxlegs_v)], ...
                'Rank', rank(v) + [0 auxlegs]);
        end
        
        function y = applyleft(O, T, B, x)
            arguments
                O MpoTensor
                T MpsTensor
                B MpsTensor
                x
            end
            
            switch num2str([nspaces(x) nspaces(T) nspaces(B)])
                case '3  3  3'
                    if isdual(space(B, 2)), twist(B, 2); end
                    
                    y = contract(x, [4 2 1], O.tensors, [2 5 -2 3], ...
                        T, [1 3 -3], B, [4 5 -1], ...
                        'Rank', [2 1]);
                    scals = reshape(O.scalars, size(O, 1), size(O, 3));
                    for j = 1:size(O, 1)
                         cols = find(scals(j, :));
                         if isempty(cols), continue; end
                         y_ = contract(x(j), [3 -2 1], T, [1 2 -3], B, [3 2 -1], ...
                             'Rank', [2 1]);
                         for i = cols
                             y(i) = y(i) + scals(j, i) * y_;
                         end
                    end
                    y2 = contract(x, [4 2 1], O, [2 5 -2 3], ...
                        T, [1 3 -3], B, [4 5 -1], ...
                        'Rank', [2 1]);
                    try
                    assert(isapprox(y, y2))
                    catch
                        bla
                    end
                    
                otherwise
                    error('not implemented.');
            end
        end
        
        function y = applyright(O, T, B, x)
            arguments
                O MpoTensor
                T MpsTensor
                B MpsTensor
                x
            end
            
            switch num2str([nspaces(x) nspaces(T) nspaces(B)])
                case '3  3  3'
                    if isdual(space(B, 2)), twist(B, 2); end
                    
                    y = contract(x, [1 2 4], O.tensors, [-2 5 2 3], ...
                        T, [-1 3 1], B, [-3 5 4], ...
                        'Rank', [2 1]);
                    
                    scals = reshape(O.scalars, size(O, 1), size(O, 3));
                    for i = 1:size(O, 3)
                        rows = find(scals(:, i));
                        if isempty(rows), continue; end
                        y_ = contract(x(i), [1 -2 4], T, [-1 2 1], B, [-3 2 4], ...
                            'Rank', [2 1]);
                        for j = rows
                            y(j) = y(j) + scals(j, i) * y_;
                        end
                    end
                    
                otherwise
                    error('not implemented.');
            end
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
            
            assert(~isa(A, 'MpoTensor') || ~isa(B, 'MpoTensor'));
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
                    B = tpermute(B, [dimB uncB], [length(dimB) length(uncB)]);
                    C = C + A * B;
                end
            else
                C = tensorprod(A, B.tensors, dimA, dimB, ca, cb);
                
                if nnz(B.scalars) > 0
                    assert(sum(dimB == 1 | dimB == 3, 'all') == 1, ...
                        'Cannot deduce output space unless leg 1 xor leg 3 is connected.');
                    assert(sum(dimB == 2 | dimB == 4, 'all') == 1, ...
                        'Cannot deduce output space unless leg 2 xor leg 4 is connected.');
                    
                    uncA = 1:nspaces(A); uncA(dimA) = [];
                    uncB = 1:nspaces(B); uncB(dimB) = [];
                    
                    A = tpermute(A, [uncA flip(dimA)], [length(uncA) length(dimA)]);
                    B = reshape(permute(B.scalars, [flip(dimB) uncB]), ...
                        [prod(size(B, dimB)) prod(size(B, uncB))]);
                    C = C + A * B;
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
%         
%         function C = contract(tensors, indices, kwargs)
%             arguments (Repeating)
%                 tensors
%                 indices (1, :) {mustBeInteger}
%             end
%             
%             arguments
%                 kwargs.Conj (1, :) logical = false(size(tensors))
%                 kwargs.Rank = []
%                 kwargs.Debug = false
%             end
%             
%             % replace mpo with tensor + scalar contribution
%             i = find(cellfun(@(x) isa(x, 'MpoTensor'), tensors), 1);
%             
%             args1 = [tensors; indices];
%             args1{1, i} = args1{1, i}.tensors;
%             conj1 = kwargs.Conj;
%             
%             args2 = [tensors([1:i-1 i+1:end]); indices([1:i-1 i+1:end])];
%             conj2 = kwargs.Conj([1:i-1 i+1:end]);
%             
%             C = contract(args1{:}, ...
%                 'Conj', conj1, 'Rank', kwargs.Rank, 'Debug', kwargs.Debug) + ...
%                 contract(args2{:}, ...
%                 'Conj', conj2, 'Rank', kwargs.Rank, 'Debug', kwargs.Debug);
%             
%             
%             for ii = 1:length(tensors)
%                 if isa(tensors{ii}, 'MpoTensor') && strcmp(tensors{ii}.info, 'dense')
%                     tensors{ii} = tensors{ii}.var;
%                 end
%             end
%             
%             % replace scalar mpo's and correct indices
%             scalar = 1;
%             mpoInds = cellfun(@(x) isa(x, 'MpoTensor'), tensors);
%             
%             % collect scalars
%             for ii = 1:length(tensors)
%                 if mpoInds(ii)
%                     if kwargs.Conj(ii)
%                         scalar = scalar * conj(tensors{ii}.var);
%                     else
%                         scalar = scalar * tensors{ii}.var;
%                     end
%                 end
%             end
%             
%             if scalar == 0
%                 warning('Endresult is 0, probably could have been more efficient.');
%             end
%             
%             % remove mpo's from tensorlist
%             tensors(mpoInds) = [];
%             kwargs.Conj(mpoInds) = [];
%             
%             % correct indices
%             toCorrect    = indices(mpoInds);
%             indices = indices(~mpoInds);
%             
%             for ii = 1:length(toCorrect)
%                 vertical = toCorrect{ii}([1 3]);
%                 horizont = toCorrect{ii}([2 4]);
%                 vertMax = max(vertical);
%                 vertMin = min(vertical);
%                 horMax  = max(horizont);
%                 horMin  = min(horizont);
%                 
%                 for jj = 1:length(indices)
%                     indices{jj}(indices{jj} == vertMax) = vertMin;
%                     indices{jj}(indices{jj} == horMax ) = horMin;
%                 end
%                 for jj = ii+1:length(toCorrect)
%                     toCorrect{jj}(toCorrect{jj} == vertMax) = vertMin;
%                     toCorrect{jj}(toCorrect{jj} == horMax ) = horMin;
%                 end
%             end
%             
%             
%             %% Perform contraction
%             args = [tensors; indices];
%             C = contract(args{:}, ...
%                 'Conj', kwargs.Conj, 'Rank', kwargs.Rank, 'Debug', kwargs.Debug);
% %             [output, med] = Contract(tensors, indices, endcenter, med);
%             if ~isequal(scalar, 1)
%                 C = C * scalar;
%             end
%         end
    end
    
    methods
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
            
            if isnumeric(v)
                t.scalars = subsasgn(t.scalars, s, v);
            elseif isa(v, 'MpoTensor')
                t.scalars = subsasgn(t.scalars, s, v.scalars);
                t.tensors = subsasgn(t.tensors, s, v.tensors);
            elseif isa(v, 'Tensor')
                t.tensors = subsasgn(t.tensors, s, v);
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

