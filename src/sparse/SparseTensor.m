classdef (InferiorClasses = {?Tensor}) SparseTensor
    % Class for multi-dimensional sparse objects.
    
    properties (Access = private)
        ind = []
        sz = []
        var (:, 1)
    end
    
    methods
        function t = SparseTensor(varargin)
            if nargin == 0 || (nargin == 1 && isempty(varargin{1}))
                return;
                
            elseif nargin == 1  % cast from existing object
                source = varargin{1};
                switch class(source)
                    case 'SparseTensor'
                        t.ind = source.ind;
                        t.sz = source.sz;
                        t.var = source.var;
                        
                    case 'Tensor'
                        t.sz = ones(1, nspaces(source(1)));
                        t.sz(1:ndims(source)) = size(source);
                        
                        t.ind = ind2sub_(t.sz, 1:numel(source));
                        t.var = source(:);
                        
                    otherwise
                        error('sparse:ArgError', 'Unknown syntax.');
                end
                
            elseif nargin == 2  % indices and values
                ind = varargin{1};
                var = reshape(varargin{2}, [], 1);
                if isscalar(var), var = repmat(var, size(ind, 1), 1); end
                assert(size(ind, 1) == size(var, 1), 'sparse:argerror', ...
                    'indices and values must be the same size.');
                t.ind = ind;
                t.var = var;
                t.sz = max(ind, [], 1);
                
            elseif nargin == 3  % indices, values and size
                ind = varargin{1};
                var = reshape(varargin{2}, [], 1);
                if isscalar(var), var = repmat(var, size(ind, 1), 1); end
                assert(size(ind, 1) == size(var, 1), 'sparse:argerror', ...
                    'indices and values must be the same size.');
                sz = reshape(varargin{3}, 1, []);
                assert(isempty(ind) || size(ind, 2) == length(sz), 'sparse:argerror', ...
                    'number of indices does not match size vector.');
                assert(isempty(ind) || all(max(ind, [], 1) <= sz), 'sparse:argerror', ...
                    'indices must not exceed size vector.');
                t.ind = ind;
                t.var = var;
                t.sz = sz;
                
            else
                error('sparse:argerror', 'unknown syntax.');
            end
        end
        
        function t = permute(t, p)
            if ~isempty(t.ind)
                t.ind = t.ind(:, p);
            end
            t.sz = t.sz(p);
        end
        
        function t = reshape(t, sz)
            assert(prod(sz) == prod(t.sz), ...
                'sparse:argerror', 'To reshape the number of elements must not change.');
            
            t.ind = sub2sub(sz, t.sz, t.ind);
            t.sz = sz;
        end
        
        function B = full(A)
            inds = ind2sub_(A.sz, 1:prod(A.sz));
            
            [lia, locb] = ismember(inds, A.ind, 'rows');
            B(lia) = A.var(locb(lia));
            
            if ~all(lia)
                s = arrayfun(@(i) space(A, i), 1:ndims(A), 'UniformOutput', false);
                r = rank(A.var(1));
                for i = find(~lia).'
                    allspace = arrayfun(@(j) s{j}(inds(i, j)), 1:length(s));
                    B(i) = Tensor.zeros(allspace(1:r(1)), allspace(r(1)+1:end)');
                end
            end
            B = reshape(B, A.sz);
        end
        
        function s = space(t, i)
            assert(isscalar(i), 'sparse:argerror', ...
                'Can only obtain spaces for single index.');
            for j = size(t, i):-1:1
                el = t.var(find(t.ind(:, i), 1));
                if isempty(el)
                    warning('cannot deduce space.');
                    continue;
                end
                s(j) = space(t.var(find(t.ind(:, i) == j, 1)), i);
            end
        end
        
        function n = ndims(A)
            n = length(A.sz);
        end
         
        function sz = size(a, i)
            if nargin == 1
                sz = a.sz;
                return
            end
            
            sz = ones(1, max(i));
            sz(1:length(a.sz)) = a.sz;
            sz = sz(i);
        end
        
        function disp(t)
            nz = nnz(t);
            if nz == 0
                fprintf('all-zero %s of size %s\n', class(t), ...
                    regexprep(mat2str(t.sz), {'\[', '\]', '\s+'}, {'', '', 'x'}));
                return
            end
            
            fprintf('%s of size %s with %d nonzeros:\n', class(t), ...
                regexprep(mat2str(t.sz), {'\[', '\]', '\s+'}, {'', '', 'x'}), nz);
            
            spc = floor(log10(max(double(t.ind), [], 1))) + 1;
            if numel(spc) == 1
                fmt = strcat("\t(%", num2str(spc(1)), "u)");
            else
                fmt = strcat("\t(%", num2str(spc(1)), "u,");
                for i = 2:numel(spc) - 1
                    fmt = strcat(fmt, "%", num2str(spc(i)), "u,");
                end
                fmt = strcat(fmt, "%", num2str(spc(end)), "u)");
            end
            
            for i = 1:nz
                fprintf('%s\t\t', compose(fmt, t.ind(i,:)));
                disp(t.var(i));
                fprintf('\n');
            end
        end
        
        function type = underlyingType(a)
            if isempty(a.var)
                type = 'double';
            else
                type = underlyingType(a.var);
            end
        end
        
        function bool = issparse(~)
            bool = true;
        end
        
        function n = nnz(t)
            n = length(t.var);
        end
    end
    
    %% Linear Algebra
    methods
        function a = conj(a)
            if ~isempty(a.var)
                a.var = conj(a.var);
            end
        end
        
        function d = dot(a, b)
            [~, ia, ib] = intersect(a.ind, b.ind, 'rows');
            if isempty(ia), d = 0; return; end
            d = dot(a.var(ia), b.var(ib));
        end
            
        function a = minus(a, b)
            assert(isequal(size(a), size(b)), ...
                'sparse:dimerror', 'input dimensions incompatible.');
            if isempty(b.ind), return; end
            if isempty(a.ind), a = -b; return; end
            
            [lia, locb] = ismember(b.ind, a.ind, 'rows');
            a.var(locb(lia)) = a.var(locb(lia)) - b.var(lia);
            if ~all(lia)
                a.var = [a.var; -b.var(~lia)];
                a.ind = [a.ind; b.ind(~lia, :)];
            end
        end
        
        function c = mtimes(a, b)
            szA = a.sz;
            szB = b.sz;
            assert(length(szA) == 2 && length(szB) == 2, 'sparse:argerror', ...
                'mtimes only defined for matrices.');
            assert(szA(2) == szB(1), 'sparse:dimerror', ...
                'incompatible sizes for mtimes.');
            
            cvar = [];
            cind = double.empty(0, 2);
            
            for k = 1:size(a, 2)
                rowlinds = a.ind(:, 2) == k;
                if ~any(rowlinds), continue; end
                
                collinds = b.ind(:, 1) == k;
                if ~any(collinds), continue; end
                
                rowinds = find(rowlinds);
                colinds = find(collinds);
                
                for i = rowinds.'
                    av = a.var(i);
                    ai = a.ind(i, 1);
                    for j = colinds.'
                        bv = b.var(j);
                        bj = b.ind(j, 2);
                        
                        
                        mask = all([ai bj] == cind, 2);
                        if any(mask)
                            cvar(mask) = cvar(mask) + av * bv;
                        else
                            cvar(end+1) = av * bv;
                            cind = [cind; ai bj];
                        end
                    end
                end
            end
            c = SparseTensor(cind, cvar, [szA(1) szB(2)]);
        end
        
        function n = norm(t, p)
            arguments
                t
                p = 'fro'
            end
            
            if isempty(t.var), n = 0; return; end
            n = norm(t.var);
        end
        
        function t = normalize(t)
            if isempty(t.var)
                warning('sparse:empty', 'cannot normalize an empty tensor.');
            end
            t = t .* (1 / norm(t));
        end
        
        function a = plus(a, b)
            assert(isequal(size(a), size(b)), ...
                'sparse:dimerror', 'input dimensions incompatible.');
            if isempty(b.ind), return; end
            if isempty(a.ind), a = b; return; end
            
            [lia, locb] = ismember(b.ind, a.ind, 'rows');
            a.var(locb(lia)) = a.var(locb(lia)) + b.var(lia);
            a.var = [a.var; b.var(~lia)];
            a.ind = [a.ind; b.ind(~lia, :)];
        end
        
        function C = tensorprod(A, B, dimA, dimB, ca, cb, options)
            arguments
                A SparseTensor
                B SparseTensor
                dimA
                dimB
                ca = false
                cb = false
                options.NumDimensionsA = ndims(A)
            end
            
            szA = size(A, 1:options.NumDimensionsA);
            szB = size(B, 1:max(ndims(B), max(dimB)));
            
            assert(length(dimA) == length(dimB) && all(szA(dimA) == szB(dimB)), ...
                'sparse:dimerror', 'incompatible contracted dimensions.');
            
            uncA = 1:length(szA); uncA(dimA) = [];
            uncB = 1:length(szB); uncB(dimB) = [];
            
            if isempty(uncA)
                if isempty(uncB)
                    szC = [1 1];
                elseif length(uncB) == 1
                    szC = [1 szB(uncB)];
                else
                    szC = szB(uncB);
                end
            elseif isempty(uncB)
                if length(uncA) == 1
                    szC = [szA(uncA) 1];
                else
                    szC = szA(uncA);
                end
            else
                szC = [szA(uncA) szB(uncB)];
            end
            
            Cvar = A.var.empty(0, 1);
            Cind = double.empty(0, length(uncA) + length(uncB));
            
            ks = unique([A.ind(:, dimA); B.ind(:, dimB)], 'rows');
            
            for k = ks.'
                rowlinds = all(A.ind(:, dimA) == k.', 2);
                if ~any(rowlinds), continue; end
                
                collinds = all(B.ind(:, dimB) == k.', 2);
                if ~any(collinds), continue; end
                
                rowinds = find(rowlinds);
                colinds = find(collinds);
                
                for i = rowinds.'
                    av = A.var(i);
                    ai = A.ind(i, uncA);
                    for j = colinds.'
                        bv = B.var(j);
                        bj = B.ind(j, uncB);
                        
                        mask = all([ai bj] == Cind, 2);
                        if any(mask)
                            Cvar(mask) = Cvar(mask) + ...
                                tensorprod(av, bv, dimA, dimB, ...
                                'NumDimensionsA', options.NumDimensionsA);
                        else
                            Cvar(end+1) = tensorprod(av, bv, dimA, dimB, ...
                                'NumDimensionsA', options.NumDimensionsA);
                            Cind = [Cind; ai bj];
                        end
                    end
                end
            end
            
            C = SparseTensor(Cind, Cvar, szC);
            if size(Cind, 1) == prod(szC), C = full(C); end
        end
        
        function t = times(t1, t2)
            if isnumeric(t1)
                if isempty(t2.var)
                    t = t2;
                    return
                end
                t = t2;
                t.var = t1 .* t.var;
                return
            end
            
            if isnumeric(t2)
                t = t2 .* t1;
                return
            end
            
            if isscalar(t1) && ~isscalar(t2)
                t1 = repmat(t1, size(t2));
            elseif isscalar(t2) && ~isscalar(t1)
                t2 = repmat(t2, size(t1));
            end
            
            assert(isequal(size(t1), size(t2)), 'sparse:dimerror', ...
                'incompatible input sizes.');
            
            if ~issparse(t1)
                if isempty(t2.var)
                    t = t2;
                    return
                end
                
                idx = sub2ind_(t2.sz, t2.ind);
                t = t2;
                t.var = t1(idx) .* t.var;
                return
            end
            
            if ~issparse(t2)
                if isempty(t1.var)
                    t = t1;
                    return
                end
                
                idx = sub2ind_(t1.sz, t1.ind);
                t = t1;
                t.var = t.var .* t2(idx);
                return
            end
            
            [inds, ia, ib] = intersect(t1.ind, t2.ind, 'rows');
            t = SparseTensor(inds, t1.var(ia) .* t2.var(ib), t1.sz);
        end
        
        function t = tpermute(t, p, r)
            for i = 1:numel(t.var)
                t.var(i) = tpermute(t.var(i), p, r);
            end
            t = permute(t, p);
        end
        
        function a = uminus(a)
            if ~isempty(a.var), a.var = -a.var; end
        end
        
        function a = uplus(a)
        end
    end
    
    
    %% Indexing
    methods
        function t = subsref(t, s)
            assert(length(s) == 1, 'sparse:index', 'only single level indexing allowed');
            assert(strcmp(s.type, '()'), 'sparse:index', 'only () indexing allowed');
            
            n = size(s.subs, 2);
            if n == 1 % linear indexing
                [s.subs{1:size(t.sz, 2)}] = ind2sub(t.sz, s.subs{1});
            else
                assert(n == size(t.sz, 2), 'sparse:index', ...
                    'number of indexing indices must match tensor size.');
            end
            f = true(size(t.ind, 1), 1);
            newsz = zeros(1, size(s.subs, 2));
            
            for i = 1:size(s.subs, 2)
                A = s.subs{i};
                if strcmp(A, ':')
                    newsz(i) = t.size(i);
                    continue;
                end
                nA = length(A);
                if nA ~= length(unique(A))
                    error("Repeated index in position %i",i);
                end
                if ~isempty(t.ind)
                    B = t.ind(:, i);
                    P = false(max(max(A), max(B)) + 1, 1);
                    P(A + 1) = true;
                    f = and(f, P(B + 1));
                    [~, ~, temp] = unique([A(:); t.ind(f, i)], 'stable');
                    t.ind(f, i) = temp(nA+1:end);
                    newsz(i) = nA;
                end
            end
            t.sz = newsz;
            if ~isempty(t.ind)
                t.ind = t.ind(f, :);
                t.var = t.var(f);
            end
        end
        
        function t = subsasgn(t, s, v)
            assert(strcmp(s(1).type, '()'), 'sparse:index', 'only () indexing allowed');
            assert(length(s(1).subs) == size(t.sz, 2), 'sparse:index', ...
                'number of indexing indices must match tensor size.');
            
            subsize = zeros(1, size(s(1).subs, 2));
            for i = 1:length(subsize)
                if strcmp(s(1).subs{i}, ':')
                    s(1).subs{i} = 1:t.sz(i);
                elseif islogical(s(1).subs{i})
                    s(1).subs{i} = find(s(1).subs{i}).';
                end
                subsize(i) = length(s(1).subs{i});
            end
            if isscalar(v), v = repmat(v, subsize); end
            subs = combvec(s(1).subs{:});
            for i = 1:size(subs, 1)
                idx = find(all(t.ind == subs(i, :), 2));
                if isempty(idx)
                    t.ind = [t.ind; subs];
                    t.var = [t.var; v(i)];
                else
                    t.var(idx) = v(i);
                end
            end
            t.sz = max(t.sz, cellfun(@max, s(1).subs));
        end
    end
end

