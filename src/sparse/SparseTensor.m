classdef (InferiorClasses = {?Tensor}) SparseTensor < AbstractTensor
    % Class for multi-dimensional sparse objects.
    
    %#ok<*PROPLC>
    
    properties
        codomain
        domain
    end
    
    properties (Access = private)
        ind = []
        sz = []
        var (:, 1) Tensor = Tensor.empty(0, 1);
    end
    
    %% Constructors
    methods
        function t = SparseTensor(varargin)
            if nargin == 0 || (nargin == 1 && isempty(varargin{1}))
                return;
                
            elseif nargin == 1  % cast from existing object
                source = varargin{1};
                
                if isa(source, 'SparseTensor')
                    t.ind = source.ind;
                    t.sz = source.sz;
                    t.var = source.var;
                    
                    t.codomain = source.codomain;
                    t.domain = source.domain;
                    
                elseif isa(source, 'Tensor')
                    t.sz = ones(1, nspaces(source(1)));
                    t.sz(1:ndims(source)) = size(source);
                    
                    t.ind = ind2sub_(t.sz, 1:numel(source));
                    t.var = source(:);
                    
                    [t.codomain, t.domain] = deduce_spaces(t);
                else
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
                [t.codomain, t.domain] = deduce_spaces(t);
                
            elseif nargin == 3  % indices, values and size
                ind = varargin{1};
                if ~isempty(ind) && ~isempty(varargin{2})
                    var = reshape(varargin{2}, [], 1);
                    if isscalar(var), var = repmat(var, size(ind, 1), 1); end
                    assert(size(ind, 1) == size(var, 1), 'sparse:argerror', ...
                        'indices and values must be the same size.');
                    sz = reshape(varargin{3}, 1, []);
                    assert(isempty(ind) || size(ind, 2) == length(sz), 'sparse:argerror', ...
                        'number of indices does not match size vector.');
                    assert(isempty(ind) || all(max(ind, [], 1) <= sz), 'sparse:argerror', ...
                        'indices must not exceed size vector.');
                    t.var = var;
                else
                    sz = reshape(varargin{3}, 1, []);
                end
                t.ind = ind;
                t.sz = sz;
                [t.codomain, t.domain] = deduce_spaces(t);
                
            elseif nargin == 5 % indices, values, size, codomain, domain
                ind = varargin{1};
                if ~isempty(ind) && ~isempty(varargin{2})
                    var = reshape(varargin{2}, [], 1);
                    if isscalar(var), var = repmat(var, size(ind, 1), 1); end
                    assert(size(ind, 1) == size(var, 1), 'sparse:argerror', ...
                        'indices and values must be the same size.');
                    sz = reshape(varargin{3}, 1, []);
                    assert(isempty(ind) || size(ind, 2) == length(sz), 'sparse:argerror', ...
                        'number of indices does not match size vector.');
                    assert(isempty(ind) || all(max(ind, [], 1) <= sz), 'sparse:argerror', ...
                        'indices must not exceed size vector.');
                    t.var = var;
                else
                    sz = reshape(varargin{3}, 1, []);
                    ind = double.empty(0, size(sz, 2));
                end
                t.ind = ind;
                t.sz = sz;
                t.codomain = varargin{4};
                t.domain = varargin{5};
            end
        end
    end
    
    methods (Static)
        function t = new(f, codomain, domain, kwargs)
            arguments
                f
                codomain SumSpace
                domain SumSpace
                kwargs.Density = 1
            end
            sz = [nsubspaces(codomain) flip(nsubspaces(domain))];
            
            inds = sort(randperm(prod(sz), round(prod(sz) * kwargs.Density)));
            subs = ind2sub_(sz, inds);
            vars = Tensor.empty(0, 1);
            for i = length(inds):-1:1
                for j = length(codomain):-1:1
                    subcodomain(j) = subspaces(codomain(j), subs(i, j));
                end
                for j = length(domain):-1:1
                    subdomain(j) = subspaces(domain(j), subs(i, end + 1 - j));
                end
                vars(i) = Tensor.new(f, subcodomain, subdomain);
            end
            t = SparseTensor(subs, vars, sz, codomain, domain);
        end
        
        function t = rand(varargin)
            t = SparseTensor.new(@rand, varargin{:});
        end
        
        function t = zeros(codomain, domain, kwargs)
            arguments
                codomain
                domain
                kwargs.Density = 0
            end
            t = SparseTensor.new(@zeros, codomain, domain, 'Density', kwargs.Density);
        end
        
        function t = eye(codomain, domain)
            t = SparseTensor.zeros(codomain, domain);
            
            sz1 = size(t, 1:length(codomain));
            sz2 = size(t, length(codomain) + (1:length(domain)));
            n = prod(sz1);
            assert(n == prod(sz2));
            
            inds = sub2ind_([n n], [1:n; 1:n].');
            for i = flip(1:n)
                t.ind(i,:) = ind2sub_(t.sz, inds(i));
                [cod, dom] = slice(codomain, domain, t.ind(i,:));
                t.var(i) = Tensor.eye(cod, dom);
            end
        end
    end
    
    
    %% Utility
    methods
        function [codomain, domain] = deduce_spaces(t)
            spaces = cell(1, ndims(t));
            for i = 1:length(spaces)
                for j = flip(1:size(t, i))
                    idx = find(t.ind(:, i) == j, 1);
                    if isempty(idx)
                        error('sparse:argerror', ...
                            'Cannot deduce %dth space at index %d.', i, j);
                    end
                    spaces{i}(j) = space(t.var(idx), i);
                end
            end
            Nout = indout(t.var(1));
            codomain = SumSpace(spaces{1:Nout});
            domain = SumSpace(spaces{(Nout+1):end})';
        end
        
        function B = full(A)
            inds = ind2sub_(A.sz, 1:prod(A.sz));
            [lia, locb] = ismember(inds, A.ind, 'rows');
            B(lia) = A.var(locb(lia));
            if ~all(lia)
                s = arrayfun(@(i) subspaces(space(A, i)), 1:ndims(A), 'UniformOutput', false);
                r = rank(A);
                for i = find(~lia).'
                    allspace = arrayfun(@(j) s{j}(inds(i, j)), 1:length(s));
                    B(i) = Tensor.zeros(allspace(1:r(1)), allspace(r(1)+1:end)');
                end
            end
            
            B = reshape(B, A.sz);
        end
        
        function A = sparse(A)
        end
        
        function sp = space(t, inds)
            sp = [t.codomain t.domain'];
            if nargin > 1
                sp = sp(inds);
            end
        end
        
        function n = nspaces(t)
            n = length(t.domain) + length(t.codomain);
        end
        
        function n = ndims(A)
            n = length(A.sz);
        end
        
        function r = rank(t, i)
            r = [length(t.codomain) length(t.domain)];
            if nargin > 1
                r = r(i);
            end
        end
        
        function varargout = size(a, i)
            if nargin == 1
                sz = a.sz;
            else
                sz = ones(1, max(i));
                sz(1:length(a.sz)) = a.sz;
                sz = sz(i);
            end
            
            if nargout <= 1
                varargout = {sz};
            else
                varargout = num2cell(sz);
            end
        end
        
        function n = numel(t)
            n = prod(t.sz);
        end
        
        function disp(t)
            r = t.rank;
            nz = nnz(t);
            if nz == 0
                fprintf('all-zero rank (%d, %d) %s:\n', r(1), r(2), class(t));
            else
                fprintf('rank (%d, %d) %s with %d nonzeros:\n', r(1), r(2), class(t), nz);
            end
            s = space(t);
            for i = 1:length(s)
                fprintf('\t%d.\t', i);
                disp(s(i));
                fprintf('\b');
            end
            fprintf('\n');
            
            if nz == 0
                fprintf('all-zero %s of size %s\n', class(t), ...
                    dim2str(t.sz));
                return
            end
            
            fprintf('%s of size %s with %d nonzeros:\n', class(t), ...
                dim2str(t.sz), nz);
            
            spc = floor(log10(max(double(t.ind), [], 1))) + 1;
            fmt_subs = sprintf("%%%du,", spc(1:end-1));
            fmt_subs = sprintf("%s%%%du", fmt_subs, spc(end));
            fmt = sprintf("\t(%s)", fmt_subs);
            for i = 1:nz
                fprintf('%s\t\t', compose(fmt, t.ind(i,:)));
                disp(t.var(i));
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
        
        function bool = isscalar(t)
            bool = prod(t.sz) == 1;
        end
        
        function n = nnz(t)
            n = length(t.var);
        end
        
        function bools = eq(a, b)
            arguments
                a SparseTensor
                b SparseTensor
            end
            
            if isscalar(a) && isscalar(b)
                bools = (isempty(a.var) && isempty(b.var)) || ...
                    (~isempty(a.var) && ~isempty(b.var) && a.var == b.var);
                return
            end
            
            if isscalar(a)
                if nnz(a) == 0
                    bools = true(size(b));
                    if nnz(b) ~= 0
                        bools(sub2ind_(b.sz, b.ind)) = false;
                    end
                else
                    bools = false(size(b));
                    bools(sub2ind_(b.sz, b.ind)) = a.var == b.var;
                end
                return
            end
            
            if isscalar(b)
                bools = b == a;
                return
            end
            
            assert(isequal(size(a), size(b)), 'sparse:dimerror', ...
                'input sizes incompatible');
            bools = true(size(a.inds));
            [inds, ia, ib] = intersect(a.ind, b.ind, 'rows');
            
            bools(sub2ind_(a.sz, a.ind)) = false;
            bools(sub2ind_(b.sz, b.ind)) = false;
            bools(sub2ind_(a.sz, inds)) = a.var(ia) == b.var(ib);
        end
    end
    
    %% Linear Algebra
    methods
        function a = conj(a)
            if ~isempty(a.var)
                a.var = conj(a.var);
            end
            a.codomain = conj(a.codomain);
            a.domain = conj(a.domain);
        end
        
        function a = ctranspose(a)
            for i = 1:nnz(a)
                a.var(i) = a.var(i)';
            end
            a.ind = fliplr(a.ind);
            a.sz = fliplr(a.sz);
            [a.codomain, a.domain] = swapvars(a.codomain, a.domain);
        end
        
        function d = dot(a, b)
            [~, ia, ib] = intersect(a.ind, b.ind, 'rows');
            if isempty(ia), d = 0; return; end
            d = dot(a.var(ia), b.var(ib));
        end
        
        function a = minus(a, b)
            a = a + (-b);
        end
        
        function c = mtimes(a, b)
            if isscalar(a) || isscalar(b)
                c = a .* b;
                return
            end
            
            szA = size(a);
            szB = size(b);
            assert(length(szA) == 2 && length(szB) == 2, 'sparse:argerror', ...
                'mtimes only defined for matrices.');
            assert(szA(2) == szB(1), 'sparse:dimerror', ...
                'incompatible sizes for mtimes.');
            
            if isnumeric(a)
                c = SparseTensor.zeros(szA(1), szB(2));
                
                for i = 1:szA(1)
                    for j = 1:szB(2)
                        c = subsasgn(c, substruct('()', {i, j}), ...
                            sum(a(i, :).' .* ...
                            subsref(b, substruct('()', {':', j}))));
                    end
                end
                return
            end
            
            if isnumeric(b)
                c = SparseTensor.zeros(szA(1), szB(2));
                
                for i = 1:szA(1)
                    for j = 1:szB(2)
                        c = subsasgn(c, substruct('()', {i, j}), ...
                            sum(subsref(a, substruct('()', {i, ':'})) .* ...
                            b(:, j).'));
                    end
                end
                return
            end
            
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
            n = max(ndims(a), ndims(b));
            assert(isequal(size(a, 1:n), size(b, 1:n)), ...
                'sparse:dimerror', 'input dimensions incompatible.');
            
            if ~issparse(a)
                if nnz(b) > 0
                    idx = sub2ind_(b.sz, b.ind);
                    a(idx) = reshape(a(idx), [], 1) + b.var;
                end
                return
            end
            
            if ~issparse(b)
                a = b + a;
                return
            end
            
            if isempty(b.ind), return; end
            if isempty(a.ind), a = b; return; end
            
            [lia, locb] = ismember(b.ind, a.ind, 'rows');
            a.var(locb(lia)) = a.var(locb(lia)) + b.var(lia);
            a.var = [a.var; b.var(~lia)];
            a.ind = [a.ind; b.ind(~lia, :)];
        end
        
        function t = rdivide(t, a)
            assert(isnumeric(a), 'method not implemented.');
            if nnz(t) > 0
                t.var = rdivide(t.var, a);
            end
        end
        
        function B = sum(A, dim)
            arguments
                A
                dim = []
            end
            
            if isscalar(A)
                B = A;
                return
            end
            
            if isempty(dim), dim = find(size(A) ~= 1, 1); end
            
            if strcmp(dim, 'all')
                if nnz(A) == 0
                    for i = flip(1:length(A.codomain))
                        cod(i) = SumSpace(subspaces(A.codomain(i), 1));
                    end
                    for i = flip(1:length(A.domain))
                        dom(i) = SumSpace(subspaces(A.domain(i), 1));
                    end
                    B = SparseTensor.zeros(cod, dom);
                    return
                end
                
                B = sum(A.var, 'all');
                return
            end
            
            if isvector(A)
                if nnz(A) == 0
                    B = SparseTensor.zeros(1, 1);
                else
                    B = sum(A.var, 'all');
                end
                return
            end
            
            if ismatrix(A)
                if dim == 1
                    B = SparseTensor.zeros(1, size(A, 2));
                    n = nnz(A);
                    for i = 1:size(A, 2)
                        if n == 0, break; end
                        idx = A.ind(:, 2) == i;
                        if ~any(idx), continue; end
                        B(1, i) = sum(A.var(idx));
                        n = n - sum(idx);
                    end
                    return
                end
                
                if dim == 2
                    B = SparseTensor.zeros(size(A, 1), 1);
                    n = nnz(A);
                    for i = 1:size(A, 2)
                        if n == 0, break; end
                        idx = A.ind(:, 1) == i;
                        if ~any(idx), continue; end
                        B(i, 1) = sum(A.var(idx));
                        n = n - sum(idx);
                    end
                    return
                end
            end
            error('TBA');
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
            
            A = reshape(permute(A, [uncA dimA]), [prod(szA(uncA)), prod(szA(dimA))]);
            B = reshape(permute(B, [flip(dimB) uncB]), [prod(szB(dimB)), prod(szB(uncB))]);
            
            if isempty(uncA) && isempty(uncB)
                C = 0;
                if nnz(A) > 0 && nnz(B) > 0
                    for i = 1:size(A, 1)
                        for j = 1:size(B, 2)
                            for k = 1:size(A, 2)
                                Aind = all(A.ind == [i k], 2);
                                if ~any(Aind), continue; end
                                Bind = all(B.ind == [k j], 2);
                                if ~any(Bind), continue; end
                                
                                C = C + ...
                                    tensorprod(A.var(Aind), B.var(Bind), dimA, dimB, ...
                                    'NumDimensionsA', options.NumDimensionsA);
                            end
                        end
                    end
                end
            else
                Cvar = A.var.empty(0, 1);
                Cind = double.empty(0, length(uncA) + length(uncB));
                
                if nnz(A) > 0 && nnz(B) > 0
                    for i = 1:size(A, 1)
                        for j = 1:size(B, 2)
                            for k = 1:size(A, 2)
                                Aind = all(A.ind == [i k], 2);
                                if ~any(Aind), continue; end
                                Bind = all(B.ind == [k j], 2);
                                if ~any(Bind), continue; end
                                if ~isempty(Cind) && all(Cind(end,:) == [i j], 2)
                                    Cvar(end) = Cvar(end) + ...
                                        tensorprod(A.var(Aind), B.var(Bind), dimA, dimB, ...
                                        'NumDimensionsA', options.NumDimensionsA);
                                else
                                    Cvar(end+1) = ...
                                        tensorprod(A.var(Aind), B.var(Bind), dimA, dimB, ...
                                        'NumDimensionsA', options.NumDimensionsA);
                                    Cind = [Cind; [i j]];
                                end
                            end
                        end
                    end
                end
                
                Ccod = space(A, uncA);
                Cdom = space(B, uncB)';
                
                C = SparseTensor(Cind, Cvar, [size(A, 1) size(B, 2)], Ccod, Cdom);
                C = reshape(C, szC);
                if size(Cind, 1) == prod(szC), C = full(C); end
            end
        end
        
        function t = times(t1, t2)
            if isscalar(t1) && ~isscalar(t2)
                t1 = repmat(t1, size(t2));
            elseif isscalar(t2) && ~isscalar(t1)
                t2 = repmat(t2, size(t1));
            end
            
            nd = max(ndims(t1), ndims(t2));
            assert(isequal(size(t1, 1:nd), size(t2, 1:nd)));
            
            if isnumeric(t1)
                t = t2;
                idx = sub2ind_(t.sz, t.ind);
                t1 = t1(idx);
                idx2 = find(t1);
                if ~isempty(idx2)
                    t.var = t.var(idx2) .* reshape(full(t1(idx2)), [],1);
                    t.ind = t.ind(idx2, :);
                else
                    t.var = t.var(idx2);
                    t.ind = t.ind(idx2, :);
                end
                return
            end
            
            if isnumeric(t2)
                t = t2 .* t1;
                return
            end
            
            if ~issparse(t1)
                idx = sub2ind_(t2.sz, t2.ind);
                t = t2;
                t.var = t.var .* t1(idx);
                return
            end
            
            if ~issparse(t2)
                t = t2 .* t1;
                return
            end
            
            t = t1;
            [t.ind, ia, ib] = intersect(t1.ind, t2.ind, 'rows');
            t.var = t1.var(ia) .* t2.var(ib);
        end
        
        function t = tpermute(t, p, r)
            for i = 1:numel(t.var)
                t.var(i) = tpermute(t.var(i), p, r);
            end
            t = permute(t, p);
            sp = space(t, p);
            t.codomain = sp(1:r(1));
            t.domain = sp(r(1) + (1:r(2)))';
        end
        
        function t = repartition(t, r)
            arguments
                t
                r (1,2) = [nspaces(t) 0]
            end
            
            if nnz(t) > 0
                t.var = arrayfun(@(x) repartition(x, r), t.var);
            end
            
            sp = space(t);
            t.codomain = sp(1:r(1));
            t.domain = sp(r(1) + (1:r(2)))';
        end
        
        function t = twist(t, i, inv)
            arguments
                t
                i
                inv = false
            end
            if nnz(t) > 0
                t.var = twist(t.var, i, inv);
            end
        end
        
        function t = twistdual(t, i, inv)
            arguments
                t
                i
                inv = false
            end
            if nnz(t) > 0
                t.var = twistdual(t.var, i, inv);
            end
        end
        
        function a = uminus(a)
            if ~isempty(a.var), a.var = -a.var; end
        end
        
        function a = uplus(a)
        end
    end
    
    methods
        function bool = ismatrix(a)
            bool = ndims(a) == 2;
        end
        
        function bool = istriu(a)
            assert(ismatrix(a), 'sparse:matrix', 'istriu is only defined for matrices');
            bool = all(a.ind(:, 1) <= a.ind(:, 2));
        end
    end
    
    
    %% Indexing
    methods
        function t = cat(dim, t, varargin)
            for i = 1:length(varargin)
                t2 = sparse(varargin{i});
                N = max(ndims(t), ndims(t2));
                dimcheck = 1:N;
                dimcheck(dim) = [];
                assert(isequal(size(t, dimcheck), size(t2, dimcheck)), ...
                    'sparse:dimagree', 'incompatible sizes for concatenation.');
                newinds = t2.ind;
                newinds(:, dim) = newinds(:, dim) + size(t, dim);
                t.var = vertcat(t.var, t2.var);
                t.ind = vertcat(t.ind, newinds);
                t.sz(dim) = t.sz(dim) + size(t2, dim);
            end
        end
        
        function i = end(t, k, n)
            if n == 1
                i = prod(t.sz);
                return
            end
            
            assert(n == length(t.sz), 'sparse:index', 'invalid amount of indices.')
            i = t.sz(k);
        end
        
        function [I, J, V] = find(t, k, which)
            arguments
                t
                k = []
                which = 'first'
            end
            
            if isempty(t.ind)
                I = [];
                J = [];
                V = [];
                return
            end
            
            [inds, p] = sortrows(t.ind, width(t.ind):-1:1);
            
            if ~isempty(k)
                if strcmp(which, 'first')
                    inds = inds(1:k, :);
                    p = p(1:k);
                else
                    inds = inds(end:-1:end-k+1, :);
                    p = p(end:-1:end-k+1);
                end
            end
            
            if nargout < 2
                I = sub2ind_(t.sz, inds);
                return
            end
            
            subs = sub2sub([t.sz(1) prod(t.sz(2:end))], t.sz, inds);
            I = subs(:,1);
            J = subs(:,2);
            
            if nargout > 2
                V = t.var(p);
            end
        end
        
        function t = horzcat(varargin)
            t = cat(2, varargin{:});
        end
        
        function t = permute(t, p)
            if ~isempty(t.ind)
                t.ind = t.ind(:, p);
            end
            t.sz = t.sz(p);
        end
        
        function t = reshape(t, varargin)
            if nargin == 1
                sz = varargin{1};
            else
                hasempty = find(cellfun(@isempty, varargin));
                if isempty(hasempty)
                    sz = [varargin{:}];
                elseif isscalar(hasempty)
                    varargin{hasempty} = 1;
                    sz = [varargin{:}];
                    sz(hasempty) = round(numel(t) / prod(sz));
                else
                    error('Can only accept a single empty size index.');
                end
            end
            assert(prod(sz) == prod(t.sz), ...
                'sparse:argerror', 'To reshape the number of elements must not change.');
            idx = sub2ind_(t.sz, t.ind);
            t.ind = ind2sub_(sz, idx);
            t.sz  = sz;
        end
        
        function t = sortinds(t)
            % Sort the non-zero entries by their index.
            %
            % Arguments
            % ---------
            % t : :class:`SparseTensor`
            %
            % Returns
            % -------
            % t : :class:`SparseTensor`
            %   tensor with sorted elements.
            
            if isempty(t), return; end
            [t.ind, p] = sortrows(t.ind, width(t.ind):-1:1);
            t.var = t.var(p);
        end
        
        function varargout = subsref(t, s)
            switch s(1).type
                case '()'
                    n = size(s(1).subs, 2);
                    if n == 1 % linear indexing
                        I = ind2sub_(t.sz, s(1).subs{1});
                        s(1).subs = arrayfun(@(x) I(:,x), 1:width(I), 'UniformOutput',false);
                    else
                        assert(n == size(t.sz, 2), 'sparse:index', ...
                            'number of indexing indices must match tensor size.');
                    end
                    f = true(size(t.ind, 1), 1);
                    newsz = zeros(1, size(s(1).subs, 2));
                    
                    for i = 1:size(s(1).subs, 2)
                        A = s(1).subs{i};
                        if strcmp(A, ':')
                            newsz(i) = t.size(i);
                            continue;
                        end
                        nA = length(A);
                        if nA ~= length(unique(A))
                            error("Repeated index in position %i",i);
                        end
                        if i > length(t.codomain)
                            t.domain(end-(i-length(t.codomain))+1) = ...
                                SumSpace(subspaces(t.domain(end-(i-length(t.codomain))+1), A));
                        else
                            t.codomain(i) = SumSpace(subspaces(t.codomain(i), A));
                        end
                        if ~isempty(t.ind)
                            B = t.ind(:, i);
                            P = false(max(max(A), max(B)) + 1, 1);
                            P(A + 1) = true;
                            f = and(f, P(B + 1));
                            [~, ~, temp] = unique([A(:); t.ind(f, i)], 'stable');
                            t.ind(f, i) = temp(nA+1:end);
                        end
                        newsz(i) = nA;
                    end
                    t.sz = newsz;
                    if ~isempty(t.ind)
                        t.ind = t.ind(f, :);
                        t.var = t.var(f);
                    end
                    if length(s) > 1
                        assert(isscalar(t))
                        t = subsref(t.var, s(2:end));
                    end
                    varargout{1} = t;
                case '.'
                    [varargout{1:nargout}] = builtin('subsref', t, s);
                otherwise
                    error('sparse:index', '{} indexing not defined');
            end
        end
        
        function n = numArgumentsFromSubscript(~, ~, ~)
            n = 1;
        end
        
        function t = subsasgn(t, s, v)
            assert(strcmp(s(1).type, '()'), 'sparse:index', 'only () indexing allowed');
            
            % Todo: check spaces when assigning
            
            if length(s(1).subs) == 1
                I = ind2sub_(t.sz, s(1).subs{1});
                s(1).subs = arrayfun(@(x) I(:,x), 1:width(I), 'UniformOutput',false);
            end
            
            assert(length(s(1).subs) == size(t.sz, 2), 'sparse:index', ...
                'number of indexing indices must match tensor size.');
            assert(all(t.sz >= cellfun(@max, s(1).subs)), 'sparse:bounds', ...
                'out of bounds assignment disallowed');
            
            subsize = zeros(1, size(s(1).subs, 2));
            for i = 1:length(subsize)
                if strcmp(s(1).subs{i}, ':')
                    s(1).subs{i} = 1:t.sz(i);
                elseif islogical(s(1).subs{i})
                    s(1).subs{i} = find(s(1).subs{i}).';
                end
                subsize(i) = length(s(1).subs{i});
            end
            if nnz(v) == 0, return; end
            if isscalar(v) && prod(subsize) > 1, v = repmat(v, subsize); end
            subs = combvec(s(1).subs{:}).';
            
            if isempty(t.ind)
                t.ind = subs;
                t.var = v(:);
            else
                for i = 1:size(subs, 1)
                    idx = find(all(t.ind == subs(i, :), 2));
                    if isempty(idx)
                        t.ind = [t.ind; subs];
                        t.var = [t.var; full(v(i))];
                    else
                        t.var(idx) = full(v(i));
                    end
                end
            end
            t.sz = max(t.sz, cellfun(@max, s(1).subs));
        end
        
        function t = vertcat(varargin)
            t = cat(1, varargin{:});
        end
    end
    
    
    
end

