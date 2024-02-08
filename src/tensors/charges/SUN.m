classdef SUN < AbstractCharge
    % Irreducible representations of the special unitary group :math:`\mathrm{SU}(N)`.
    %
    % .. todo::
    %   Explain irrep labeling and give some references.
    %
    % Properties
    % ----------
    % I : (1, :) :class:`uint8`
    %   integer vector representation label
    
    properties
        I (1,:) uint8
    end
    
    methods
        function style = braidingstyle(~)
            style = BraidingStyle.Bosonic;
        end
        
        function a = conj(a)
            for i = 1:numel(a)
                a(i).I = a(i).I(1) - fliplr(a(i).I);
            end
        end
        
        function bools = eq(a, b)
            if isscalar(a)
                bools = reshape(all(a.I == vertcat(b.I), 2), size(b));
            elseif isscalar(b)
                bools = reshape(all(b.I == vertcat(a.I), 2), size(a));
            else
                assert(isequal(size(a), size(b)));
                bools = reshape(all(vertcat(a.I) == vertcat(b.I), 2), size(a));
            end
        end
        
        function F = Fsymbol(a, b, c, d, e, f)
            persistent cache
            if isempty(cache)
                cache = LRU;
            end
            
            if Options.CacheEnabled()
                key = GetMD5([a.I; b.I; c.I; d.I; e.I; f.I], 'Array', 'hex');
                F = get(cache, key);
                if isequal(F, [])
                    F = Fsymbol_(a, b, c, d, e, f);
                    cache = set(cache, key, F);
                end
            else
                F = Fsymbol_(a, b, c, d, e, f);
            end
        end
        
        function style = fusionstyle(~)
            style = FusionStyle.Generic;
        end
        
        function C = fusiontensor(a, b, c)
            persistent cache
            if isempty(cache)
                cache = LRU;
            end
            
            if Options.CacheEnabled()
                key = GetMD5([a.I; b.I; c.I], 'Array', 'hex');
                C = get(cache, key);
                if isempty(C)
                    C = clebschgordan(a, b, c);
                    cache = set(cache, key, C);
                end
            else
                C = clebschgordan(a, b, c);
            end
        end
        
        function c = mtimes(a, b)
            c = unique(directproduct(a, b));
        end
        
        function bools = ne(a, b)
            if isscalar(a)
                bools = reshape(any(a.I ~= vertcat(b.I), 2), size(b));
            elseif isscalar(b)
                bools = reshape(any(b.I ~= vertcat(a.I), 2), size(a));
            else
                assert(isequal(size(a), size(b)));
                bools = reshape(any(vertcat(a.I) ~= vertcat(b.I), 2), size(a));
            end
        end
        
        function N = Nsymbol(a, b, c)
            if numel(a) > 1
                N = arrayfun(@Nsymbol, a, b, c);
                return
            end
            N = sum(c == directproduct(a, b));
        end
        
        function e = one(a)
            e = SUN(zeros(1, length(a(1).I), 'uint8'));
        end
        
        function d = qdim(a)
            I = double(vertcat(a.I));
            d = 1;
            for k2 = 2:size(I, 2)
                for k1 = 1:k2-1
                    d = d .* (k2 - k1 + I(:, k1) - I(:, k2)) / (k2 - k1);
                end
            end
            d = reshape(double(d), size(a));
        end
        
        function R = Rsymbol(a, b, c, ~)
            persistent cache
            if isempty(cache)
                cache = LRU;
            end
            
            if Options.CacheEnabled()
                key = GetMD5([a.I; b.I; c.I], 'Array', 'hex');
                R = get(cache, key);
                if isempty(R)
                    R = Rsymbol_(a, b, c);
                    cache = set(cache, key, R);
                end
            else
                R = Rsymbol_(a, b, c);
            end
        end
        
        function [b, I] = sort(a)
            [~, I] = sortrows(vertcat(a.I));
            b = a(I);
        end
        
        function charge = SUN(varargin)
            if nargin == 0
                labels = [];
            else
                labels = vertcat(varargin{:});
            end
            for i = size(labels, 1):-1:1
                charge(i).I = labels(i, :);
            end
        end
    end
    
    methods %(Access = private)
        function J = annihilation(a)
            J = cellfun(@ctranspose, creation(a), 'UniformOutput', false);
        end
        
        function b = basis(a)
            N = length(a.I);
            if N == 1
                b = GtPattern(N, a.I);
                return
            end
            
            if N == 2
                Is = (a.I(2):a.I(1));
                b = GtPattern(2, [repmat(a.I(:), 1, length(Is)); Is]);
                return
            end     
            
            if N == 3
                M = [];
                I = a.I;
                is = I(2):I(1);
                js = I(3):I(2);
                for i = is
                    for j = js
                        M = [M [repmat([i; j], 1, i - j + 1); j:i]];
                    end
                end
                b = GtPattern(N, [repmat(I(:), 1, size(M, 2)); M]);
                return
            end
            
            allIs = arrayfun(@(i) a.I(i+1):a.I(i), 1:N-1, 'UniformOutput', false);
            M = [];
            for v = combvec(allIs{end:-1:1})
                b_ = basis(SUN(v(end:-1:1).'));
                M = [M b_.M];
            end
            b = GtPattern(N, [repmat(a.I(:), 1, size(M, 2)); M]);
            
        end
        
        function C = clebschgordan(a, b, c)
            DISPLAY = false;
            t_elapsed = tic;
            d1 = qdim(a);   d2 = qdim(b);   d3 = qdim(c);
            N = length(a.I);
            
            Jplus1 = creation(a);
            Jplus2 = creation(b);
            
            eqs = zeros(N-1, d1, d2, d1, d2);
            rows = false(N-1, d1, d2);
            cols = false(d1, d2);
            
            basis_a = basis(a);
            weigth_a = pWeight(basis_a);
            
            basis_b = basis(b);
            weigth_b = pWeight(basis_b);
            
            w3 = pWeight(highest_weight(c));
            wshift = (sum(a.I) + sum(b.I) - sum(c.I)) / N;
            
            for m1 = 1:length(basis_a)
                w1 = weigth_a(m1, :);
                w2 = w3 - w1 + wshift;
                
                m2list = all(w2 == weigth_b, 2);
                cols(m1, m2list) = true;
                for m2 = find(m2list).'
                    for l = 1:length(Jplus1)
                        m2_ = m2;
                        [m1_, ~, v] = find(Jplus1{l}(:, m1));
                        eqs(l, m1_, m2_, m1, m2) = eqs(l, m1_, m2_, m1, m2) + ...
                            reshape(v, 1, []);
                        rows(l, m1_, m2_) = true;
                        m1_ = m1;
                        [m2_, ~, v] = find(Jplus2{l}(:, m2));
                        eqs(l, m1_, m2_, m1, m2) = eqs(l, m1_, m2_, m1, m2) + ...
                            reshape(v, 1, 1, []);
                        rows(l, m1_, m2_) = true;
                    end
                end
            end
            
            eqs = reshape(eqs, (N-1) * d1 * d2, d1 * d2);
            reduced_eqs = eqs(rows(:), cols(:));
            
            solutions = null(reduced_eqs);
            N123 = size(solutions, 2);
            assert(N123 == Nsymbol(a, b, c));
            
            solutions = leftorth(rref(solutions', 1e-12)', 'qrpos');
            
            C = zeros(d1 * d2, d3, N123);
            C(cols, d3, :) = reshape(solutions, [], 1, N123);
            C = reshape(C, d1, d2, d3, N123);
            
            %% lower weight CGC
            basis_c = basis(c);
            weigth_c = pWeight(basis_c);
            
            Jmin1 = cellfun(@ctranspose, Jplus1, 'UniformOutput', false);
            Jmin2 = cellfun(@ctranspose, Jplus2, 'UniformOutput', false);
            Jmin3 = annihilation(c);
            
            w1list = unique(weigth_a, 'rows');
            w3list = sortrows(unique(weigth_c, 'rows'), 'descend');
            
            for alpha = 1:N123
                known = false(1, d3);
                known(d3) = true;
                
                for w3 = w3list(2:end, :).'
                    m3list = find(all(w3.' == weigth_c, 2));
                    jmax = length(m3list);
                    imax = 0;
                    for l = 1:N-1
                        w3_ = w3.'; w3_(l:l+1) = w3_(l:l+1) + [1 -1];
                        imax = imax + sum(all(w3_ == weigth_c, 2));
                    end
                    eqs = zeros(imax, jmax);
                    rhs = zeros(imax, d1, d2);
                    
                    i = 0;
                    rows = false(d1, d2);
                    for l = 1:length(Jmin1)
                        w3_ = w3.'; w3_(l:l+1) = w3_(l:l+1) + [1 -1];
                        m3_list = find(all(w3_ == weigth_c, 2));
                        for m3_ = m3_list.'
                            i = i + 1;
                            eqs(i, :) = reshape(Jmin3{l}(m3list, m3_), 1, jmax);
                            assert(known(m3_));
                            for w1_ = w1list.'
                                m1_list = find(all(w1_.' == weigth_a, 2)).';
                                w1 = w1_.'; w1(l:l+1) = w1(l:l+1) + [-1 1];
                                m1list = find(all(w1 == weigth_a, 2)).';
                                
                                w2_ = w3_ - w1_.' + wshift;
                                m2_list = find(all(w2_ == weigth_b, 2)).';
                                w2 = w2_; w2(l:l+1) = w2(l:l+1) + [-1 1];
                                m2list = find(all(w2 == weigth_b, 2)).';
                                for m2_ = m2_list
                                    for m1_ = m1_list
                                        CGCcoeff = C(m1_, m2_, m3_, alpha);
                                        for m1 = m1list
                                            m2 = m2_;
                                            Jm1coeff = Jmin1{l}(m1, m1_);
                                            rhs(i, m1, m2) = rhs(i, m1, m2) + ...
                                                Jm1coeff * CGCcoeff;
                                        end
                                        for m2 = m2list
                                            m1 = m1_;
                                            Jm2coeff = Jmin2{l}(m2, m2_);
                                            rhs(i, m1, m2) = rhs(i, m1, m2) + ...
                                                Jm2coeff * CGCcoeff;
                                        end
                                    end
                                end
                                rows(m1_list, m2list) = true;
                                rows(m1list, m2_list) = true;
                            end
                        end
                    end
                    ieqs = pinv(eqs);
                    [rows1, rows2] = find(rows);
                    
                    for j = 1:length(m3list)
                        m3 = m3list(j);
                        assert(~known(m3));
                        for k = 1:length(rows1)
                            C(rows1(k), rows2(k), m3, alpha) = ...
                                ieqs(j, :) * rhs(:, rows1(k), rows2(k));
                        end
                        known(m3) = true;
                    end
                end
            end
            
            t_elapsed = toc(t_elapsed);
            if DISPLAY
                fprintf('Generated SU(%d) Clebsch-Gordan Coefficient:\n', N);
                fprintf('\t(%s) x (%s) -> (%s)\n', num2str(a.I), num2str(b.I), num2str(c.I));
                fprintf('\ttime = %.2fs\n', t_elapsed);
                [bytes, unit] = memsize(C);
                fprintf('\tsize = %.2f%s\n', bytes, unit);
                fprintf('\n');
            end
        end
        
        function J = creation(a)
            N = length(a.I);
            
            J = cell(1, N-1);
            for i = 1:length(J)
                J{i} = zeros(qdim(a));
            end
            
            b = basis(a);
            for i = 1:length(b)
                m = b(i);
                
                for l = 1:N-1
                    n = [reshape(get(m, 1:l+1, l+1), 1, l+1) - (1:l+1) ...
                        reshape(get(m, 1:l-1, l-1), 1, l-1) - (1:l-1) - 1];
                    d = reshape(get(m, 1:l, l), 1, l) - (1:l);
                    mkl = get(m, 1:l, l);
                    for k = 1:l
                        numerator = -prod(n - mkl(k) + k);
                        if numerator == 0, continue; end
                        d2 = d;
                        d2(k) = [];
                        denominator = prod((d2 - mkl(k) + k) .* (d2 - mkl(k) + k - 1));
                        if denominator == 0, continue; end
                        m2 = set(m, k, l, mkl(k) + 1);
                        j = m2 == b;
                        J{l}(j, i) = sqrt(numerator / denominator);
                    end
                end
            end
            J = cellfun(@sparse, J, 'UniformOutput', false);
        end
        
        function c = directproduct(a, b)
            if qdim(a) > qdim(b)
                c = directproduct(b, a);
                return
            end
            
            newI = [];
            N = size(a.I, 2);
            
            for p = basis(a)
                t = b.I;
                bad_pattern = false;
                for k = 1:N
                    if bad_pattern, break; end
                    for l = N:-1:k
                        if bad_pattern, break; end
                        
                        b_kl = get(p, k, l);
                        if checkbounds(p, k, l-1)
                            b_kl = b_kl - get(p, k, l-1);
                        end
                        t(l) = t(l) + b_kl;
                        if l > 1
                            bad_pattern = t(l-1) < t(l);
                        end
                    end
                end
                
                if ~bad_pattern
                    newI = [newI; t];
                end
            end
            c = SUN(newI - newI(:, end));
        end
        
        function F = Fsymbol_(a, b, c, d, e, f)
            N1 = Nsymbol(a, b, e);
            N2 = Nsymbol(e, c, d);
            N3 = Nsymbol(b, c, f);
            N4 = Nsymbol(a, f, d);
            
            if N1 == 0 || N2 == 0 || N3 == 0 || N4 == 0
                F = zeros(N1, N2, N3, N4);
                return
            end
            
            A = fusiontensor(a, b, e);
            B = fusiontensor(e, c, d);
            B = reshape(B(:, :, 1, :), size(B, [1 2 4]));
            C = fusiontensor(b, c, f);
            D = fusiontensor(a, f, d);
            D = reshape(D(:, :, 1, :), size(D, [1 2 4]));
            
            F = contract(conj(D), [1 5 -4], conj(C), [2 4 5 -3], ...
                A, [1 2 3 -1], B, [3 4 -2]);
        end
        
        function p = highest_weight(a)
            N = length(a.I);
            if N == 1
                p = GtPattern(1, a.I);
                return
            end
            
            mask = fliplr(triu(true(N)));
            M = repmat(a.I(:), 1, N);
            p = GtPattern(N, reshape(M(mask), [], 1));
%             p_ = highest_weight(SUN(a.I(1:end-1)));
%             p = GtPattern(N, [a.I(:); p_.M]);
        end
        
        function s = normalize(s)
            for i = 1:numel(s)
                s(i).I = s(i).I - s(i).I(end);
            end
        end
        
        function R = Rsymbol_(a, b, c)
            N1 = Nsymbol(a, b, c);
            N2 = Nsymbol(b, a, c);
            
            if N1 == 0 || N2 == 0
                R = zeros(N1, N2);
                return
            end
            
            A = fusiontensor(a, b, c);
            B = fusiontensor(b, a, c);
            
            R = contract(conj(B(:, :, 1, :)), [1 2 -2], A(:, :, 1, :), [2 1 -1]);
        end
    end
end
