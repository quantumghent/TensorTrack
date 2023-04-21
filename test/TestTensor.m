classdef TestTensor < matlab.unittest.TestCase
    % TestTensor - Unit tests for tensors.
    
    properties
        tol = 1e-12
    end
    
    methods (TestClassSetup)
        function classSetup(tc)
            orig = Options.CacheEnabled;
            Options.CacheEnabled(true);
            tc.addTeardown(@Options.CacheEnabled, orig); 
        end
    end
    
    properties (TestParameter)
        spaces = struct(...
            'cartesian', CartesianSpace(3, [], 4, [], 5, [], 6, [], 7, []), ...
            'complex', ComplexSpace(3, false, 4, true, 5, false, 6, false, 7, true), ...
            'Z2', GradedSpace.new(Z2(0, 1), [1 1], false, Z2(0, 1), [1 2], true, ...
                Z2(0, 1), [3 2], true, Z2(0, 1), [2 3], false, Z2(0, 1), [2 5], false), ...
            'fZ2', GradedSpace.new(fZ2(0, 1), [1 1], false, fZ2(0, 1), [1 2], true, ...
                fZ2(0, 1), [3 2], true, fZ2(0, 1), [2 3], false, fZ2(0, 1), [2 5], false), ...
            'U1', GradedSpace.new(U1(0, 1, -1), [1 2 2], false, U1(0, 1, -1), [3 1 1], false, ...
                U1(0, 1, -1), [2 2 1], true, U1(0, 1, -1), [1, 2, 3], false, ...
                U1(0, 1, -1), [1 3 3], true), ...
            'SU2', GradedSpace.new(SU2(1, 2), [3 1], false, SU2(1, 3), [2 1], false, ...
                SU2(2, 3), [1 1], true, SU2(1, 2), [2 2], false, SU2(1, 2, 4), [1 1 1], true), ...
            'A4', GradedSpace.new(A4(1:2), [2 1], false, A4(1:4), [2 1 2 1], false, ...
                A4([1 4]), [1 2], true, A4(1:4), [2 2 1 1], false, A4(2:4), [1 2 2], true), ...
            'U1xSU2', GradedSpace.new(...
                ProductCharge(U1(-1:1), SU2(2,1,2)), [1 2 1], false, ...
                ProductCharge(U1(-2:2), SU2(1,2,1,2,1)), [2 1 2 2 2], true, ...
                ProductCharge(U1(-1:1), SU2(2,1,2)), [1 2 2], false, ...
                ProductCharge(U1(-1:1), SU2(2,1,2)), [1 2 1], false, ...
                ProductCharge(U1(-1:1), SU2(2,1,2)), [1 2 2], false), ...
            'Hubbard', GradedSpace.new(...
                ProductCharge(U1(0, 1, 2), SU2(1, 2, 1), fZ2(0, 1, 0)), [1 1 1], false, ...
                ProductCharge(U1(-3:1), SU2(1, 2, 1, 2, 1), fZ2(0, 1, 0, 1, 0)), [1 1 3 1 1], false, ...
                ProductCharge(U1(0, 1), SU2(1, 2), fZ2(0, 1)), [1 1], true, ...
                ProductCharge(U1(0, 1, 2), SU2(1, 2, 1), fZ2(0, 1, 0)), [1 1 1], false, ...
                ProductCharge(U1(0, 1), SU2(2, 1), fZ2(1, 0)), [1 1], true) ...
            )
    end
    
    methods (Test)
        %% General properties
        function basic_linear_algebra(tc, spaces)
            t1 = Tensor.rand(spaces(1:3), spaces(4:5));
            
            tc.verifyTrue(isapprox(norm(t1)^2, dot(t1, t1), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol), 'Norm and dot incompatible.');
            
            a = rand();
            tc.verifyTrue(isapprox(norm(a .* t1), abs(a) * norm(t1), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                'Norm and scalar multiplication incompatible.');
            tc.verifyTrue(isapprox(t1 + t1, 2 .* t1, ...
                'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                '2*t and t+t incompatible.')
            
            tc.verifyTrue(isapprox(-t1, t1 .* (-1), 'AbsTol', tc.tol), ...
                '-t and t .* (-1) incompatible.');
            
            tc.verifyTrue(isapprox(t1 .* (1/a), t1 ./ a, 'AbsTol', tc.tol), ...
                '.* and ./ are incompatible');
            
            tc.verifyTrue(isapprox(norm(normalize(t1)), 1), ...
                'normalize should result in unit norm.');
            
            t2 = Tensor.rand(spaces(1:3), spaces(4:5));
            b = rand();
            tc.verifyTrue(isapprox(dot(b .* t2, a .* t1), conj(b) * a * dot(t2, t1), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol) && ...
                isapprox(dot(t2, t1), conj(dot(t1, t2)), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                'Dot should be sesquilinear.');
        end
        
        function transpose_via_conversion(tc, spaces)
            tc.assumeTrue(istwistless(braidingstyle(spaces)));
            
            t = Tensor.ones(spaces(1:3), spaces(4:5));
            tdagger = t';
            tc.assertTrue(isequal(t.domain, tdagger.codomain));
            tc.assertTrue(isequal(t.codomain, tdagger.domain));
            
            tc.assertEqual(flip(dims(tdagger)), dims(t));
            tc.assertEqual(conj(double(t)), double(conj(t)), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol, ...
                sprintf('conj(double(t)) should be double(conj(t)). (%e)', distance(conj(double(t)), double(conj(t)))));
        end
        
        function matrix_functions(tc, spaces)
            for i = 1:3
                t = Tensor.randnc(spaces(1:i), spaces(1:i));
                assertTrue(tc, isapprox(t*t, t^2, 'AbsTol', tc.tol, 'RelTol', tc.tol));
                assertTrue(tc, isapprox(t*t*t, t^3, 'AbsTol', tc.tol, 'RelTol', tc.tol));
                
                assertTrue(tc, isapprox((t^(1/2))^2, t, ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol));
                assertTrue(tc, isapprox((t^(1/3))^3, t, ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol));
                
                assertTrue(tc, isapprox(sqrtm(t)^2, t, 'AbsTol', tc.tol, 'RelTol', tc.tol));
                
                at = 0.01 * normalize(t);
                assertTrue(tc, isapprox(expm(at), ...
                    1 + at + at^2/2 + at^3/6 + at^4/24 + at^5/120, ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol));
                
                assertTrue(tc, isapprox(t * inv(t), inv(t) * t, ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol) && ...
                    isapprox(t * inv(t), t.eye(t.codomain, t.domain), ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol));
            end
        end
        
        
        %% Contractions
        function permute_via_inner(tc, spaces)
            rng(213);
            t1 = Tensor.rand(spaces, []);
            t2 = Tensor.rand(spaces, []);
            
            inner = dot(t1, t2);
            
            for i = 0:5
                ps = perms(1:nspaces(t1)).';
                for p = ps(:, randperm(size(ps, 2), min(size(ps, 2), 5)))
                    t3 = tpermute(t1, p.', [i 5-i]);
                    tc.assertTrue(all(dims(t1, p.') == dims(t3)), ...
                        'Incorrect size after permutation.');
                    tc.assertTrue(...
                        isapprox(norm(t1), norm(t3), 'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                        'Permute should preserve norms.')
                    
                    t4 = tpermute(t2, p.', [i 5-i]);
                    tc.assertTrue(all(dims(t2, p.') == dims(t4)), ...
                        'Incorrect size after permutation.');
                    tc.assertTrue(...
                        isapprox(dot(t3, t4), inner, 'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                        'Permute should preserve inner products.');
                end
            end
        end
        
        function permute_via_conversion(tc, spaces)
            tc.assumeTrue(istwistless(braidingstyle(spaces)));
            t = Tensor.rand(spaces, []);
            a = double(t);
            rng(123);
            tc.assertTrue(all(dims(t) == size(a, 1:nspaces(t))));
            for k = 0:nspaces(t)
                ps = perms(1:nspaces(t)).';
                for p = ps(:, randperm(size(ps, 2), min(size(ps, 2), 5)))
                    t2 = tpermute(t, p.', [k nspaces(t)-k]);
                    a2 = double(t2);
                    tc.assertTrue(all(dims(t2) == size(a2, 1:nspaces(t))));
                    tc.assertTrue(all(dims(t2) == size(a, p.')));
                    tc.assertTrue(...
                        isapprox(permute(a, p.'), a2, 'Abstol', tc.tol, 'RelTol', tc.tol), ...
                        'Permute should be compatible with conversion.');
                end
            end
        end
        
        function multiplication_via_conversion(tc, spaces)
            tc.assumeTrue(istwistless(braidingstyle(spaces)));
            t1 = Tensor.randnc(spaces(1), spaces(2));
            t2 = Tensor.randnc(spaces(2), spaces(3));
            
            t3 = t1 * t2;
            tc.assertTrue(isapprox(double(t3), double(t1) * double(t2)));
            
            t1 = Tensor.randnc(spaces(1), spaces(2:3));
            t2 = Tensor.randnc(spaces(2:3), spaces(4));
            
            tc.assertTrue(isapprox(double(t1 * t2), ...
                tensorprod(double(t1), double(t2), [2 3], [2 1], 'NumDimensionsA', 3)));
            
            l1 = contract(t1, [1 2 3], conj(t1), [1 2 3]);
            l2 = contract(double(t1), [1 2 3], conj(double(t1)), [1 2 3]);
            tc.assertTrue(isapprox(l1, l2));
            
            W1 = spaces(1:3);
            W2 = spaces(4:5);
            
            t1 = Tensor.randnc(W1, W1);
            t2 = Tensor.randnc(W1, W2);
            
            t1_array = double(t1);
            t2_array = double(t2);
            
            tc.assertTrue(isapprox(double(t1 * t2), ...
                tensorprod(t1_array, t2_array, [4 5 6], [3 2 1], ...
                'NumDimensionsA', 6), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol));
        end
        
        function tensorprod_via_conversion(tc, spaces)
            tc.assumeTrue(istwistless(braidingstyle(spaces)));
            t1 = Tensor.randnc([], spaces(1:2));
            t2 = Tensor.randnc(spaces(1:2), []);
            
            tc.assertTrue(isapprox(tensorprod(t1, t2, [1 2], [2 1]), ...
                tensorprod(double(t1), double(t2), [1 2], [2 1]), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol));
            
            A = Tensor.randnc(spaces(1), spaces(1:2));
            C = Tensor.randnc(spaces(1), spaces(1));
            AC = contract(A, [-1 -2 1], C, [1 -3]);
            
            tc.assertTrue(isapprox(double(AC), contract(double(A), [-1 -2 1], double(C), [1 -3])));
        end
        
        function tensortrace(tc, spaces)
            t1 = Tensor.randnc(spaces(1:3), spaces(1:3));
            t2 = contract(t1, [-1 -2 1 1 -3 -4], 'Rank', [2 2]);
            t3 = contract(t2, [-1 1 1 -2], 'Rank', [1 1]);
            t4 = contract(t1, [-1 1 2 2 1 -2], 'Rank', [1 1]);
            tc.assertTrue(isapprox(t3, t4, 'AbsTol', tc.tol, 'RelTol', tc.tol));
            
            % issue with fermionic traces:
            Nl = Tensor.randnc(spaces(1), spaces(1:2));
            Nr = Tensor.randnc(spaces(1:2), spaces(1));
            result1 = contract(Nl, [-1 1 2], Nr, [2 1 -2]);
            result2 = contract(contract(Nl, [-1 1 -4], Nr, [-2 1 -3]), [-1 1 -2 1]);
            tc.assertTrue(isapprox(result1, result2, 'AbsTol', tc.tol, 'RelTol', tc.tol));
            
            t5 = contract(t1, [1 2 3 3 2 1]);
            t6 = contract(t4, [1 1]);
            tc.assertTrue(isapprox(t5, t6, 'AbsTol', tc.tol, 'RelTol', tc.tol));
            if istwistless(braidingstyle(spaces))
                t7 = contract(double(t1), [1 2 3  3 2 1]);
                tc.assertTrue(isapprox(t6, t7, 'AbsTol', tc.tol, 'RelTol', tc.tol));
            end
        end
        
        function contract_order(tc, spaces)
            A = Tensor.randnc(spaces(1:2), spaces(1)');
            r = Tensor.randnc(spaces(1)', spaces(1)');
            
            args = {A, r, conj(A); [-1 2 1], [1 3], [-2 2 3]};
            
            r1 = contract(args{:}, 'Rank', rank(r));
            for p = perms(1:3)'
                args2 = args(:, p);
                r2 = contract(args2{:}, 'Rank', rank(r));
                tc.assertTrue(isapprox(r1, r2), ...
                    'Contraction order should leave result invariant.');
            end
        end
        
        
        %% Factorizations
        function orthogonalize(tc, spaces)
            t = Tensor.randnc(spaces, []);
            
            %% Left orthogonalize
            p1 = [3 4 2];
            p2 = [1 5];
            tc.assumeTrue(spaces(3) * spaces(4) * spaces(2) >= spaces(1)' * spaces(5)')
            
            for alg = ["qr", "qrpos", "polar", "svd", "ql", "qlpos"]
                [Q, R] = leftorth(t, p1, p2, alg);
                
                assertTrue(tc, ...
                    isapprox(Q * R, tpermute(t, [p1 p2], [length(p1) length(p2)]), ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                    sprintf('Q and R not a valid %s factorization.', alg));
                
                assertTrue(tc, isisometry(Q, 'left', ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                    'Q should be a left isometry.');
                
                switch alg
                    case 'polar'
                        tc.assertTrue(isposdef(R), 'R should be positive definite.');
                        
                    case {'qr', 'qrpos'}
                        tc.assertTrue(istriu(R), 'R should be upper triangular.');
                        
                    case {'ql', 'qlpos'}
                        tc.assertTrue(istril(R), 'R should be lower triangular.');
                        
                end
            end
            
            
            %% Right orthogonalize
            tc.assumeTrue(spaces(3) * spaces(4) <= spaces(1)' * spaces(2)' * spaces(5)');
            p1 = [3 4];
            p2 = [2 1 5];
            for alg = ["lq", "lqpos", "polar", "svd", "rq", "rqpos"]
                [L, Q] = rightorth(t, p1, p2, alg);
                
                assertTrue(tc, ...
                    isapprox(L * Q, tpermute(t, [p1 p2], [length(p1) length(p2)]), ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                    sprintf('Q and R not a valid %s factorization.', alg));
                
                assertTrue(tc, isisometry(Q, 'right', ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                    'Q should be a right isometry.');
                
                switch alg
                    case 'polar'
                        assertTrue(tc, isposdef(L), 'L should be positive definite.');
                        
                    case {'rq', 'rqpos'}
                        assertTrue(tc, istriu(L), 'L should be upper triangular.');
                        
                    case {'lq', 'lqpos'}
                        assertTrue(tc, istril(L) ,'L should be lower triangular.');
                        
                end
            end
        end
        
        function nullspace(tc, spaces)
            t = Tensor.randnc(spaces, []);
            tc.assumeTrue(spaces(3) * spaces(4) * spaces(2) >= spaces(1)' * spaces(5)', ...
                'tensor not full rank')
            %% Left nullspace
            for alg = ["qr", "svd"]
                N = leftnull(t, [3 4 2], [1 5], alg);
                dimN = dims(N, nspaces(N));
                dimW = prod(dims(t, [3 4 2]));
                dimV = prod(dims(t, [1 5]));
                tc.assertEqual(dimW, dimN + dimV, 'Nullspace should be full rank');
                assertTrue(tc, norm(N' * tpermute(t, [3 4 2 1 5], [3 2])) < ...
                    100 * eps(norm(t)), ...
                    'N should be a left nullspace.');
                assertTrue(tc, isisometry(N, 'left', ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                    'N should be a left isometry.');
            end
            
            
            %% Right nullspace
            tc.assumeTrue(spaces(3) * spaces(4) <= spaces(1)' * spaces(2)' * spaces(5)', ...
                'tensor not full rank');
            for alg = ["lq", "svd"]
                N = rightnull(t, [3 4], [2 1 5], alg);
                dimN = dims(N, 1);
                dimW = prod(dims(t, [3 4]));
                dimV = prod(dims(t, [2 1 5]));
                tc.assertEqual(dimV, dimW + dimN, 'Nullspace should be full rank');
                assertTrue(tc, norm(tpermute(t, [3 4 2 1 5], [2 3]) * N') < ...
                    100 * eps(norm(t)), ...
                    'N should be a right nullspace.');
                assertTrue(tc, isisometry(N, 'right', ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                    'N should be a right isometry.');
            end
        end
        
        function singularvalues(tc, spaces)
            t = normalize(Tensor.randc(spaces, []));
            [U, S, V] = tsvd(t, [3 4 2], [1 5]);
            assertTrue(tc, isapprox(tpermute(t, [3 4 2 1 5], [3 2]), U * S * V), ...
                'USV should be a factorization.');
            assertTrue(tc, isisometry(U), ...
                'U should be an isometry.');
            assertTrue(tc, isisometry(V), ...
                'V should be an isometry.');
            
            %% truncation
            d = max(cellfun(@(x) min(size(x, 1), size(x, 2)), matrixblocks(S)));
            [Utrunc, Strunc, Vtrunc, eta] = tsvd(t, [3 4 2], [1 5], 'TruncDim', d-1);
            assertTrue(tc, isapprox(norm(tpermute(t, [3 4 2 1 5], [3 2]) - ...
                Utrunc * Strunc * Vtrunc), eta, 'AbsTol', 1e-10, 'RelTol', 1e-6));
            assertTrue(tc, isisometry(U, 'left'));
            assertTrue(tc, isisometry(V, 'right'));
            d2 = max(cellfun(@(x) max(size(x, 1), size(x, 2)), matrixblocks(Strunc)));
            assertTrue(tc, d2 <= ceil(0.95*d));
            
            d = min(dims(S, 1:2));
            [Utrunc, Strunc, Vtrunc, eta] = tsvd(t, [3 4 2], [1 5], 'TruncTotalDim', ceil(0.9*d));
            assertTrue(tc, isapprox(norm(tpermute(t, [3 4 2 1 5], [3 2]) - ...
                Utrunc * Strunc * Vtrunc), eta, 'AbsTol', 1e-10, 'RelTol', 1e-6));
            assertTrue(tc, isisometry(U, 'left'));
            assertTrue(tc, isisometry(V, 'right'));
            d2 = max(dims(Strunc, 1:2));
            assertTrue(tc, d2 <= ceil(0.9*d));
            
            s = min(cellfun(@(x) min(diag(x), [], 'all'), matrixblocks(S)));
            [Utrunc, Strunc, Vtrunc, eta] = tsvd(t, [3 4 2], [1 5], 'TruncBelow', s * 1.2);
            assertTrue(tc, isapprox(norm(tpermute(t, [3 4 2 1 5], [3 2]) - ...
                Utrunc * Strunc * Vtrunc), eta, 'AbsTol', 1e-10, 'RelTol', 1e-6));
            assertTrue(tc, isisometry(U, 'left'));
            assertTrue(tc, isisometry(V, 'right'));
            s2 = min(cellfun(@(x) min(diag(x)), matrixblocks(Strunc)));
            assertTrue(tc, s * 1.2 <= s2);
            
        end
        
        function eigenvalues(tc, spaces)
            for i = 1:3
                t = Tensor.randnc(spaces(1:i), spaces(1:i));
                [V, D] = eig(t);
                tc.assertTrue(isapprox(t * V, V * D, 'AbsTol', tc.tol, 'RelTol', tc.tol));
                [V, D, W] = eig(t);
                tc.assertTrue(isapprox(W' * t, D * W', 'AbsTol', tc.tol, 'RelTol', tc.tol));
            end
        end
    end
    
end
