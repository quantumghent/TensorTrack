classdef TestTensor < matlab.unittest.TestCase
    % TestTensor - Unit tests for tensors.
    
    properties
        tol = 1e-12
    end
    
    methods (TestClassSetup)
        function classSetup(tc)
            orig = Options.CacheEnabled;
            Options.CacheEnabled(false);
            tc.addTeardown(@Options.CacheEnabled, orig); 
        end
    end
    
    properties (TestParameter)
        spaces = struct(...
            'cartesian', CartesianSpace(3, [], 4, [], 5, [], 6, [], 7, []), ...
            'complex', ComplexSpace(3, false, 4, true, 5, false, 6, false, 7, true), ...
            'Z2', GradedSpace.new(Z2(0, 1), [1 1], false, Z2(0, 1), [1 2], true, ...
                Z2(0, 1), [3 2], true, Z2(0, 1), [2 3], false, Z2(0, 1), [2 5], false), ...
            'U1', GradedSpace.new(U1(0, 1, -1), [1 2 2], false, U1(0, 1, -1), [3 1 1], false, ...
                U1(0, 1, -1), [2 2 1], true, U1(0, 1, -1), [1, 2, 3], false, ...
                U1(0, 1, -1), [1 3 3], true), ...
            'SU2', GradedSpace.new(SU2(1, 2), [3 1], false, SU2(1, 3), [2 1], false, ...
                SU2(2, 3), [1 1], true, SU2(1, 2), [2 2], false, SU2(1, 2, 4), [1 1 1], true), ...
            'A4', GradedSpace.new(A4(1:2), [2 1], false, A4([1 4]), [1 2], false, ...
                A4(1:4), [2 1 3 1], true, A4(1:4), [2 2 1 2], false, A4(2:3), [1 2], true) ...
            )
    end
    
    methods (Test)
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
        
        
        function permute_via_inner(tc, spaces)
            rng(213);
            t1 = Tensor.rand(spaces, []);
            t2 = Tensor.rand(spaces, []);
            
            inner = dot(t1, t2);
            
            for i = 0:5
                ps = perms(1:nspaces(t1)).';
                for p = ps(:, randperm(size(ps, 2), min(size(ps, 2), 20)))
                    t3 = permute(t1, p.', [i 5-i]);
                    tc.assertTrue(all(dims(t1, p.') == dims(t3)), ...
                        'Incorrect size after permutation.');
                    tc.assertTrue(...
                        isapprox(norm(t1), norm(t3), 'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                        'Permute should preserve norms.')
                    
                    t4 = permute(t2, p.', [i 5-i]);
                    tc.assertTrue(all(dims(t2, p.') == dims(t4)), ...
                        'Incorrect size after permutation.');
                    tc.assertTrue(...
                        isapprox(dot(t3, t4), inner, 'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                        'Permute should preserve inner products.');
                end
            end
        end
        
        function permute_via_conversion(tc, spaces)
            t = Tensor.rand(spaces, []);
            a = double(t);
            rng(123);
            tc.assertTrue(all(dims(t) == size(a, 1:nspaces(t))));
            for k = 0:nspaces(t)
                ps = perms(1:nspaces(t)).';
                for p = ps(:, randperm(size(ps, 2), min(size(ps, 2), 10)))
                    t2 = permute(t, p.', [k nspaces(t)-k]);
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
            t1 = Tensor.randnc(spaces(1), spaces(2));
            t2 = Tensor.randnc(spaces(2), spaces(3));
            
            t3 = t1 * t2;
            tc.assertTrue(isapprox(double(t3), double(t1) * double(t2)));
            
            t1 = Tensor.randnc(spaces(1), spaces(2:3));
            t2 = Tensor.randnc(spaces(2:3), spaces(4));
            
            tc.assertTrue(isapprox(double(t1 * t2), ...
                tensorprod(double(t1), double(t2), [2 3], [2 1], 'NumDimensionsA', 3)));
            
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
            t1 = Tensor.randnc([], spaces(1:2));
            t2 = Tensor.randnc(spaces(1:2), []);
            
            tc.assertTrue(isapprox(tensorprod(t1, t2, [1 2], [2 1]), ...
                tensorprod(double(t1), double(t2), [1 2], [2 1]), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol));
        end
        
        function orthogonalize(tc, spaces)
            t = Tensor.randnc(spaces, []);
            
            %% Left orthogonalize
            p1 = [3 4 2];
            p2 = [1 5];
            tc.assumeTrue(spaces(3) * spaces(4) * spaces(2) >= spaces(1)' * spaces(5)')
            
            for alg = ["qr", "qrpos", "polar", "svd" "ql" "qlpos"]
                [Q, R] = leftorth(t, p1, p2, alg);
                
                assertTrue(tc, ...
                    isapprox(Q * R, permute(t, [p1 p2], [length(p1) length(p2)]), ...
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
            for alg = ["lq", "lqpos", "polar", "svd" "rq" "rqpos"]
                [L, Q] = rightorth(t, p1, p2, alg);
                
                assertTrue(tc, ...
                    isapprox(L * Q, permute(t, [p1 p2], [length(p1) length(p2)]), ...
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
            
            %% Left nullspace
            for alg = ["qr", "svd"]
                N = leftnull(t, [3 4 2], [1 5], alg);
                
                assertTrue(tc, norm(N' * permute(t, [3 4 2 1 5], [3 2])) < ...
                    100 * eps(norm(t)), ...
                    'N should be a left nullspace.');
                assertTrue(tc, isisometry(N, 'left', ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                    'N should be a left isometry.');
            end
            
            
            %% Right nullspace
            for alg = ["lq", "svd"]
                N = rightnull(t, [3 4], [2 1 5], alg);
                assertTrue(tc, norm(permute(t, [3 4 2 1 5], [2 3]) * N') < ...
                    100 * eps(norm(t)), ...
                    'N should be a right nullspace.');
                assertTrue(tc, isisometry(N, 'right', ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                    'N should be a right isometry.');
            end
        end
        
        function singularvalues(tc, spaces)
            t = Tensor.randnc(spaces, []);
            [U, S, V] = tsvd(t, [3 4 2], [1 5]);
            assertTrue(tc, isapprox(permute(t, [3 4 2 1 5], [3 2]), U * S * V), ...
                'USV should be a factorization.');
            assertTrue(tc, isisometry(U), ...
                'U should be an isometry.');
            assertTrue(tc, isisometry(V), ...
                'V should be an isometry.');
            
            %% truncation
            
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

