classdef TestSolvers < matlab.unittest.TestCase
    % TestLinsolve - Unit tests for linsolve.
    
    properties
        tol = 1e-10
    end
    
    properties (TestParameter)
        spaces = struct(...
            'cartesian', CartesianSpace.new([3 4]), ...
            'complex', ComplexSpace.new(3, false, 4, true, 2, false), ...
            'Z2', GradedSpace.new(Z2(0, 1), [1 1], false, Z2(0, 1), [1 2], true, ...
                Z2(0, 1), [3 2], true), ...
            'U1', GradedSpace.new(U1(0, 1, -1), [1 2 2], false, U1(0, 1, -1), [3 1 1], false), ...
            'SU2', GradedSpace.new(SU2(1, 2), [3 1], false, SU2(1, 3), [2 1], false) ...
            )
    end
    
    methods (Test)
        function linsolve(tc, spaces)
            A = Tensor.randnc(spaces, spaces);
            A = normalize((A + A') / 2);
            while cond(A) > 20
                A = Tensor.randnc(spaces, spaces);
                A = normalize((A + A') / 2);
            end
            
            xin = normalize(Tensor.randnc(spaces, []));
            b = A * xin;
            
            for alg = ["bicgstab", "bicgstabl", "gmres"]
                [x, flag, relres] = linsolve(A, b, 'Algorithm', alg);
                tc.assertTrue(isapprox(norm(A * x - b) / norm(b), relres, ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol));
                tc.assertTrue(flag == 0);

                f = @(x) A * x;
                [x2, flag, relres] = linsolve(f, b, 'Algorithm', alg);
                tc.assertTrue(isapprox(norm(f(x2) - b) / norm(b), relres, ...
                    'AbsTol', tc.tol, 'RelTol', tc.tol));
                tc.assertTrue(flag == 0);
            end
        end
        
        function eigsolve(tc, spaces)
            rng(123);
            A = Tensor.randc(spaces, spaces) - Tensor.eye(spaces, spaces) .* (1 + 1i);
            A = (A + A') ./ 2;
            
            x0 = Tensor.randc(spaces, []);
            
            d1 = eigsolve(A, x0, 1, 'IsSymmetric', true);
            [v, d, flag] = eigsolve(A, x0, 1, 'IsSymmetric', true);
            tc.verifyTrue(isapprox(d, d1));
            tc.verifyTrue(all(flag == 0));
            tc.verifyTrue(isapprox(A * v, v * d));
            
            [v, d, flag] = eigsolve(A, x0, 3, 'IsSymmetric', true);
            tc.verifyTrue(all(flag == 0));
            tc.verifyTrue(isapprox(A * v, v * d));
        end
    end
end

