classdef TestFiniteMpo < matlab.unittest.TestCase
    % Unit tests for finite matrix product operators.
    
    properties (TestParameter)
        mpo = struct('trivial0', ...
            FiniteMpo.randnc(ComplexSpace.new(2, false, 3, false), ComplexSpace.new(4, false)), ...
            'trivial1', ...
                FiniteMpo.randnc(ComplexSpace.new(2, false, 4, true, 2, true), ...
                ComplexSpace.new(2, true, 3, false)) ...
            )
        
    end
    
    methods (Test)
        function testFixedpoints(tc, mpo)
            [V, D] = eigsolve(mpo);
            tc.assertTrue(isapprox(mpo.apply(V), D * V));
            
            v2 = insert_onespace(V);
            v3 = apply(mpo, v2);
        end
        
        function testProperties(tc, mpo)
            mpo_tensor = Tensor(mpo);
            tc.assertTrue(isequal(mpo.domain, mpo_tensor.domain), ...
                'domain should remain fixed after conversion.');
            tc.assertTrue(isequal(mpo.codomain, mpo_tensor.codomain), ...
                'codomain should remain fixed after conversion.');
            
            v = initialize_fixedpoint(mpo);
            tc.assertTrue(isapprox(mpo.apply(v), mpo_tensor * v));
        end
        
        function testTransposes(tc, mpo)
            tc.assertTrue(isapprox(Tensor(mpo)', Tensor(mpo')), ...
                'ctranspose should not change mpo');
            tc.assertTrue(isequal(domain(mpo'), codomain(mpo)));
            tc.assertTrue(isequal(codomain(mpo'), domain(mpo)));
            tc.assertTrue(isapprox(Tensor(mpo).', Tensor(mpo.')), ...
                'transpose should not change mpo');
            tc.assertTrue(isequal(domain(mpo.'), codomain(mpo)'));
            tc.assertTrue(isequal(codomain(mpo.'), domain(mpo)'));
        end
        
        function test2dIsing(tc)
            L = 8;
            D = 64;
            alg = Dmrg('maxiter', 10, 'which', 'largestabs');
            
            mpo = statmech2dIsing('beta', 2, 'L', L);
            vspace_max = CartesianSpace.new(D);
            mps = initialize_mps(mpo, 'MaxVspace', vspace_max);
            
            [mps, envs, eta] = fixedpoint(alg, mpo, mps);
            
            mpo = statmech2dIsing('beta', 2, 'L', L, 'Symmetry', 'Z2');
            vspace_max = GradedSpace.new(Z2(0, 1), D ./ [2 2], false);
            mps = initialize_mps(mpo, 'MaxVspace', vspace_max);
            [mps, envs, eta] = fixedpoint(alg, mpo, mps);
        end
        
        function test1dIsing(tc)
            L = 8;
            D = 64;
            alg = Dmrg('miniter', 2, 'maxiter', 5, 'which', 'smallestreal');
            
            mpo = quantum1dIsing('L', L);
            vspace_max = CartesianSpace.new(D);
            mps = initialize_mps(mpo, 'MaxVspace', vspace_max);
            
            [mps, envs, eta] = fixedpoint(alg, mpo, mps);
            
            mpo = quantum1dIsing('L', L, 'Symmetry', 'Z2');
            vspace_max = GradedSpace.new(Z2(0, 1), D ./ [2 2], false);
            mps = initialize_mps(mpo, 'MaxVspace', vspace_max);
            [mps, envs, eta] = fixedpoint(alg, mpo, mps);
        end
    end
end

