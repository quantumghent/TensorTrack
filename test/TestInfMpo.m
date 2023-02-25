classdef TestInfMpo < matlab.unittest.TestCase
    % Unit tests for infinite matrix product operators.
    
    properties (TestParameter)
        mpo = struct(...
            'trivial', InfMpo.Ising(), ...
            'Z2', InfMpo.Ising('Symmetry', 'Z2'), ...
            'fermion', block(InfMpo.fDimer()) ...
            )
        mps = struct(...
            'trivial', UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(4)), ...
            'Z2', UniformMps.randnc(GradedSpace.new(Z2(0,1), [1 1], false), ...
                GradedSpace.new(Z2(0,1), [4 4], false)), ...
            'fermion', UniformMps.randnc(GradedSpace.new(fZ2(0, 1), [1 1], false), ...
                GradedSpace.new(fZ2(0, 1), [4 4], false)) ...
            )
    end
    
    methods (Test, ParameterCombination='sequential')
        function testEnvironments(tc, mpo, mps)
            [GL, lambdaL] = leftenvironment(mpo, mps, mps);
            T = transfermatrix(mpo, mps, mps, 'Type', 'LL');
            tc.assertTrue(isapprox(apply(T, GL{1}), lambdaL * GL{1}));
            
            [GR, lambdaR] = rightenvironment(mpo, mps, mps);
            T = transfermatrix(mpo, mps, mps, 'Type', 'RR');
            tc.assertTrue(isapprox(apply(T.', GR{1}), lambdaR * GR{1}));
            
            tc.assertTrue(isapprox(lambdaL, lambdaR));
        end
        
        function testDerivatives(tc, mpo, mps)
            [GL, GR] = environments(mpo, mps, mps);
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_AC)
                AC_ = mps.AC(i);
                [AC_.var, lambda] = eigsolve(H_AC{i}, mps.AC(i).var, 1, 'largestabs');
                AC_2 = apply(H_AC{i}, AC_);
                tc.assertTrue(isapprox(AC_2, lambda * AC_.var));
            end
            
            H_C = C_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_C)
                [C_, lambda] = eigsolve(H_C{i}, mps.C(i), 1, 'largestabs');
                tc.assertTrue(isapprox(apply(H_C{i}, C_), lambda * C_));
            end
        end
        
        function test2dIsing(tc)
            beta = 0.9 * log(1 + sqrt(2)) / 2;
            
            f = statmech2dIsing_free_energy(beta);
            
            D = 16;
            alg = Vumps('MaxIter', 10);
            mpo = statmech2DIsing('beta', beta, 'Symmetry', 'Z1');
            mps = UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(D));
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.assertEqual(-log(lambda) / beta, f, 'RelTol', 1e-5);
            
            mps = UniformMps.randnc(GradedSpace.new(Z2(0, 1), [1 1], false), ...
                GradedSpace.new(Z2(0, 1), [D D] ./ 2, false));
            mpo = statmech2DIsing('beta', beta, 'Symmetry', 'Z2');
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.assertEqual(-log(lambda) / beta, f, 'RelTol', 1e-5);
            
            mps = [mps mps];
            mpo = [mpo mpo];
            
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.assertEqual(-log(lambda) / 2 / beta, f, 'RelTol', 1e-5);
            
            mps = [mps; mps];
            mpo = [mpo; mpo];
            
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.assertEqual(-log(lambda)/ 4 / beta, f, 'RelTol', 1e-5);
        end
        
        function test2dfDimer(tc)
            D = 32;
            mpo = block(InfMpo.fDimer());
            mps = UniformMps.randnc(GradedSpace.new(fZ2(0, 1), [1 1], false), ...
                GradedSpace.new(fZ2(0, 1), [D D], false));
            [mps2, lambda] = fixedpoint(Vumps('tol', 1e-4, 'maxiter', 25), mpo, mps);
            tc.assertEqual(log(abs(lambda)) / 2, 0.29156, 'RelTol', 1e-4);
        end
    end
end

