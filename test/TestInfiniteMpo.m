classdef TestInfiniteMpo < matlab.unittest.TestCase
    % Unit tests for infinite matrix product operators.
    
    properties (TestParameter)
        mpo = struct(...
            'trivial', InfiniteMpo.Ising(), ...
            'Z2', InfiniteMpo.Ising('Symmetry', 'Z2') ...
            )
        mps = struct(...
            'trivial', UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(4)), ...
            'Z2', UniformMps.randnc(GradedSpace.new(Z2(0,1), [1 1], false), ...
                GradedSpace.new(Z2(0,1), [4 4], false)) ...
            )
    end
    
    methods (Test, ParameterCombination='sequential')
        function testEnvironments(tc, mpo, mps)
            
            T = transfermatrix(mpo, mps);
            [GL, lambdaL] = eigsolve(T, [], 1, 'largestabs');
            tc.assertTrue(isapprox(apply(T, GL), lambdaL * GL));
            
            [GR, lambdaR] = eigsolve(T', [], 1, 'largestabs');
            tc.assertTrue(isapprox(apply(T', GR), lambdaR * GR));
            
            N = period(mps);
        end
        
        function testDerivatives(tc, mpo, mps)
            [GL, GR] = environments(mpo, mps, mps);
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            H_C = C_hamiltonian(mpo, mps, GL, GR);
            
            [AC_, lambda] = eigsolve(H_AC, mps.AC, 1, 'largestabs');
            tc.assertTrue(isapprox(apply(H_AC, AC_), lambda * AC_));
            
            [C_, lambda] = eigsolve(H_C, mps.C, 1, 'largestabs');
            tc.assertTrue(isapprox(apply(H_C, C_), lambda * C_));
        end
        
        function test2dIsing(tc)
            beta = 0.95 * log(1 + sqrt(2)) / 2;
            theta = 0:1e-6:pi/2;
            x = 2 * sinh(2 * beta) / cosh(2 * beta)^2;
            freeEnergyExact = -1 / beta * (log(2 * cosh(2 * beta)) + 1 / pi * ...
                trapz(theta, log(1/2 * (1 + sqrt(1 - x^2 * sin(theta).^2)))));
            
            D = 64;
            
            mpo = InfiniteMpo.Ising(beta);
            mps = UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(D));
            [mps2, lambda] = fixedpoint(Vumps(), mpo, mps);
            tc.assertEqual(-log(lambda) / beta, freeEnergyExact, 'RelTol', 1e-5);
            
            mps = UniformMps.randnc(GradedSpace.new(Z2(0, 1), [1 1], false), ...
                GradedSpace.new(Z2(0, 1), [D D] ./ 2, false));
            mpo = InfiniteMpo.Ising(beta, 'Symmetry', 'Z2');
            [mps2, lambda] = fixedpoint(Vumps(), mpo, mps);
            tc.assertEqual(-log(lambda) / beta, freeEnergyExact, 'RelTol', 1e-5);
        end
        
        function test2dfDimer(tc)
            D = 32;
            mpo = block(InfiniteMpo.fDimer());
            mps = UniformMps.randnc(GradedSpace.new(fZ2(0, 1), [1 1], false), ...
                GradedSpace.new(fZ2(0, 1), [D D], false));
            [mps2, lambda] = fixedpoint(Vumps('tol', 1e-4, 'maxiter', 25), mpo, mps);
            tc.assertEqual(log(abs(lambda)) / 2, 0.29156, 'RelTol', 1e-4);
        end
    end
end

