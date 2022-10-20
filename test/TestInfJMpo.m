classdef TestInfJMpo < matlab.unittest.TestCase
    % Unit tests for infinite matrix product operators.
    
    properties (TestParameter)
        mpo = struct(...
            'trivial', InfJMpo.Ising() ...
            )
        mps = struct(...
            'trivial', UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(4)) ...
            )
    end
    
    methods (Test, ParameterCombination='sequential')
        function testEnvironments(tc, mpo, mps)
            [GL, lambdaL] = leftenvironment(mpo, mps, mps);
            tc.verifyTrue(isapprox(abs(lambdaL), abs(real(lambdaL)), 'RelTol', 1e-6), ...
                sprintf('lambda should be real. (%-g i)', imag(lambdaL)));
            
            [GR, lambdaR] = rightenvironment(mpo, mps, mps);
            tc.verifyTrue(isapprox(abs(lambdaR), abs(real(lambdaR)), 'RelTol', 1e-6), ...
                sprintf('lambda should be real. (%gi)', imag(lambdaR)));
            
            tc.verifyTrue(isapprox(lambdaL, lambdaR), 'lambdas should be equal.');
        end
        
        function testDerivatives(tc, mpo, mps)
            [GL, GR] = environments(mpo, mps, mps);
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_AC)
                [AC_, lambda] = eigsolve(H_AC{i}, mps.AC(i), 1, 'largestabs');
                tc.assertTrue(isapprox(apply(H_AC{i}, AC_), lambda * AC_));
            end
            
            H_C = C_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_C)
                [C_, lambda] = eigsolve(H_C{i}, mps.C(i), 1, 'largestabs');
                tc.assertTrue(isapprox(apply(H_C{i}, C_), lambda * C_));
            end
        end
        
        function test1dIsing(tc)            
            alg = Vumps('which', 'smallestreal', 'maxiter', 5);
            D = 16;
            mpo = InfJMpo.Ising(1, 1);
            mps = UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(D));
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyTrue(isapprox(lambda, -1.27, 'RelTol', 1e-2))
            
            mpo = InfJMpo.Ising(1, 1, 'Symmetry', 'Z2');
            mps = UniformMps.randnc(pspace(mpo), ...
                GradedSpace.new(Z2(0, 1), [D D] ./ 2, false));
            [mps2, lambda2] = fixedpoint(alg, mpo, mps);
            tc.verifyTrue(isapprox(lambda, -1.27, 'RelTol', 1e-2))
            
            mpo = [mpo mpo];
            mps = [mps mps];
            [mps2, lambda2] = fixedpoint(alg, mpo, mps);
            tc.verifyTrue(isapprox(lambda2/2, -1.27, 'RelTol', 5e-2))
        end
    end
end

