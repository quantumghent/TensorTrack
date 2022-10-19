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
            T = transfermatrix(mpo, mps, mps, 'Type', 'LL');
            
%             for i = 1:size(mpo.O{1}, 1)
%                 tc.assertTrue(isapprox(apply(slice(T, i, 1:i-1), GL(1:i-1)), ...
%                     apply(slice(T, i, i), GL(i))));
%             end
            GL_ = apply(T, GL);
            for i = 1:numel(GL_)
                tc.verifyTrue(isapprox(GL_(i), GL(i)), ...
                    sprintf('GL(%d) disagrees', i));
            end
%             tc.assertTrue(isapprox(apply(T, GL), GL));
            
            [GR, lambdaR] = rightenvironment(mpo, mps, mps);
            tc.verifyTrue(isapprox(abs(lambdaR), abs(real(lambdaR)), 'RelTol', 1e-6), ...
                sprintf('lambda should be real. (%gi)', imag(lambdaR)));
            T = transfermatrix(mpo, mps, mps, 'Type', 'RR');
            GR_ = apply(T.', GR);
            for i = 1:numel(GR_)
                tc.verifyTrue(isapprox(GR_(i), GR(i)), ...
                    sprintf('GR(%d) disagrees', i));
            end
            
            tc.verifyTrue(isapprox(lambdaL, lambdaR), 'lambdas should be equal.')
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
        
        function test1dIsing(tc)            
            alg = Vumps('which', 'smallestreal', 'maxiter', 20, 'verbosity', Verbosity.iter);
            D = 16;
            mpo = InfJMpo.Ising(1, 1);
            mps = UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(D));
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            profile on
            mpo = InfJMpo.Ising(1, 1, 'Symmetry', 'Z2');
            mps = UniformMps.randnc(pspace(mpo), ...
                GradedSpace.new(Z2(0, 1), [D D] ./ 2, false));
            [mps2, lambda2] = fixedpoint(alg, mpo, mps);
            profile viewer
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

