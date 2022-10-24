classdef TestAlgorithms < matlab.unittest.TestCase
    % Unit tests for algorithms
    
    
    methods (Test)
        function test2dIsing(tc)
            mpo = InfMpo.Ising();
            mps = initialize_mps(mpo, CartesianSpace.new(12));
            
            alg = Vumps('which', 'largestabs', 'maxiter', 5);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533, 'RelTol', 1e-2);
            tc.verifyEqual(lambda, expectation_value(gs, mpo, gs), 'RelTol', 1e-2);
            
            alg = IDmrg('which', 'largestabs', 'maxiter',10);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.5337, 'RelTol', 1e-2);
            tc.verifyEqual(lambda, expectation_value(gs, mpo, gs), 'RelTol', 1e-2);
            
            mps2 = approximate(Vomps(), ...
                mpo, gs, initialize_mps(mpo, CartesianSpace.new(24)));
            tc.verifyEqual(fidelity(gs, mps2), 1, 'RelTol', 1e-2);
            
            mps = [mps mps];
            mpo = [mpo mpo];
            alg = IDmrg2('which', 'largestabs', 'maxiter', 10);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533^2, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs, mpo, gs), [2.533 2.533], ...
                'RelTol', 1e-2);
            
            alg = IDmrg('which', 'largestabs', 'maxiter', 10);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533^2, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs, mpo, gs), [2.533 2.533], ...
                'RelTol', 1e-2);
            
            alg = Vumps('which', 'largestabs', 'maxiter', 5);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533^2, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs, mpo, gs), [2.533 2.533], ...
                'RelTol', 1e-2);
        end
        
        function test1dIsing(tc)
            mpo = InfJMpo.Ising();
            mps = initialize_mps(mpo, CartesianSpace.new(12));
            
            alg = Vumps('which', 'smallestreal', 'maxiter', 5);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, -1.27, 'RelTol', 1e-2);
            tc.verifyEqual(lambda, expectation_value(gs, mpo, gs), 'RelTol', 1e-2);
            
            alg = IDmrg('which', 'smallestreal', 'maxiter',10);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.5337, 'RelTol', 1e-2);
            tc.verifyEqual(lambda, expectation_value(gs, mpo, gs), 'RelTol', 1e-2);
            
            mps2 = approximate(Vomps(), ...
                mpo, gs, initialize_mps(mpo, CartesianSpace.new(24)));
            tc.verifyEqual(fidelity(gs, mps2), 1, 'RelTol', 1e-2);
            
            mps = [mps mps];
            mpo = [mpo mpo];
            alg = IDmrg2('which', 'smallestreal', 'maxiter', 10);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533^2, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs, mpo, gs), [2.533 2.533], ...
                'RelTol', 1e-2);
            
            alg = IDmrg('which', 'smallestreal', 'maxiter', 10);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533^2, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs, mpo, gs), [2.533 2.533], ...
                'RelTol', 1e-2);
            
            alg = Vumps('which', 'smallestreal', 'maxiter', 5);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533^2, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs, mpo, gs), [2.533 2.533], ...
                'RelTol', 1e-2);
        end
    end
end

