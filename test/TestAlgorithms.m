classdef TestAlgorithms < matlab.unittest.TestCase
    % Unit tests for algorithms
    
    
    methods (Test)
        function test2dIsing(tc)
            mpo = InfMpo.Ising();
            mps = initialize_mps(mpo, CartesianSpace.new(12));
            
            alg = Vumps('which', 'largestabs', 'maxiter', 10);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533, 'RelTol', 1e-2);
            tc.verifyEqual(lambda, expectation_value(gs, mpo, gs), 'RelTol', 1e-2);
            
            alg = IDmrg('which', 'largestabs', 'maxiter', 10);
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
            
            mpo = InfMpo.Ising('Symmetry', 'Z2');
            mps = initialize_mps(mpo, GradedSpace.new(Z2(0, 1), [6 6], false));
            
            alg = Vumps('which', 'largestabs', 'maxiter', 5);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.533, 'RelTol', 1e-2);
            tc.verifyEqual(lambda, expectation_value(gs, mpo, gs), 'RelTol', 1e-2);
            
            alg = IDmrg('which', 'largestabs', 'maxiter',10);
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, 2.5337, 'RelTol', 1e-2);
            tc.verifyEqual(lambda, expectation_value(gs, mpo, gs), 'RelTol', 1e-2);
            
            mps2 = approximate(Vomps(), ...
                mpo, gs, initialize_mps(mpo, GradedSpace.new(Z2(0, 1), [12 12], false)));
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
            E0 = -1.273;
            
            %% no symmetry
            mpo1 = InfJMpo.Ising();
            mps1 = initialize_mps(mpo1, CartesianSpace.new(12));
            
            alg = Vumps('which', 'smallestreal', 'maxiter', 5);
            [gs1, lambda] = fixedpoint(alg, mpo1, mps1);
            tc.verifyEqual(lambda, E0, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs1, mpo1, gs1), E0, 'RelTol', 1e-2);
            
            alg = IDmrg('which', 'smallestreal', 'maxiter', 5);
            [gs1, lambda] = fixedpoint(alg, mpo1, mps1);
            tc.verifyEqual(expectation_value(gs1, mpo1, gs1), E0, 'RelTol', 1e-2);
            
            mpo2 = [mpo1 mpo1];
            mps2 = [mps1 mps1];
            
            alg = Vumps('which', 'smallestreal', 'maxiter', 5);
            [gs2, lambda] = fixedpoint(alg, mpo2, mps2);
            tc.verifyEqual(lambda, 2 * E0, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs2, mpo2, gs2), 2 * E0, 'RelTol', 1e-2);
            
            alg = Vumps2('which', 'smallestreal', 'maxiter', 5);
            [gs2, lambda] = fixedpoint(alg, mpo2, mps2);
            tc.verifyEqual(lambda, E0 * 2, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs2, mpo2, gs2), 2 * E0, 'RelTol', 1e-2);
            
            alg = IDmrg('which', 'smallestreal', 'maxiter', 5);
            [gs2, lambda] = fixedpoint(alg, mpo2, mps2);
            tc.verifyEqual(expectation_value(gs2, mpo2, gs2), 2 * E0, ...
                'RelTol', 1e-2);
            
            alg = IDmrg2('which', 'smallestreal', 'maxiter', 5);
            [gs2, lambda] = fixedpoint(alg, mpo2, mps2);
            tc.verifyEqual(expectation_value(gs2, mpo2, gs2), 2 * E0, 'RelTol', 1e-2);
            
            %% Z2 symmetry
            mpo1 = InfJMpo.Ising('Symmetry', 'Z2');
            mps1 = initialize_mps(mpo1, GradedSpace.new(Z2(0, 1), [6 6], false));
            
            alg = Vumps('which', 'smallestreal', 'maxiter', 5);
            [gs1, lambda] = fixedpoint(alg, mpo1, mps1);
            tc.verifyEqual(lambda, E0, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs1, mpo1, gs1), E0, 'RelTol', 1e-2);
            
            alg = IDmrg('which', 'smallestreal', 'maxiter', 5);
            [gs1, lambda] = fixedpoint(alg, mpo1, mps1);
            tc.verifyEqual(expectation_value(gs1, mpo1, gs1), E0, 'RelTol', 1e-2);
            
            mpo2 = [mpo1 mpo1];
            mps2 = [mps1 mps1];
            
            alg = Vumps('which', 'smallestreal', 'maxiter', 5);
            [gs2, lambda] = fixedpoint(alg, mpo2, mps2);
            tc.verifyEqual(lambda, 2 * E0, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs2, mpo2, gs2), 2 * E0, 'RelTol', 1e-2);
            
            alg = Vumps2('which', 'smallestreal', 'maxiter', 5);
            [gs2, lambda] = fixedpoint(alg, mpo2, mps2);
            tc.verifyEqual(lambda, E0 * 2, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs2, mpo2, gs2), 2 * E0, 'RelTol', 1e-2);
            
            alg = IDmrg('which', 'smallestreal', 'maxiter', 5);
            [gs2, lambda] = fixedpoint(alg, mpo2, mps2);
            tc.verifyEqual(expectation_value(gs2, mpo2, gs2), 2 * E0, ...
                'RelTol', 1e-2);
            
            alg = IDmrg2('which', 'smallestreal', 'maxiter', 5);
            [gs2, lambda] = fixedpoint(alg, mpo2, mps2);
            tc.verifyEqual(expectation_value(gs2, mpo2, gs2), 2 * E0, 'RelTol', 1e-2);
        end
    end
end

