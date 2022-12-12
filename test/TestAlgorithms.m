classdef TestAlgorithms < matlab.unittest.TestCase
    % Unit tests for algorithms
    
    properties (TestParameter)
        unitcell = {1, 2, 3, 4}
        alg = {Vumps('which', 'smallestreal', 'maxiter', 5), ...
            ...IDmrg('which', 'smallestreal', 'maxiter', 5), ...
            Vumps2('which', 'smallestreal', 'maxiter', 6), ...
            IDmrg2('which', 'smallestreal', 'maxiter', 5) ...
            }
        symm = {'Z1', 'Z2'}
    end
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
        
        function test1dIsing(tc, unitcell, alg, symm)
            
            tc.assumeTrue(unitcell > 1 || isa(alg, 'IDmrg') || isa(alg, 'Vumps'))
            E0 = -1.273 * unitcell;
            
            if strcmp(symm, 'Z1')
                mpo1 = InfJMpo.Ising();
                mps1 = initialize_mps(mpo1, CartesianSpace.new(12));
            else
                mpo1 = InfJMpo.Ising('Symmetry', 'Z2');
                mps1 = initialize_mps(mpo1, GradedSpace.new(Z2(0, 1), [6 6], false));
            end
            
            mpo = mpo1;
            mps = mps1;
            for i = 2:unitcell
                mpo = [mpo mpo1];
                mps = [mps mps1];
            end
            
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(lambda, E0, 'RelTol', 1e-2);
            tc.verifyEqual(expectation_value(gs, mpo), E0, 'RelTol', 1e-2);
        end
    end
end

