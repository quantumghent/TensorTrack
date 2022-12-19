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
        function test2dIsing(tc, alg, unitcell, symm)
            alg.which = 'largestabs';
            tc.assumeTrue(unitcell > 1 || isa(alg, 'IDmrg') || isa(alg, 'Vumps'))
            
            E0 = 2.5337 ^ unitcell;
            if strcmp(symm, 'Z1')
                mpo1 = InfMpo.Ising();
                mps1 = initialize_mps(mpo1, CartesianSpace.new(12));
            else
                mpo1 = InfMpo.Ising('Symmetry', 'Z2');
                mps1 = initialize_mps(mpo1, GradedSpace.new(Z2(0, 1), [6 6], false));
            end
            
            mpo = mpo1;
            mps = mps1;
            for i = 2:unitcell
                mpo = [mpo mpo1];
                mps = [mps mps1];
            end
            
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyEqual(expectation_value(gs, mpo, gs), E0, 'RelTol', 1e-2);
        end
        
        function test1dIsing(tc, unitcell, alg, symm)
            alg.which = 'smallestreal';
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
            tc.verifyEqual(expectation_value(gs, mpo), E0, 'RelTol', 1e-2);
        end
    end
end

