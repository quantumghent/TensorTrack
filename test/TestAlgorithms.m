classdef TestAlgorithms < matlab.unittest.TestCase
    % Unit tests for algorithms
    
    properties (TestParameter)
        alg = {Vumps('which', 'smallestreal', 'maxiter', 5), ...
            ...IDmrg('which', 'smallestreal', 'maxiter', 5), ...
            Vumps2('which', 'smallestreal', 'maxiter', 6), ...
            IDmrg2('which', 'smallestreal', 'maxiter', 5) ...
            }
        symm = {'Z1', 'Z2'}
    end
    methods (Test)
        function test2dIsing(tc, alg, symm)
            alg.which = 'largestabs';
            for unitcell = 1:3
                if unitcell == 1 && (isa(alg, 'Vumps2') || isa(alg, 'IDmrg2'))
                    continue;
                end
                
                E0 = 2.5337 ^ unitcell;
                mpo1 = statmech2dIsing('Symmetry', symm);
                
                if strcmp(symm, 'Z1')
                    mps1 = initialize_mps(mpo1, CartesianSpace.new(12));
                else
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
                
                alg_expand = Expand('which', alg.which, 'schmidtcut', 1e-7, 'finalize', Vomps('maxiter', 10));
                gs2 = changebonds(alg_expand, mpo, mps);
                tc.verifyGreaterThan(abs(fidelity(gs, gs2)), 0.85^unitcell)
            end
        end
        
        function test1dIsing(tc, alg, symm)
            alg.which = 'smallestreal';
            for unitcell = 1:3
                if unitcell == 1 && (isa(alg, 'IDmrg2') || isa(alg, 'Vumps2'))
                    continue;
                end
                E0 = -1.273 * unitcell;
                
                mpo1 = quantum1dIsing('Symmetry', symm);
                if strcmp(symm, 'Z1')
                    mps1 = initialize_mps(mpo1, CartesianSpace.new(12));
                else
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
end

