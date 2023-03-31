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
            end
        end
        
        function test1dIsing(tc)
            E0 = -0.2525;
            momenta = [0 0.5 pi];
            Delta0 = [0.805833 0.813325 1.00029];
            
            for L = 1:3
                for symm = ["Z1", "Z2"]
                    H = quantum1dIsing('Symmetry', symm, 'h', 0.1);
                    H = repmat(H, 1, L);
                    if strcmp(symm, 'Z1')
                        vspace = CartesianSpace.new(12);
                    else
                        vspace = GradedSpace.new(Z2(0, 1), [6 6], false);
                    end
                    gs = initialize_mps(H, vspace);
                    
                    %% Groundstate algorithms
                    gs = fixedpoint(Vumps('which', 'smallestreal', 'maxiter', 5), ...
                        H, gs);
                    tc.assertEqual(expectation_value(gs, H), E0 * L, 'RelTol', 1e-2);
%                     gs = fixedpoint(IDmrg('which', 'smallestreal', 'maxiter', 5), ...
%                         H, gs);
%                     tc.verifyEqual(expectation_value(gs, H), E0 * L, 'RelTol', 1e-2);
                    if L > 1
                        gs2 = fixedpoint(IDmrg2('which', 'smallestreal', 'maxiter', 5), ...
                            H, gs);
                        tc.assertEqual(expectation_value(gs2, H), E0 * L, 'RelTol', 1e-2);
                        gs2 = fixedpoint(Vumps2('which', 'smallestreal', 'maxiter', 6), ...
                            H, gs);
                        tc.assertEqual(expectation_value(gs2, H), E0 * L, 'RelTol', 1e-2);
                    end
                    
                    %% Excitations
                    for i = 1:length(momenta)
                        qp = InfQP.randnc(gs, gs, momenta(i));
                        [qp, mu] = excitations(QPAnsatz(), H, qp);
                        tc.verifyEqual(mu, Delta0(i), 'RelTol', 1e-3);
                    end
                    
                end
            end
        end
    end
end

