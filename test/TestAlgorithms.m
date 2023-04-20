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
                
                alg_expand = Expand('which', alg.which, 'schmidtcut', 1e-7, 'finalize', Vomps('maxiter', 5));
                gs2 = changebonds(alg_expand, mpo, mps);
                tc.verifyGreaterThan(abs(fidelity(gs, gs2)), 0.9^unitcell)
            end
        end
        
        function test1dIsing_ordered(tc)
            J = 1;
            h = 1;
            e0 = quantum1dIsing_energy(J, h);
            d0 = @(k) quantum1dIsing_dispersion(k, 'J', J, 'h', h);
            D = 12;
            momenta = [0 pi 0.5];
            
            for L = 3
                %% No symmetry
                H = repmat(quantum1dIsing('h', h, 'J', J), 1, L);
                
                vspace = CartesianSpace.new(D, D+1, D+2);
                gs = initialize_mps(H, vspace);
                
                % Groundstate algorithms
                gs = fixedpoint(Vumps('which', 'smallestreal', 'maxiter', 5), ...
                    H, gs);
                tc.assertEqual(expectation_value(gs, H), e0 * L, 'RelTol', 1e-2);
                
                % Excitation algorithms
                for k = momenta
                    qp = InfQP.randnc(gs, gs, k);
                    [~, mu] = excitations(QPAnsatz(), H, qp);
                    tc.assertEqual(mu, d0(k/L), 'RelTol', 1e-3, ...
                        sprintf('qp failed at momentum %.2f', k));
                end
                
                %% Symmetry
                H = repmat(quantum1dIsing('h', h, 'J', J, 'Symmetry', 'Z2'), 1, L);
                
                vspace = GradedSpace.new(Z2(0, 1), [D D] / 2, false);
                gs = initialize_mps(H, vspace);
                
                % Groundstate algorithms
                gs = fixedpoint(Vumps('which', 'smallestreal', 'maxiter', 5), ...
                    H, gs);
                tc.assertEqual(expectation_value(gs, H), e0 * L, 'RelTol', 1e-2);
                
                % Excitation algorithms
                for k = momenta
                    qp = InfQP.randnc(gs, gs, k, Z2(1));
                    [~, mu] = excitations(QPAnsatz(), H, qp);
                    tc.assertEqual(mu, d0(k / L), 'RelTol', 1e-3, ...
                        sprintf('qp failed at momentum %.2f', k));
                end
            end
        end
    end
end

