classdef TestInfMpo < matlab.unittest.TestCase
    % Unit tests for infinite matrix product operators.
    
    properties (TestParameter)
        mpo = struct(...
            'trivial', statmech2dIsing('Symmetry', 'Z1'), ...
            'Z2', statmech2dIsing('Symmetry', 'Z2'), ...
            'fermion', block(InfMpo.fDimer()) ...
            )
        mps = struct(...
            'trivial', UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(4)), ...
            'Z2', UniformMps.randnc(GradedSpace.new(Z2(0,1), [1 1], false), ...
                GradedSpace.new(Z2(0,1), [4 4], false)), ...
            'fermion', UniformMps.randnc(GradedSpace.new(fZ2(0, 1), [1 1], false), ...
                GradedSpace.new(fZ2(0, 1), [4 4], false)) ...
            )
    end
    
    methods (Test, ParameterCombination='sequential')
        function testEnvironments(tc)
            beta = 0.5;
            D = 16;
            
            %% testset 1: environments
            mpo = statmech2dIsing('beta', beta, 'Symmetry', 'Z2');
            mps = initialize_mps(mpo, GradedSpace.new(Z2(0, 1), D ./ [2 2], false));
            
            [GL, GR, lambda] = environments(mpo, mps);
            
            % left environment
            TL = transfermatrix(mpo, mps, mps, 'Type', 'LL');
            tc.assertTrue(isapprox(apply(TL, GL{1}), lambda * GL{1}), ...
                'left environment fixed point equation unfulfilled.');
            for i = 1:length(GL)
                tc.assertTrue(...
                    isapprox(apply(TL(i), GL{i}), lambda^(1/length(GL)) * GL{next(i, length(GL))}), ...
                    'left environment partial fixed point equation unfulfilled.');
            end
            
            % right environment
            TR = transfermatrix(mpo, mps, mps, 'Type', 'RR').';
            tc.assertTrue(isapprox(apply(TR, GR{1}), lambda * GR{1}), ...
                'right environment fixed point equation unfulfilled.');
            for i = length(GR):-1:1
                tc.assertTrue(...
                    isapprox(apply(TR(i), GR{i}), lambda^(1/length(GR)) * GR{prev(i, length(GR))}), ...
                    'right environment partial fixed point equation unfulfilled.');
            end
            
            % normalization
            for i = 1:length(GL)
                gl = multiplyright(MpsTensor(GL{i}), mps.C(i));
                gr = multiplyright(MpsTensor(GR{i}), mps.C(i)');
                tc.assertEqual(overlap(gl, gr), 1, 'AbsTol', 1e-10, ...
                    'environment normalization incorrect.');
            end
            
            %% testset 2: quasiparticle environments
            mpo = mpo / lambda;
            [GL, GR, lambda] = environments(mpo, mps, mps, GL, GR);
            tc.assertEqual(lambda, 1, 'AbsTol', 1e-10);
            
            for charge = Z2([0 1])
                for p = [0 pi 0.5]
                    qp = InfQP.randnc(mps, mps, p, charge);
                    
                    % left environments
                    GBL = leftquasienvironment(mpo, qp, GL, GR);
                    T_R = transfermatrix(mpo, qp, qp, 'Type', 'RL');
                    T_B = transfermatrix(mpo, qp, qp, 1, 'Type', 'BL');
                    if istrivial(qp)
                        C = qp.mpsleft.C(1);
                        FL_L = insert_onespace(multiplyright(MpsTensor(GL{1}), C), ...
                            nspaces(GL{1}) + 1, isdual(auxspace(qp, 1)));
                        FR_L = insert_onespace(multiplyright(MpsTensor(GR{1}), C'), ...
                            1, ~isdual(auxspace(qp, 1)));
                        
                        tc.verifyTrue(isapprox(...
                            apply_regularized(T_B, FL_L, FR_L, GL{1}) + ...
                            apply_regularized(T_R, FL_L, FR_L, GBL{1}), ...
                            exp(1i*p) * GBL{1}), ...
                            sprintf('left quasi environment fixed point equation unfulfilled for p=%e, c=%d', p, charge));
                    else
                        tc.verifyTrue(isapprox(apply(T_B, GL{1}) + apply(T_R, GBL{1}), ...
                            exp(1i*p) * GBL{1}), ...
                            sprintf('left quasi environment fixed point equation unfulfilled for p=%e, c=%d', p, charge));
                    end
                    
                    % right environments
                    GBR = rightquasienvironment(mpo, qp, GL, GR);
                    T_L = transfermatrix(mpo, qp, qp, 'Type', 'LR').';
                    T_B = transfermatrix(mpo, qp, qp, 'Type', 'BR').';
                    if istrivial(qp)
                        C = qp.mpsright.C(1);
                        FL_R = insert_onespace(multiplyleft(MpsTensor(GL{1}), C'), ...
                            1, ~isdual(auxspace(qp, 1)));
                        FR_R = insert_onespace(multiplyleft(MpsTensor(GR{1}), C), ...
                            nspaces(GR{1}) + 1, isdual(auxspace(qp, 1)));
                        
                        tc.verifyTrue(isapprox(...
                            apply_regularized(T_B, FR_R, FL_R, GR{1}) + ...
                            apply_regularized(T_L, FR_R, FL_R, GBR{1}), ...
                            exp(-1i*p) * GBR{1}), ...
                            sprintf('right quasi environment fixed point equation unfulfilled for p=%e, c=%d.', p, charge));
                    else
                        tc.verifyTrue(isapprox(apply(T_B, GR{1}) + apply(T_L, GBR{1}), ...
                            exp(-1i*p) * GBR{1}), ...
                            sprintf('right quasi environment fixed point equation unfulfilled for p=%e, c=%d.', p, charge));
                    end
                end
            end
        end
        
        function testDerivatives(tc, mpo, mps)
            [GL, GR] = environments(mpo, mps, mps);
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_AC)
                AC_ = mps.AC(i);
                [AC_.var, lambda] = eigsolve(H_AC{i}, mps.AC(i).var, 1, 'largestabs');
                AC_2 = apply(H_AC{i}, AC_);
                tc.assertTrue(isapprox(AC_2, lambda * AC_.var));
            end
            
            H_C = C_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_C)
                [C_, lambda] = eigsolve(H_C{i}, mps.C(i), 1, 'largestabs');
                tc.assertTrue(isapprox(apply(H_C{i}, C_), lambda * C_));
            end
        end
        
        function test2dIsing(tc)
            % compute exact solution
            beta = 0.9 * log(1 + sqrt(2)) / 2;
            
            f = statmech2dIsing_free_energy(beta);
            
            % compute fixedpoint
            D = 16;
            alg = Vumps('MaxIter', 10);
            mpo = statmech2dIsing('beta', beta, 'Symmetry', 'Z1');
            mps = UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(D));
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.assertEqual(-log(lambda) / beta, f, 'RelTol', 1e-5);
            
            % compute excitations
            qp = InfQP.new([], mps2, mps2);
            GBL = leftquasienvironment(mpo, qp);
            
            
            mps = UniformMps.randnc(GradedSpace.new(Z2(0, 1), [1 1], false), ...
                GradedSpace.new(Z2(0, 1), [D D] ./ 2, false));
            mpo = statmech2dIsing('beta', beta, 'Symmetry', 'Z2');
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.assertEqual(-log(lambda) / beta, f, 'RelTol', 1e-5);
            
            mps = [mps mps];
            mpo = [mpo mpo];
            
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.assertEqual(-log(lambda) / 2 / beta, f, 'RelTol', 1e-5);
            
            mps = [mps; mps];
            mpo = [mpo; mpo];
            
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.assertEqual(-log(lambda)/ 4 / beta, f, 'RelTol', 1e-5);
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

