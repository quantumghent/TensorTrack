classdef TestInfJMpo < matlab.unittest.TestCase
    % Unit tests for infinite matrix product operators.
    
    properties (TestParameter)
        mpo = struct(...
            'trivial', quantum1dIsing() ...
            )
        mps = struct(...
            'trivial', UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(4)) ...
            )
    end
    
    methods (Test, ParameterCombination='sequential')
        function testEnvironments(tc)
            D = 16;
            mpo = quantum1dIsing('Symmetry', 'Z2');
            mps = initialize_mps(mpo, GradedSpace.new(Z2(0, 1), D ./ [2 2], false));
            
            [GL, GR, lambda] = environments(mpo, mps);
            
            % left environment
            fp_left  = repartition(insert_onespace(fixedpoint(mps, 'l_LL'), ...
                2, ~isdual(leftvspace(mpo, 1))), rank(GL{1}));
            fp_right = insert_onespace(fixedpoint(mps, 'r_LL'), ...
                2, ~isdual(rightvspace(mpo, 1)));
            
            TL = transfermatrix(mpo, mps, mps, 'Type', 'LL');
            GL_2 = apply(TL, GL{1});
            lambda2 = overlap(GL_2(end), fp_right);
            GL_2(end) = GL_2(end) - lambda2 * fp_left;
            tc.assertTrue(isapprox(GL_2, GL{1}), ...
                'left environment fixed point equation unfulfilled.');
            tc.assertEqual(lambda, lambda2, 'RelTol', 1e-10);
            
            % right environment
            fp_left  = insert_onespace(fixedpoint(mps, 'l_RR'), ...
                2, ~isdual(leftvspace(mpo, 1)));
            fp_right = repartition(insert_onespace(fixedpoint(mps, 'r_RR'), ...
                2, ~isdual(rightvspace(mpo, 1))), rank(GR{1}));
            
            TR = transfermatrix(mpo, mps, mps, 'Type', 'RR').';
            GR_2 = apply(TR, GR{1});
            lambda3 = overlap(GR_2(1), fp_left);
            GR_2(1) = GR_2(1) - lambda3 * fp_right;
            tc.assertTrue(isapprox(GR_2, GR{1}), ...
                'right environment fixed point equation unfulfilled.');
            tc.assertEqual(lambda, lambda3, 'RelTol', 1e-10);
        end
        
        function testQuasiEnvironments(tc)
            D = 16;
            mpo = quantum1dIsing('Symmetry', 'Z2');
            mps = initialize_mps(mpo, GradedSpace.new(Z2(0, 1), D ./ [2 2], false));
            [~, ~, lambda] = environments(mpo, mps);
            
            mpo = mpo - lambda;
            [GL, GR, lambda] = environments(mpo, mps);
            tc.assertEqual(lambda, 0, 'AbsTol', 1e-10);
            
            for charge = Z2([1 0])
                for p = [0 pi 0.5]
                    qp = InfQP.randnc(mps, mps, p, charge);
                    
%                     [GBL, GBR] = quasienvironments(mpo, qp, GL, GR);
                    
                    % left environments
                    GBL = leftquasienvironment(mpo, qp, GL, GR);
                    tc.assertEqual(norm(GBL{1}(1)), 0, 'AbsTol', 1e-10, ...
                        'left gauge quasiparticle violation.');
                    
                    T_R = transfermatrix(mpo, qp, qp, 'Type', 'RL');
                    T_B = transfermatrix(mpo, qp, qp, 'Type', 'BL');
                    
                    GBL2 = apply(T_R, GBL{1}) + apply(T_B, GL{1});
                    if istrivial(qp)
                        fp_left = insert_onespace(insert_onespace(...
                            fixedpoint(qp, 'l_RL'), ...
                            2, ~isdual(leftvspace(mpo, 1))), ...
                            4, isdual(auxspace(qp, 1)));
                        fp_left = repartition(fp_left, rank(GBL2(end)));
                        fp_right = insert_onespace(insert_onespace(...
                            fixedpoint(qp, 'r_RL'), ...
                            2, ~isdual(rightvspace(mpo, 1))), ...
                            1, ~isdual(auxspace(qp, 1)));
                        GBL2(end) = GBL2(end) - overlap(GBL2(end), fp_right) * fp_left;
                    end
                    tc.assertTrue(isapprox(GBL2, GBL{1} * exp(+1i*p)), ...
                        sprintf('left environment fixed point equation unfulfilled for p=%e, c=%d', p, charge));
                    
                    % right environments
                    GBR = rightquasienvironment(mpo, qp, GL, GR);
                    
                    T_L = transfermatrix(mpo, qp, qp, 'Type', 'LR').';
                    T_B = transfermatrix(mpo, qp, qp, 'Type', 'BR').';
                    
                    GBR2 = apply(T_L, GBR{1}) + apply(T_B, GR{1});
                    if istrivial(qp)
                        fp_left = insert_onespace(insert_onespace(...
                            fixedpoint(qp, 'l_LR'), ...
                            2, ~isdual(leftvspace(mpo, 1))), ...
                            1, ~isdual(auxspace(qp, 1)));
                        fp_right = insert_onespace(insert_onespace(...
                            fixedpoint(qp, 'r_LR'), ...
                            2, ~isdual(rightvspace(mpo, 1))), ...
                            4, isdual(auxspace(qp, 1)));
                        fp_right = repartition(fp_right, rank(GBR2(1)));
                        GBR2(1) = GBR2(1) - overlap(GBR2(1), fp_left) * fp_right;
                    end
                    tc.assertTrue(isapprox(GBR2, GBR{1} * exp(-1i*p)), ...
                        sprintf('right environment fixed point equation unfulfilled for p=%e, c=%d', p, charge));
                end
            end
            
        end
        
        function testDerivatives(tc, mpo, mps)
            [GL, GR] = environments(mpo, mps, mps);
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_AC)
                AC_ = mps.AC(i);
                [AC_.var, lambda] = eigsolve(H_AC{i}, mps.AC(i).var, 1, 'largestabs');
                tc.assertTrue(isapprox(apply(H_AC{i}, AC_), lambda * AC_.var));
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
            mpo = quantum1dIsing('J', 1, 'h', 1, 'L', Inf);
            mps = initialize_mps(mpo, CartesianSpace.new(D));
            [mps2, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyTrue(isapprox(lambda, -1.27, 'RelTol', 1e-2))
            
            mpo = quantum1dIsing('J', 1, 'h', 1, 'L', Inf, 'Symmetry', 'Z2');
            mps = initialize_mps(mpo, GradedSpace.new(Z2(0, 1), [D D] ./ 2, false));
            [mps2, lambda2] = fixedpoint(alg, mpo, mps);
            tc.verifyTrue(isapprox(lambda, -1.27, 'RelTol', 1e-2))
            
            mpo = [mpo mpo];
            mps = [mps mps];
            [mps2, lambda2] = fixedpoint(alg, mpo, mps);
            tc.verifyTrue(isapprox(lambda2/2, -1.27, 'RelTol', 5e-2))
        end
        
        function test1dHeisenberg(tc)
            alg = Vumps('which', 'smallestreal', 'maxiter', 5);
            
            mpo = quantum1dHeisenberg('Spin', 1, 'Symmetry', 'SU2');
            mpo = [mpo mpo];
            
            vspace1 = GradedSpace.new(SU2(1:2:5), [5 5 1], false);
            vspace2 = GradedSpace.new(SU2(1:2:5), [5 5 1], false);
            mps = initialize_mps(mpo, vspace1, vspace2);
            
            [gs_mps] = fixedpoint(alg, mpo, mps);
            lambda = expectation_value(gs_mps, mpo);
            tc.verifyEqual(lambda / period(mps), -1.40, 'RelTol', 1e-2);
        end
    end
end

