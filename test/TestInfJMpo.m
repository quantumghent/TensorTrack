classdef TestInfJMpo < matlab.unittest.TestCase
    % Unit tests for infinite matrix product operators.
    
    properties (TestParameter)
        mpo = struct(...
            'trivial', quantum1dIsing(), ...
            'fermion', quantum1d_randomfermion() ...
            )
        mps = struct(...
            'trivial', UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(4)), ...
            'fermion', UniformMps.randnc(fZ2Space([0 1], [1 1], false), fZ2Space([0 1], [2 2], false)) ...
            )
    end
    
    methods (Test, ParameterCombination='sequential')
        function testEnvironments(tc, mpo, mps)
%             D = 16;
%             mpo = quantum1dIsing('Symmetry', 'Z2');
%             mps = initialize_mps(mpo, GradedSpace.new(Z2(0, 1), D ./ [2 2], false));
            
            [GL, GR, lambda] = environments(mpo, mps);
            
            % left environment
            fp_left  = repartition(insert_onespace(fixedpoint(mps, 'l_LL'), ...
                2, ~isdual(leftvspace(mpo, 1))), rank(GL{1}));
            fp_right = insert_onespace(fixedpoint(mps, 'r_LL'), ...
                2, ~isdual(rightvspace(mpo, 1)));
            
            TL = transfermatrix(mpo, mps, mps, 'Type', 'LL');
            GL_2 = apply(TL, GL{1});
            lambda2 = overlap(slice(GL_2, length(GL_2)), fp_right);
            GL_2.var(end) = GL_2.var(end) - lambda2 * fp_left;
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
            lambda3 = overlap(slice(GR_2, 1), fp_left);
            GR_2.var(1) = GR_2.var(1) - lambda3 * fp_right;
            tc.assertTrue(isapprox(GR_2, GR{1}), ...
                'right environment fixed point equation unfulfilled.');
            tc.assertEqual(lambda, lambda3, 'RelTol', 1e-10);
        end
        
        function testQuasiEnvironments(tc)
            %% Ising
            D = 16;
            mpo = quantum1dIsing('Symmetry', 'Z2');
            mps = initialize_mps(mpo, GradedSpace.new(Z2(0, 1), D ./ [2 2], false));
            [GL, GR] = environments(mpo, mps);
            
            for charge = Z2([1 0])
                for p = [0 pi 0.5]
                    qp = InfQP.randnc(mps, mps, p, charge);
                    
                    % left environments
                    GBL = leftquasienvironment(mpo, qp, GL, GR);
                    tc.assertEqual(norm(GBL{1}(1)), 0, 'AbsTol', 1e-10, ...
                        'left gauge quasiparticle violation.');
                    
                    T_R = transfermatrix(mpo, qp, qp, 'Type', 'RL');
                    T_B = transfermatrix(mpo, qp, qp, 'Type', 'BL');
                    
                    GBL2 = apply(T_R, GBL{1}) + apply(T_B, GL{1});
                    
                    if istrivial(qp)
                        r = rank(GBL2(end));
                        fp_left = repartition(fixedpoint(mpo, qp, 'l_RL_0'), r);
                        fp_right = fixedpoint(mpo, qp, 'r_RL_1');
                        GBL2(end) = GBL2(end) - overlap(GBL2(end), fp_right) * fp_left;
                    end
                    
                    tc.verifyTrue(isapprox(GBL2, GBL{1} * exp(+1i*p)), ...
                        sprintf('left environment fixed point equation unfulfilled for p=%e, c=%d', p, charge));

                    % right environments
                    GBR = rightquasienvironment(mpo, qp, GL, GR);
                    
                    T_L = transfermatrix(mpo, qp, qp, 'Type', 'LR').';
                    T_B = transfermatrix(mpo, qp, qp, 'Type', 'BR').';
                    
                    GBR2 = apply(T_L, GBR{1}) + apply(T_B, GR{1});
                    if istrivial(qp)
                        r = rank(GBR2(1));
                        fp_left  = fixedpoint(mpo, qp, 'l_LR_1');
                        fp_right = repartition(fixedpoint(mpo, qp, 'r_LR_0'), r);
                        GBR2(1) = GBR2(1) - overlap(GBR2(1), fp_left) * fp_right;
                    end
                    tc.assertTrue(isapprox(GBR2, GBR{1} * exp(-1i*p)), ...
                        sprintf('right environment fixed point equation unfulfilled for p=%e, c=%d', p, charge));
                end
            end
        end
        
        function testQuasiEnvironments2(tc)
            D = 16;
            %% Heisenberg
            mpo = quantum1dHeisenberg('Spin', 0.5, 'Symmetry', 'SU2');
            mpo = [mpo mpo];
            vspace = GradedSpace.new(...
                SU2(1:2:5), [2 2 1] * (D / 4), false, ...
                SU2(2:2:6), [2 2 1] * (D / 4), false);
            mps = initialize_mps(mpo, vspace(1), vspace(2));
            alg = Vumps('maxiter', 3, 'which', 'smallestreal');
            mps = fixedpoint(alg, mpo, mps);
            [GL, GR] = environments(mpo, mps);
            
            for charge = SU2(1, 3)
                for p = [0 pi 0.5]
                    qp = InfQP.randnc(mps, mps, p, charge);
                    
                    % left environments
                    GBL = leftquasienvironment(mpo, qp, GL, GR);
                    T_R = transfermatrix(mpo, qp, qp, 'Type', 'RL');
                    T_B = transfermatrix(mpo, qp, qp, 'Type', 'BL');
                    
                    for w = 1:period(mpo)
                        tc.assertEqual(norm(GBL{w}.var(1)), 0, 'AbsTol', 1e-10, ...
                            'left gauge quasiparticle violation.');
                        GBL2 = apply(T_R(w), GBL{w}) + apply(T_B(w), GL{w});
                        
                        if istrivial(qp) && w == prev(1, period(mpo))
                            r = rank(GBL2);
                            fp_left = repartition(fixedpoint(mpo, qp, 'l_RL_0', next(w, period(mpo))), r);
                            fp_right = fixedpoint(mpo, qp, 'r_RL_1', w);
                            GBL2.var(end) = GBL2.var(end) - overlap(GBL2.var(end), fp_right) * fp_left;
                        end
                        
                        tc.verifyTrue(isapprox(GBL2, GBL{next(w, period(mpo))} * exp(+1i*p/period(mpo))), ...
                            sprintf('left environment fixed point equation unfulfilled for p=%e, c=%d', p, charge));
                    end
                    
                    % right environments
                    GBR = rightquasienvironment(mpo, qp, GL, GR);
                    T_L = transfermatrix(mpo, qp, qp, 'Type', 'LR').';
                    T_B = transfermatrix(mpo, qp, qp, 'Type', 'BR').';
                    
                    for w = 1:period(mpo)
                        GBR2 = apply(T_L(w), GBR{w}) + apply(T_B(w), GR{w});
                        
                        if istrivial(qp) && w == next(1, period(mpo))
                            r = rank(GBR2);
                            fp_left = fixedpoint(mpo, qp, 'l_LR_1', prev(w, period(mpo)));
                            fp_right = repartition(fixedpoint(mpo, qp, 'r_LR_0', w), r);
                            GBR2.var(1) = GBR2.var(1) - overlap(GBR2.var(1), fp_left) * fp_right;
                        end
                        
                        tc.verifyTrue(isapprox(GBR2, GBR{prev(w, period(mpo))} * exp(-1i*p/period(mpo))), ...
                            sprintf('right environment fixed point equation unfulfilled for p=%e, c=%d', p, charge));
                    end
                end
            end
            
        end
        
        function testDerivatives(tc, mpo, mps)
            [GL, GR] = environments(mpo, mps, mps);
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_AC)
                AC_ = mps.AC{i};
                [AC_.var, lambda] = eigsolve(H_AC{i}, mps.AC{i}.var, 1, 'smallestreal');
                tc.assertTrue(isapprox(apply(H_AC{i}, AC_), lambda * AC_.var));
            end
            
            H_C = C_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_C)
                [C_, lambda] = eigsolve(H_C{i}, mps.C{i}, 1, 'smallestreal');
                tc.assertTrue(isapprox(apply(H_C{i}, C_), lambda * C_));
            end
        end
        
        function test1dIsing(tc)
            alg = Vumps('which', 'smallestreal', 'maxiter', 10);
            D = 16;
            mpo = quantum1dIsing('J', 1, 'h', 1, 'L', Inf);
            mps = initialize_mps(mpo, CartesianSpace.new(D));
            [gs, lambda] = fixedpoint(alg, mpo, mps);
            tc.verifyTrue(isapprox(lambda, -0.53, 'RelTol', 1e-2))
            
            p = 0;
            qp = InfQP.randnc(gs, gs, p);
            [qp, mu] = excitations(QPAnsatz(), mpo, qp);
            tc.verifyEqual(mu, 0.5, 'AbsTol', 1e-8);
            
            gs2 = [gs gs];
            mpo2 = [mpo mpo];
            gs2 = fixedpoint(alg, mpo2, gs2);
            
            qp2 = InfQP.randnc(gs2, gs2);
            [qp2, mu2] = excitations(QPAnsatz(), mpo2, qp2);
            
            tc.verifyEqual(mu, mu2, 'AbsTol', 1e-8);
            
            mpo = quantum1dIsing('J', 1, 'h', 1, 'L', Inf, 'Symmetry', 'Z2');
            mps = initialize_mps(mpo, GradedSpace.new(Z2(0, 1), [D D] ./ 2, false));
            [mps2, lambda2] = fixedpoint(alg, mpo, mps);
            tc.verifyTrue(isapprox(lambda2, -0.53, 'RelTol', 1e-2))
            
            p = 0;
            qp = InfQP.randnc(mps2, mps2, p, Z2(0));
            [qp, mu] = excitations(QPAnsatz(), mpo, qp);
            tc.verifyEqual(mu, 0.5, 'AbsTol', 1e-8);
            
            mpo = [mpo mpo];
            mps2 = [mps2 mps2];
            [mps2, lambda2] = fixedpoint(alg, mpo, mps2);
            tc.verifyTrue(isapprox(lambda2/2, -0.53, 'RelTol', 1e-2))
            
            p = 0;
            qp = InfQP.randnc(mps2, mps2, p, Z2(0));
            [qp, mu] = excitations(QPAnsatz(), mpo, qp);
            tc.verifyEqual(mu, 0.5, 'AbsTol', 1e-8);
        end
        
        function test1dHeisenberg(tc)
            alg = Vumps('which', 'smallestreal', 'maxiter', 100);
            
            mpo = quantum1dHeisenberg('Spin', 1, 'Symmetry', 'SU2');
            vspace1 = GradedSpace.new(SU2(2:2:6), [5 5 1], false);
            mps = initialize_mps(mpo, vspace1);
            
            [gs_mps] = fixedpoint(alg, mpo, mps);
            lambda = expectation_value(gs_mps, mpo);
            tc.verifyEqual(lambda, -1.401, 'RelTol', 1e-2);
            
            p = pi;
            charge = SU2(3);
            qp = InfQP.randnc(gs_mps, gs_mps, p, charge);
            
            [qp, mu] = excitations(QPAnsatz(), mpo, qp);
            mu
            
        end
    end
end

function H = quantum1d_randomfermion()

pspace = fZ2Space([0 1], [1 1], false);
D = Tensor.randnc([pspace pspace], [pspace pspace]);
D = normalize(D + D');

H = InfJMpo.twosite(D);

end

