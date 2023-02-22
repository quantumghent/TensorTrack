classdef TestFiniteMps < matlab.unittest.TestCase
    % Unit tests for uniform matrix product states.
    
    properties (TestParameter)
        A = struct(...
            'trivial', {{Tensor.randnc(CartesianSpace.new([4 2]), CartesianSpace.new(4))}}, ...
            'trivial2', {{Tensor.randnc(CartesianSpace.new([4 2]), CartesianSpace.new(5)), ...
                        Tensor.randnc(CartesianSpace.new([5 2]), CartesianSpace.new(4))}}, ...
            'trivial3', {{Tensor.randnc(CartesianSpace.new([4 2]), CartesianSpace.new(5)), ...
                        Tensor.randnc(CartesianSpace.new([5 2]), CartesianSpace.new(6)), ...
                        Tensor.randnc(CartesianSpace.new([6 2]), CartesianSpace.new(4))}}, ...
            'fermion1', {{Tensor.randnc(...
                GradedSpace.new(fZ2(0,1), [2 2], false, fZ2(0,1), [1 1], false), ...
                GradedSpace.new(fZ2(0,1), [2 2], false))}}, ...
            'fermion2', {{Tensor.randnc(...
                GradedSpace.new(fZ2(0,1), [2 2], true, fZ2(0,1), [1 1], false), ...
                GradedSpace.new(fZ2(0,1), [2 2], true))}}, ...
            'fermion3', {{Tensor.randnc(...
                GradedSpace.new(fZ2(0,1), [2 2], false, fZ2(0,1), [1 1], true), ...
                GradedSpace.new(fZ2(0,1), [2 2], false))}}, ...
            'fermion4', {{Tensor.randnc(...
                GradedSpace.new(fZ2(0,1), [2 2], true, fZ2(0,1), [1 1], true), ...
                GradedSpace.new(fZ2(0,1), [2 2], true))}}, ...
            'haldane', {{Tensor.randnc(GradedSpace.new(SU2(1:2:5), [5 3 2], false, SU2(2), 1, false), ...
                                    GradedSpace.new(SU2(2:2:6), [5 2 1], false)), ...
                        Tensor.randnc(GradedSpace.new(SU2(2:2:6), [5 2 1], false, SU2(2), 1, false), ...
                                    GradedSpace.new(SU2(1:2:5), [5 3 2], false))}} ...
            )
    end
    
    properties
        tol = 1e-10
    end
    
    methods (Test)
        function testFullMps(tc)
            L = 12;
            P = CartesianSpace.new(2);
            mps = FiniteMps.new([], P, 'L', L);
            
            % test conversion
            psi = Tensor(mps);
            psi2 = Tensor(FiniteMps(psi));
            tc.verifyTrue(isapprox(psi, psi2, 'RelTol', tc.tol));
            
            % test norms
            tc.verifyEqual(norm(psi), norm(mps), 'RelTol', tc.tol);
            tc.verifyEqual(sqrt(abs(overlap(mps, mps))), norm(psi), 'RelTol', tc.tol);
            mps_normed = normalize(mps);
            tc.verifyEqual(norm(mps_normed), 1, 'RelTol', tc.tol);
            tc.verifyEqual(norm(Tensor(mps_normed)), 1, 'RelTol', tc.tol);
            
            % test moving orthogonality center
            for i = 1:length(mps)
                tc.verifyEqual(overlap(movegaugecenter(mps_normed, i), mps_normed), 1, ...
                    'RelTol', tc.tol);
            end
        end
        
        function testDiagonalC(tc, A)
            mps = UniformMps(A);
            mps2 = diagonalizeC(mps);
            f = fidelity(mps, mps2);
            tc.assertTrue(isapprox(f, 1), 'Diagonalizing C should not alter the state.');
        end
        
        function testFixedpoints(tc, A)
            mps = UniformMps(A);
            for top = ["L" "R"]
                for bot = ["L" "R"]
                    T = transfermatrix(mps, mps, 'Type', sprintf('%c%c', top, bot));
                    rhoL = fixedpoint(mps, sprintf('l_%c%c', top, bot));
                    rhoR = fixedpoint(mps, sprintf('r_%c%c', top, bot));
                    tc.verifyTrue(isapprox(rhoL, T.apply(rhoL)), ...
                        sprintf('rho_left should be a %c%c fixed point.', top, bot));
                    tc.verifyTrue(isapprox(rhoR, apply(T', rhoR)), ...
                        sprintf('rho_right should be a %c%c fixed point.', top, bot));
                end
            end
        end
        
        function testTransferEigs(tc, A)
           mps = UniformMps(A);
           
           [V, D] = transfereigs(mps, mps, 1, 'largestabs');
           [~, charges] = matrixblocks(V);
           [V2, D2] = transfereigs(mps, mps, 1, 'largestabs', 'Charge', one(charges));
           
        end
    end
end

