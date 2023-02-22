classdef TestFiniteMps < matlab.unittest.TestCase
    % Unit tests for finite matrix product states.
    
    properties
        tol = 1e-10
    end
    
    methods (Test)
        function testFullMps(tc)
            L = 12;
            P = ComplexSpace.new(2, false);
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
    end
end

