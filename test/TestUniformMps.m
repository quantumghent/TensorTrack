classdef TestUniformMps < matlab.unittest.TestCase
    % Unit tests for uniform matrix product states.
    
    properties (TestParameter)
        mps = struct(...
            'trivial', UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(4)), ...
            'trivial2', UniformMps.randnc(CartesianSpace.new([2 2]), CartesianSpace.new([12 12])), ...
            'fermion1', UniformMps.randnc(GradedSpace.new(fZ2(0,1), [1 1], false), ...
                GradedSpace.new(fZ2(0,1), [2 2], false)), ...
            'fermion2', UniformMps.randnc(GradedSpace.new(fZ2(0,1), [1 1], false), ...
                GradedSpace.new(fZ2(0,1), [2 2], true)), ...
            'fermion3', UniformMps.randnc(GradedSpace.new(fZ2(0,1), [1 1], true), ...
                GradedSpace.new(fZ2(0,1), [2 2], false)), ...
            'fermion4', UniformMps.randnc(GradedSpace.new(fZ2(0,1), [1 1], true), ...
                GradedSpace.new(fZ2(0,1), [2 2], true)), ...
            'haldane', UniformMps.randnc(GradedSpace.new(SU2(2), 1, false, SU2(2), 1, false), ...
                GradedSpace.new(SU2(1:2:5), [5 3 2], false, SU2(2:2:6), [5 2 1], false)) ...
            )
    end
    
    methods (Test)
        function testCanonical(tc, mps)
            mps = canonicalize(mps);
            AL = mps.AL; AR = mps.AR; C = mps.C; AC = mps.AC;
            tc.assertTrue(all(isisometry(mps.AL, 'left')), ...
                'AL should be a left isometry.');
            tc.assertTrue(all(isisometry(mps.AR, 'right')), ...
                'AR should be a right isometry.');
            
            for w = 1:period(mps)
                ALC = repartition(multiplyright(AL(w), C(w)), rank(AC(w)));
                CAR = multiplyleft(AR(w), C(prev(w, period(mps))));
                tc.assertTrue(isapprox(ALC, AC(w)) && isapprox(AC(w), CAR), ...
                    'AL, AR, C and AC should be properly related.');
            end
        end
        
        function testDiagonalC(tc, mps)
            mps2 = diagonalizeC(mps);
            f = fidelity(mps, mps2, 'Verbosity', Verbosity.diagnostics);
            tc.assertTrue(isapprox(f, 1), 'Diagonalizing C should not alter the state.');
        end
        
        function testFixedpoints(tc, mps)
            for top = ["L" "R"]
                if strcmp(top, "L")
                    T = mps.AL;
                else
                    T = mps.AR;
                end
                
                for bot = ["L" "R"]
                    if strcmp(bot, "L")
                        B = mps.AL;
                    else
                        B = mps.AR;
                    end
                    
                    rhoL = fixedpoint(mps, sprintf('l_%c%c', top, bot));
                    tc.assertTrue(isapprox(rhoL, applyleft(T, conj(B), rhoL)), ...
                        'rho_left should be a fixed point.');
                    rhoR = fixedpoint(mps, sprintf('r_%c%c', top, bot));
                    tc.assertTrue(isapprox(rhoR, applyright(T, conj(B), rhoR)));
                end
            end 
        end
    end
end

