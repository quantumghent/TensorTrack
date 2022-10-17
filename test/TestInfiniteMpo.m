classdef TestInfiniteMpo < matlab.unittest.TestCase
    % Unit tests for infinite matrix product operators.
    
    properties (TestParameter)
        mpo = struct(...
            'trivial', InfiniteMpo.Ising(), ...
            'Z2', InfiniteMpo.Ising('Symmetry', 'Z2') ...
            )
        mps = struct(...
            'trivial', UniformMps.randnc(CartesianSpace.new(2), CartesianSpace.new(4)), ...
            'Z2', UniformMps.randnc(GradedSpace.new(Z2(0,1), [1 1], false), ...
                GradedSpace.new(Z2(0,1), [4 4], false)) ...
            )
    end
    
    methods (Test, ParameterCombination='sequential')
        function testEnvironments(tc, mpo, mps)
            [GL, lambdaL] = leftenvironment(mpo, mps, mps);
            tc.assertTrue(isapprox(transferleft(mpo, mps, mps, GL{1}), GL{1} * lambdaL));
            [GR, lambdaR] = rightenvironment(mpo, mps, mps);
            tc.assertTrue(isapprox(transferright(mpo, mps, mps, GR{1}), GR{1} * lambdaR));
            
            N = period(mps);
            for w = 1:N
                tc.assertTrue(isapprox(...
                    applyleft(mpo.O{w}, mps.AL(w), conj(mps.AL(w)), GL{w}), ...
                    lambdaL^(1/N) * GL{next(w, N)}));
                tc.assertTrue(isapprox(...
                    applyright(mpo.O{w}, mps.AR(w), conj(mps.AR(w)), GR{prev(w, N)}), ...
                    lambdaR^(1/N) * GR{w}));
            end
        end
        
        function testDerivatives(tc, mpo, mps)
            
        end
    end
end

