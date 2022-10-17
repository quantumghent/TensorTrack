classdef TestFiniteMpo < matlab.unittest.TestCase
    % Unit tests for finite matrix product operators.
    
    properties (TestParameter)
        mpo = struct('trivial0', ...
            FiniteMpo(Tensor.randnc(ComplexSpace.new(2, false, 4, true, 2, false), []), ...
            {}, Tensor.randnc(ComplexSpace.new(2, false, 4, false, 2, true), []))..., ... [2 4 2], [false true true]), {}, Tensor.randnc([2 4 2], [true false false])), ...
            ...'trivial1', ...
            ...    FiniteMpo(Tensor.randnc([2 4 2]), {Tensor.randnc([4 2 4 2])}, ...
            ...    Tensor.randnc([2 4 2])) ...
                )
    end
    
    methods (Test)
        function testFixedpoints(tc, mpo)
            
        end
        
        function testProperties(tc, mpo)
            mpo_tensor = Tensor(mpo);
            tc.assertTrue(isequal(mpo.domain, mpo_tensor.domain));
            tc.assertTrue(isequal(mpo.codomain, mpo_tensor.codomain));
            
            v = initialize_fixedpoint(mpo);
            tc.assertTrue(isapprox(mpo.apply(v), mpo_tensor * v));
        end
        
        function testTransposes(tc, mpo)
            
        end
    end
end

