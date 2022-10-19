classdef TestFiniteMpo < matlab.unittest.TestCase
    % Unit tests for finite matrix product operators.
    
    properties (TestParameter)
        mpo = struct('trivial0', ...
            FiniteMpo.randnc(ComplexSpace.new(2, false, 3, false), ComplexSpace.new(4, false)), ...
            'trivial1', ...
                FiniteMpo.randnc(ComplexSpace.new(2, false, 4, true, 2, true), ...
                ComplexSpace.new(2, true, 3, false)) ...
            )
    end
    
    methods (Test)
        function testFixedpoints(tc, mpo)
            [V, D] = eigsolve(mpo);
            tc.assertTrue(isapprox(mpo.apply(V), D * V));
        end
        
        function testProperties(tc, mpo)
            mpo_tensor = Tensor(mpo);
            tc.assertTrue(isequal(mpo.domain, mpo_tensor.domain), ...
                'domain should remain fixed after conversion.');
            tc.assertTrue(isequal(mpo.codomain, mpo_tensor.codomain), ...
                'codomain should remain fixed after conversion.');
            
            v = initialize_fixedpoint(mpo);
            tc.assertTrue(isapprox(mpo.apply(v), mpo_tensor * v));
        end
        
        function testTransposes(tc, mpo)
            tc.assertTrue(isapprox(Tensor(mpo)', Tensor(mpo')), ...
                'ctranspose should not change mpo');
            tc.assertTrue(isequal(domain(mpo'), codomain(mpo)));
            tc.assertTrue(isequal(codomain(mpo'), domain(mpo)));
            tc.assertTrue(isapprox(Tensor(mpo).', Tensor(mpo.')), ...
                'transpose should not change mpo');
            tc.assertTrue(isequal(domain(mpo.'), codomain(mpo)'));
            tc.assertTrue(isequal(codomain(mpo.'), domain(mpo)'));
            
        end
    end
end

