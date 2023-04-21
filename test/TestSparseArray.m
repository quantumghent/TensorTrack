classdef TestSparseArray < matlab.unittest.TestCase
    % Unit tests for sparse arrays.
    
    properties
        tol = 1e-12
    end
    
    methods (Test)
        
        function test_subsref(tc)
            a = SparseArray([4, 4, 4; 2, 2, 1; 2, 3, 2], [3; 5; 1], [4, 4, 4]);
            tc.verifyTrue(isequal(a(1, 2, 1), 0));
            tc.verifyTrue(isequal(a(4, 4, 4), 3));
            tc.verifyTrue(isequal(size(a(2, :, :)), [1, 4, 4]));
            tc.verifyTrue(isequal(size(a([0, 4, 0])), [4, 1, 4]));
            tc.verifyTrue(isequal(size(a([4, 4, 0])), [1, 1, 4]));
        end

        function test_basics(tc)
            a = SparseArray.zeros([2, 2, 2]);
            sz1 = size(a(:, :, 1));
            b = SparseArray.zeros([2, 2, 2]);
            b(1, 1, 1) = 1;
            sz2 = size(b(:, :, 1));
            tc.verifyTrue(isequal(sz1, sz2));
        end

    end
end
