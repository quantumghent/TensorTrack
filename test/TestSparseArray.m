classdef TestSparseArray < matlab.unittest.TestCase
    % Unit tests for sparse arrays.
    
    properties
        tol = 1e-12
    end
    
    methods (Test)
        
        function test_basics(tc)
            a = SparseArray.zeros([2, 2, 2]);
            sz1 = size(a(:, :, 1));
            b = SparseArray.zeros([2, 2, 2]);
            b(1, 1, 1) = 1;
            sz2 = size(b(:, :, 1));
            tc.verifyTrue(isequal(sz1, sz2));
        end
        
        function test_doctests(tc)
            % test examples in docstrings
            
            % constructor
            subs = [1 1 1; 1 1 3; 2 2 2; 4 4 4; 1 1 1; 1 1 1];
            vals = [0.5; 1.5; 2.5; 3.5; 4.5; 5.5];
            sz = [4, 4, 4];
            a = SparseArray(subs, vals, sz);
            tc.verifyTrue(isequal(size(a), [4, 4, 4]));
            tc.verifyTrue(isequal(a(1, 1, 3), 1.5));
            tc.verifyTrue(isequal(a(4, 4, 4), 3.5));
            tc.verifyTrue(isequal(a(1, 1, 1), 10.5));
            
            % groupind
            a = SparseArray.random([2, 3, 4, 5, 6], .1);
            b = groupind(a, [3, 2]);
            tc.verifyTrue(isequal(size(b), [24, 30]));
            
            % minus
            a = SparseArray.random([4 3 2], .1);
            b = SparseArray.random([4 3 2], .1);
            tc.verifyTrue(isa(a - b, 'SparseArray'));
            tc.verifyTrue(isa(a - 5, 'double'));
            tc.verifyTrue(isa(a - 0, 'double'));
            tc.verifyTrue(isa(a - full(a), 'double'));
            
            % mrdivide
            b = a / 3;
            tc.verifyTrue(isa(b, 'SparseArray'));
            tc.verifyTrue(isequal(nonzeros(b), nonzeros(a) / 3));
            
            % ndims
            tc.verifyTrue(isequal(ndims(a), 3));
            
            % plus
            tc.verifyTrue(isa(a + b, 'SparseArray'));
            tc.verifyTrue(isa(a + 5, 'double'));
            tc.verifyTrue(isa(a + 0, 'double'));
            tc.verifyTrue(isa(a + full(a), 'double'));
            
            % rdivide
            tc.verifyTrue(isa(a ./ 5, 'SparseArray'));
            tc.verifyTrue(isa(5 ./ a, 'double'));
            tc.verifyTrue(isa(a ./ full(a), 'SparseArray'));
            tc.verifyTrue(isa(full(a) ./ a, 'double'));
            
            % squeeze
            c = squeeze(SparseArray.random([2, 1, 3], 0.5));
            tc.verifyTrue(isa(c, 'SparseArray'));
            tc.verifyTrue(isequal(size(c), [2, 3]));
            d = squeeze(SparseArray([1, 1, 1], 1, [1, 1, 1]));
            tc.verifyTrue(isa(d, 'double'));
            tc.verifyTrue(isequal(size(d), [1, 1]));
            
            % subsref
            a = SparseArray([4, 4, 4; 2, 2, 1; 2, 3, 2], [3; 5; 1], [4, 4, 4]);
            tc.verifyTrue(isequal(a(1, 2, 1), 0));
            tc.verifyTrue(isequal(a(4, 4, 4), 3));
            tc.verifyTrue(isequal(size(a(2, :, :)), [1, 4, 4]));
            tc.verifyTrue(isequal(size(a([0, 4, 0])), [4, 1, 4]));
            tc.verifyTrue(isequal(size(a([4, 4, 0])), [1, 1, 4]));
            
            % times
            a = SparseArray.random([4 3 2], .1);
            b = SparseArray.random([4 3 2], .1);
            tc.verifyTrue(isa(a .* b, 'SparseArray'));
            tc.verifyTrue(isa(a .* 5, 'SparseArray'));
            tc.verifyTrue(isa(a .* 0, 'SparseArray'));
            tc.verifyTrue(isa(a .* full(a), 'SparseArray'));
        end

    end
end
