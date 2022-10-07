classdef TestSparseTensor < matlab.unittest.TestCase
    % Unit tests for sparse arrays of tensors.
    
    properties
        tol = 1e-12
    end
    
    methods (Test)
        function basic_linear_algebra(tc)
            spaces = CartesianSpace.new([2 3 4 5]);
            szs = cell(1, 4);
            for i = 1:length(szs)
                szs{i} = spaces(randperm(length(spaces), randi([2 3])));
            end
            t1 = generatesparsetensor([2 2], szs, 0.3);
            a = rand();
            
            tc.verifyTrue(isapprox(norm(t1)^2, dot(t1, t1), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol), 'norm and dot incompatible.');
            
            tc.verifyTrue(isapprox(norm(a .* t1), abs(a) * norm(t1), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                'norm and scalar multiplication incompatible.');
            tc.verifyTrue(isapprox(t1 + t1, 2 .* t1, ...
                'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                '2t and t+t incompatible.');
            
            tc.verifyTrue(isapprox(-t1, t1 .* (-1), 'AbsTol', tc.tol), ...
                '=t and t * (-1) incompatible.');
            tc.verifyTrue(isapprox(norm(normalize(t1)), 1), ...
                'normalize should result in unit norm.');
            
            t2 = generatesparsetensor([2 2], szs, 0.3);
            b = rand();
            
            tc.verifyTrue(isapprox(dot(b .* t2, a .* t1), conj(b) * a * dot(t2, t1), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol) && ...
                isapprox(dot(t2, t1), conj(dot(t1, t2)), ...
                'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                'Dot should be sesquilinear.');
        end
        
        function permute_via_inner(tc)
            spaces = CartesianSpace.new([2 3 4 5]);
            szs = cell(1, 4);
            for i = 1:length(szs)
                szs{i} = spaces(randperm(length(spaces), randi([2 3])));
            end
            t1 = generatesparsetensor([2 2], szs, 0.3);
            t2 = generatesparsetensor([2 2], szs, 0.3);
            
            inner = dot(t1, t2);
            for i = 0:4
                ps = perms(1:ndims(t1)).';
                for p = ps(:, randperm(size(ps, 2), min(size(ps, 2), 20)))
                    t3 = tpermute(t1, p.', [i 4-i]);
                    for j = 1:ndims(t1)
                        tc.assertTrue(isequal(space(t1, p(j)), space(t3, j)), ...
                            'incorrect spaces after permutation.');
                    end
                    tc.assertTrue(...
                        isapprox(norm(t1), norm(t3), 'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                        'Permute should preserve norms.')
                    
                    t4 = tpermute(t2, p.', [i 4-i]);
                    for j = 1:ndims(t1)
                        tc.assertTrue(isequal(space(t2, p(j)), space(t4, j)), ...
                            'incorrect spaces after permutation.');
                    end
                    tc.assertTrue(...
                        isapprox(dot(t3, t4), inner, 'AbsTol', tc.tol, 'RelTol', tc.tol), ...
                        'Permute should preserve inner products.');
                end
            end
        end
    end
end

function t = generatesparsetensor(rank, spaces, sparsity)

sz = cellfun(@length, spaces);
ind = ind2sub_(sz, 1:prod(sz));

for i = prod(sz):-1:1
    for j = length(sz):-1:1
        localsz(j) = spaces{j}(ind(i,j));
    end
    var(i) = Tensor.rand(localsz(1:rank(1)), localsz(rank(1)+1:end)', 'Rank', rank);
end

li1 = rand([1, length(var)]) < sparsity;
t = SparseTensor(ind(li1, :), var(li1), sz);

end