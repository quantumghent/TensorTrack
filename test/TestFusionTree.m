classdef TestFusionTree < matlab.unittest.TestCase
    % TestFusionTree - Unit tests for fusion trees.
    
    %#ok<*PROPLC>
    
    properties (ClassSetupParameter)
        weight = {'small', 'medium'}%, 'large'}
        charge = {'Z1', 'Z2', 'U1', 'O2', 'SU2', 'Z2xU1'}
    end
    
    methods (TestClassSetup)
        function generateTrees(tc, charge, weight)
            rng(123);
            
            %% Setup charges
            switch charge
                case 'Z1'
                    chargeset = Z1;
                case 'Z2'
                    chargeset = Z2([0 1]);
                case 'U1'
                    switch weight
                        case 'small'
                            chargeset = U1(-1:1);
                        case 'medium'
                            chargeset = U1(-2:2);
                        case 'large'
                            chargeset = U1(-2:2);
                    end
                case 'O2'
                    switch weight
                        case 'small'
                            chargeset = O2([0 0 1], [0 1 2]);
                        case 'medium'
                            chargeset = O2([0 0 1 2], [0 1 2 2]);
                        case 'large'
                            chargeset = O2([0 0 1 2 3], [0 1 2 2 2]);
                    end
                case 'SU2'
                    switch weight
                        case 'small'
                            chargeset = SU2(1:2);
                        case 'medium'
                            chargeset = SU2(1:4);
                        case 'large'
                            chargeset = SU2(1:6);
                    end
                case 'Z2xU1'
                    switch weight
                        case 'small'
                            chargeset = ProductCharge(Z2([0 0 0 1 1 1]), U1([-1:1 -1:1]));
                        case 'medium'
                            chargeset = ProductCharge(Z2([0 0 0 1 1 1]), U1([-1:1 -1:1]));
                        case 'large'
                            chargeset = ProductCharge(Z2([0 0 0 0 0 1 1 1 1 1]), ...
                                U1([-2:2 -2:2]));
                    end
            end
            
            %% Setup weight
            switch weight
                case 'small'
                    legs = 0:3;
                    maxTrees = 25;
                    maxLegs = 4;
                    tc.testWeight = 0.5;
                case 'medium'
                    legs = 0:4;
                    maxTrees = 100;
                    maxLegs = 5;
                    tc.testWeight = 0.25;
                case 'large'
                    legs = 0:4;
                    maxTrees = Inf;
                    maxLegs = 6;
                    tc.testWeight = 0.5;
            end
            
            %% Setup trees
            tc.trees = cell(length(legs));
            for rank1 = legs
                for rank2 = legs
                    if rank1 + rank2 == 0
                        tc.trees{rank1+1, rank2+1} = FusionTree.new([0 0]);
                        continue;
                    end
                    if rank1 + rank2 > maxLegs
                        continue;
                    end
                    args = cell(2, rank1+rank2);
                    args(1, :) = {chargeset};
                    args(2, :) = num2cell(randi([0 1], 1, rank1+rank2));
                    trees = FusionTree.new([rank1 rank2], args{:});
                    nTrees = length(trees);
                    while nTrees < 1 || nTrees > maxTrees
                        if all(cellfun(@length, args(1,:)) == 1)
                            break;
                        end
                        ind = randi([1 rank1+rank2]);
                        l = length(args{1, ind});
                        if l == 1,  continue;   end
                        args{1, ind}(randi(l)) = [];
                        newtrees = FusionTree.new([rank1 rank2], args{:});
                        nTrees = length(trees);
                        if nTrees == 0
                            break;
                        end
                        trees = newtrees;
                    end
                    assert(nTrees > 0);
                    tc.trees{rank1+1, rank2+1} = trees;
                end
            end
        end
    end
    
    properties
        tol = 1e-14
        trees
        testWeight
    end
    
    methods (Test)
        function trees_properties(tc)
            for i = 1:size(tc.trees, 1)
                for j = 1:size(tc.trees, 2)
                    if i == 1 && j == 1 || isempty(tc.trees{i,j})
                        continue;
                    end
                    verifyTrue(tc, all(isallowed(tc.trees{i,j})), ...
                        'Generated invalid fusion trees.');
                end
            end
        end
        
        function braiding(tc)
            for i = 1:numel(tc.trees)
                if isempty(tc.trees{i}) || tc.trees{i}.legs < 2
                    continue;
                end
                if rand < tc.testWeight
                    f = tc.trees{i};
                    ps = perms(1:f.legs);
                    for p = ps(randperm(size(ps, 2), min(size(ps, 2), 5)), :).'
                        lvl = randperm(f.legs);
                        indout = randi([0 f.legs]);
                        
                        [c1, f1] = braid(f, p.', lvl, [indout, f.legs - indout]);
                        
                        % basic properties
                        verifyEqual(tc, size(c1), [length(f) length(f1)], ...
                            'Coefficients have the wrong size.');
                        verifyEqual(tc, full(abs(c1).^2 * qdim(f1.coupled)), ...
                            qdim(f.coupled), 'AbsTol', tc.tol, ...
                            'Braiding must preserve centernorm.');
                        verifyEqual(tc, isallowed(f1), true(length(f1), 1), ...
                            'Output trees are not allowed.');
                        
                        % compatible with fusiontensor
                        if issymmetric(braidingstyle(f))
                            a1_cell = fusiontensor(f);
                            a2_cell = fusiontensor(f1);
                            for j = 1:length(f)
                                a1 = permute(a1_cell{j}, p.');
                                a2 = zeros(size(a1));
                                [~, col, val] = find(c1(j, :));
                                for k = 1:length(val)
                                    a2 = a2 + val(k) * a2_cell{col(k)};
                                end
                                verifyEqual(tc, a2, a1, 'AbsTol', tc.tol, ...
                                    'Braiding should be compatible with fusiontensors.');
                            end
                        end
                        
                        % invertible
                        [c2, f2] = braid(f1, invperm(p.'), lvl(p), f.rank);
                        verifyEqual(tc, c1 * c2, speye(length(f)), ...
                            'AbsTol', tc.tol, 'RelTol', tc.tol, ...
                            'Braiding should be invertible.');
                        verifyEqual(tc, f2, f, 'Braiding should be invertible.');
                    end
                end
            end
        end
        %       function repartitioning(tc)
        %          for i = 1:numel(tc.trees)
        %             if isempty(tc.trees{i}) || tc.trees{i}.legs < 2
        %                continue;
        %             end
        %             if rand < tc.testWeight
        %                f = tc.trees{i};
        %                for n = 0:f.legs
        %                   [f1, c1] = repartition(f, [n f.legs-n]);
        %                   verifyEqual(tc, full(abs(c1).^2 * qdim(f1.coupled)), ...
        %                      qdim(f.coupled), 'AbsTol', tc.tol, 'RelTol', tc.tol, ...
        %                      'Repartition must preserve centernorm.');
        %                   verifyEqual(tc, isallowed(f1), true(size(f1), 1), ...
        %                      'Output trees are not allowed.');
        %
        %                   [f2, c2] = repartition(f1, f.rank);
        %                   verifyEqual(tc, c1 * c2, speye(size(f)), 'AbsTol', tc.tol, ...
        %                      'RelTol', tc.tol, 'Repartition should be invertible.');
        %                   verifyEqual(tc, f, f2, 'Repartition should be invertible.');
        %
        %                   if issymmetric(braidingstyle(f))
        %                      a1_cell = arrayfun(@double, f, 'UniformOutput', false);
        %                      a2_cell = arrayfun(@double, f1, 'UniformOutput', false);
        %                      for j = 1:length(f)
        %                         a1 = a1_cell{j};
        %                         a2 = zeros(size(a1));
        %                         [~, col, val] = find(c1(j, :));
        %                         for k = 1:length(val)
        %                            a2 = a2 + val(k) * a2_cell{col(k)};
        %                         end
        %                         verifyEqual(tc, a2, a1, 'AbsTol', tc.tol, ...
        %                            'Repartition should be compatible with fusiontensors.');
        %                      end
        %                   end
        %                end
        %             end
        %          end
        %       end
        
        %       end
        %
        %       function permute_fusiontensor(testCase, smallset, legs)
        %          trees = generateFusionTrees(testCase, smallset, legs);
        %          ps = perms(1:legs);
        %          for p = ps(randperm(size(ps, 2), min(size(ps, 2), 5)), :).'
        %             [t1, c1] = permute(trees, p.');
        %             for j = 1:length(trees)
        %                array1 = permute(cast(slice(trees, j), 'double'), ...
        %                   [p.' legs+1]);
        %                array2 = zeros(size(array1));
        %                [~, col, val] = find(c1(j,:));
        %                for k = 1:length(val)
        %                   array2 = array2 + val(k) * cast(slice(t1, col(k)), 'double');
        %                end
        %                verifyEqual(testCase, array1, array2, 'AbsTol', testCase.tol, ...
        %                   'Permute should be compatible with fusiontensors.');
        %             end
        %          end
        %       end
        %
        %       function permute_properties(testCase, smallset, legs)
        %          trees = generateFusionTrees(testCase, smallset, legs);
        %
        %          treeargs = cell(2, legs);
        %          for i = 1:legs
        %             treeargs{1,i} = unique(trees.uncoupled(:, i));
        %             treeargs{2,i} = trees.arrows(i);
        %          end
        %          coupled = unique(trees.coupled);
        %
        %          ps = perms(1:legs);
        %          for p = ps(randperm(size(ps, 2), min(size(ps, 2), 5)), :).'
        %             [t1, c1] = permute(trees, p.');
        %
        %             verifyEqual(testCase, size(c1), [length(trees) length(t1)], ...
        %                'Coefficients have the wrong size.');
        %             verifyEqual(testCase, full(vecnorm(c1)), ones(1, length(trees)), ...
        %                'AbsTol', testCase.tol, 'Norm of coefficients should be 1.');
        %             verifyEqual(testCase, c1 * c1', speye(length(trees)), ...
        %                'AbsTol', testCase.tol, 'Permute should be unitary.');
        %             verifyEqual(testCase, c1' * c1, speye(length(t1)), ...
        %                'AbsTol', testCase.tol, 'Permute should be unitary.');
        %
        %             verifyTrue(testCase, issorted(t1), ...
        %                'Permute should return sorted trees.');
        %
        %             tempargs = treeargs(:, p.');
        %             t3 = FusionTree.new(tempargs{:});
        %             lia = ismember(t3.coupled, coupled);
        %             t3.charges = t3.charges(lia, :);
        %             if ~isempty(t3.vertices)
        %                t3.vertices = t3.vertices(lia,:);
        %             end
        %             verifyEqual(testCase, t1, t3, ...
        %                'Permute should transform full basis to full basis.');
        %
        %             [t2, c2] = permute(t1, invperm(p));
        %             verifyEqual(testCase, c2, c1', 'AbsTol', testCase.tol, ...
        %                'RelTol', testCase.tol, 'Permute should be invertible.');
        %             verifyEqual(testCase, t2, trees, ...
        %                'Permute should be invertible.');
        %          end
        %       end
        %
        %       function yangbaxter(testCase, smallset, legs)
        %          trees = generateFusionTrees(testCase, smallset, legs);
        %          for i = 1:legs-2
        %             [t1, c1] = artinbraid(trees, i);
        %             [t1, c2] = artinbraid(t1, i+1);
        %             [t1, c3] = artinbraid(t1, i);
        %
        %             [t2, c4] = artinbraid(trees, i+1);
        %             [t2, c5] = artinbraid(t2, i);
        %             [t2, c6] = artinbraid(t2, i+1);
        %
        %             verifyEqual(testCase, ...
        %                c1 * c2 * c3, c4 * c5 * c6, 'AbsTol', testCase.tol, ...
        %                'Yang-Baxter equation not satisfied.');
        %             verifyEqual(testCase, t1, t2, ...
        %                'Yang-Baxter equation not satisfied.');
        %          end
        %       end
    end
    
    %    methods
    %       function trees = generateFusionTrees(testCase, smallset, legs)
    %          % generateFusionTrees - Generate some trees for testing.
    %          %   trees = generateFusionTrees(testCase, smallset, legs)
    %
    %          rng(123);
    %          args = cell(2, legs);
    %          for i = 1:legs
    %             args{1, i} = smallset(randperm(length(smallset), ...
    %                randi([1 length(smallset)])));
    %             args{2, i} = randi([0 1]);
    %          end
    %          trees = FusionTree.new(args{:});
    %
    %          % get reasonable amount of trees:
    %          while length(trees) > testCase.maxTrees
    %             coupleds = trees.coupled;
    %             coupled = unique(trees.coupled);
    %             coupled = coupled(randi([1 length(coupled)]));
    %             lia = ismember(coupleds, coupled);
    %             if all(lia)
    %                break
    %             end
    %             trees.charges = trees.charges(~lia,:);
    %             if ~isempty(trees.vertices)
    %                trees.vertices = trees.vertices(~lia,:);
    %             end
    %             trees = sort(trees);
    %          end
    %       end
    %    end
    
end
