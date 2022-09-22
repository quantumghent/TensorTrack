classdef TestFusionTree < matlab.unittest.TestCase
    % TestFusionTree - Unit tests for fusion trees.
    
    %#ok<*PROPLC>
    
    properties (ClassSetupParameter)
        weight = {'small', 'medium'}%, 'large'}
        charge = {'A4', 'Z1', 'Z2', 'fZ2', 'U1', 'O2', 'SU2', 'Z2xU1'}
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
                case 'fZ2'
                    chargeset = fZ2([0 1]);
                case 'U1'
                    switch weight
                        case 'small'
                            chargeset = U1(-1:1);
                        case {'medium', 'large'}
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
                    
                case 'A4'
                    switch weight
                        case 'small'
                            chargeset = A4([1 4]);
                        case {'medium', 'large'}
                            chargeset = A4(1:4);
                    end
                    
            end
            
            %% Setup weight
            switch weight
                case 'small'
                    legs = 0:3;
                    maxTrees = 15;
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
%                     args(2, :) = num2cell(randi([0 1], 1, rank1+rank2));
                    args(2, 1:2:end) = {false};
                    args(2, 2:2:end) = {true};
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
                    assertTrue(tc, all(isallowed(tc.trees{i,j})), ...
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
                    for p = ps(randperm(size(ps, 1), min(size(ps, 1), 5)), :).'
                        lvl = randperm(f.legs);
                        indout = randi([0 f.legs]);
                        
                        [c1, f1] = braid(f, p.', lvl, [indout, f.legs - indout]);
                        
                        % basic properties
                        assertEqual(tc, size(c1), [length(f) length(f1)], ...
                            'Coefficients have the wrong size.');
                        assertEqual(tc, full(abs(c1).^2 * qdim(f1.coupled)), ...
                            qdim(f.coupled), 'AbsTol', tc.tol, ...
                            'Braiding must preserve centernorm.');
                        assertEqual(tc, isallowed(f1), true(length(f1), 1), ...
                            'Output trees are not allowed.');
                        
                        % compatible with fusiontensor
                        if istwistless(braidingstyle(f))
                            a1_cell = fusiontensor(f);
                            a2_cell = fusiontensor(f1);
                            for j = 1:length(f)
                                a1 = permute(a1_cell{j}, p.');
                                a2 = zeros(size(a1));
                                [~, col, val] = find(c1(j, :));
                                for k = 1:length(val)
                                    a2 = a2 + val(k) * a2_cell{col(k)};
                                end
                                assertEqual(tc, a2, a1, 'AbsTol', tc.tol, ...
                                    'Braiding should be compatible with fusiontensors.');
                            end
                        end
                        
                        % invertible
                        [c2, f2] = braid(f1, invperm(p.'), lvl(p), f.rank);
                        assertEqual(tc, c1 * c2, speye(length(f)), ...
                            'AbsTol', tc.tol, 'RelTol', tc.tol, ...
                            'Braiding should be invertible.');
                        assertEqual(tc, f2, f, 'Braiding should be invertible.');
                    end
                end
            end
        end
        
        function traces(tc)
            for i = 1:numel(tc.trees)
                if isempty(tc.trees{i}) || tc.trees{i}.legs < 2 || tc.trees{i}.rank(2) ~= 0
                    continue;
                end
                
                f = tc.trees{i};
                for j = 1:f.legs-1
                    [c1, f1] = elementary_trace(f, j);
                    
                    assertEqual(tc, size(c1), [length(f) length(f1)], ...
                        'Coefficients have the wrong size.');
                    assertEqual(tc, isallowed(f1), true(length(f1), 1), ...
                        'Output trees are not allowed.');
                    assertEqual(tc, unique(f1), f1, ...
                        'Output trees are not unique.');
                    
                    % compatible with fusiontensor
                    if issymmetric(braidingstyle(f))
                        a1_cell = fusiontensor(f);
                        a2_cell = fusiontensor(f1);
                        
                        for k = 1:length(f)
                            if f.uncoupled(k, j) ~= conj(f.uncoupled(k, j+1))
                                assertEqual(tc, norm(c1(k, :)), 0);
                            else
                                a1 = tensortrace(a1_cell{k}, [1:j-1 j j j+2:f.legs+1]);
                                a2 = zeros(size(a1));
                                [~, col, val] = find(c1(k, :));
                                for l = 1:length(val)
                                    a2 = a2 + val(l) * a2_cell{col(l)};
                                end
                                
                                assertEqual(tc, a2, a1, 'AbsTol', tc.tol, ...
                                    'Tracing should be compatible with fusiontensors.');
                            end
                        end
                    end
                end
            end
        end
        
        function repartitioning(tc)
            for i = 1:numel(tc.trees)
                if isempty(tc.trees{i}) || tc.trees{i}.legs < 2
                    continue;
                end
                if rand < tc.testWeight
                    f = tc.trees{i};
                    for n = 0:f.legs
                        [c1, f1] = repartition(f, [n f.legs-n]);
                        
                        assertEqual(tc, size(c1), [length(f) length(f1)], ...
                            'Coefficients have the wrong size.');
                        assertEqual(tc, full(abs(c1).^2 * qdim(f1.coupled)), ...
                            qdim(f.coupled), 'AbsTol', tc.tol, 'RelTol', tc.tol, ...
                            'Repartition must preserve centernorm.');
                        assertEqual(tc, isallowed(f1), true(length(f1), 1), ...
                            'Output trees are not allowed.');
                        
                        [c2, f2] = repartition(f1, f.rank);
                        assertEqual(tc, c1 * c2, speye(length(f)), 'AbsTol', tc.tol, ...
                            'RelTol', tc.tol, 'Repartition should be invertible.');
                        assertEqual(tc, f, f2, 'Repartition should be invertible.');
                        
                        if issymmetric(braidingstyle(f))
                            a1_cell = fusiontensor(f);
                            a2_cell = fusiontensor(f1);
                            for j = 1:length(f)
                                a1 = a1_cell{j};
                                a2 = zeros(size(a1));
                                [~, col, val] = find(c1(j, :));
                                for k = 1:length(val)
                                    a2 = a2 + val(k) * a2_cell{col(k)};
                                end
                                assertEqual(tc, a2, a1, 'AbsTol', tc.tol, ...
                                    'Repartition should be compatible with fusiontensors.');
                            end
                        end
                    end
                end
            end
        end
    end
end
