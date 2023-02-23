classdef TestPeps < matlab.unittest.TestCase
    % Unit tests for PEPS transfer matrices.
    
    properties (TestParameter)
        spaces = struct(...
            'cartesian', CartesianSpace.new([2, 3, 4, 5]), ...
            'complex', ComplexSpace.new(2, false, 3, false, 4, false, 5, false), ...
            'Z2', GradedSpace.new(Z2(0, 1), [1 1], false, Z2(0, 1), [1 2], false, ...
                Z2(0, 1), [2 1], false, Z2(0, 1), [3 3], false), ...
            'fZ2', GradedSpace.new(fZ2(0, 1), [1 1], false, fZ2(0, 1), [1 2], false, ...
                fZ2(0, 1), [2 1], false, fZ2(0, 1), [5 5], false), ...
            'U1', GradedSpace.new(U1(0, 1, -1), [1 1 1], false, U1(0, 1, -1), [1 2 2], false, ...
                U1(0, 1, -1), [2 1 1], false, U1(0, 1, -1), [2 1 1], false), ...
            'SU2', GradedSpace.new(SU2(1, 2), [1 1], false, SU2(2), [2], false, ...
                SU2(2), [1], false, SU2(1, 2), [2 1], false) ...
            )
        depth = {1}
        width = {1, 2}
        dualdepth = {false, true}
        dualwidth = {false, true}
        samebot = {true, false}
    end
    
    methods (Test, ParameterCombination='exhaustive')
        
        function testFiniteMpo(tc, spaces, depth, width, dualdepth, dualwidth, samebot)
            tc.assumeTrue(braidingstyle(spaces) ~= BraidingStyle.Fermionic, 'fZ2 test broken')
            tc.assumeTrue(depth == 1, 'Test ill-defined for multiline PEPS.')
            
            mpo = tc.random_mpo(spaces, depth, width, dualdepth, dualwidth, samebot);
            
            % test transpose (co)domain
            tc.assertTrue(isequal(domain(mpo.'), codomain(mpo)'));
            tc.assertTrue(isequal(codomain(mpo.'), domain(mpo)'));
            
            % test ctranspose(co)domain
            tc.assertTrue(isequal(domain(mpo'), codomain(mpo)));
            tc.assertTrue(isequal(codomain(mpo'), domain(mpo)));
            
            
            % test tensor behavior
            tc.assumeTrue(width < 3 || (fusionstyle(spaces) <= FusionStyle.Unique || width < 2), ...
                'FiniteMpo -> Tensor contraction not reasonable.')
            
            % test domain and codomain
            mpo_tensor = tc.mpo_to_tensor(mpo);
            tc.assertTrue(isequal(mpo.domain, mpo_tensor.domain), ...
                'domain should remain fixed after conversion.');
            tc.assertTrue(isequal(mpo.codomain, mpo_tensor.codomain), ...
                'codomain should remain fixed after conversion.');
            
            % test apply
            v = initialize_fixedpoint(mpo);
            tc.assertTrue(isapprox(mpo.apply(v), mpo_tensor * v));
            
            % test transpose
            tc.assertTrue(isapprox(tc.mpo_to_tensor(mpo).', tc.mpo_to_tensor(mpo.')), ...
                'transpose should not change mpo');
            
            % test ctranspose
            tc.assertTrue(isapprox(tc.mpo_to_tensor(mpo)', tc.mpo_to_tensor(mpo')), ...
                'ctranspose should not change mpo');
        end
        
        function testFixedpoints(tc, spaces, depth, width, dualdepth, dualwidth, samebot)
            tc.assumeTrue(braidingstyle(spaces) ~= BraidingStyle.Fermionic, 'fZ2 test broken')

            mpo = tc.random_mpo(spaces, depth, width, dualdepth, dualwidth, samebot);

            [V, D] = eigsolve(mpo, 'MaxIter', 1e3, 'KrylovDim', 32);
            tc.assertTrue(isapprox(mpo.apply(V), D * V));
            
            v2 = insert_onespace(V);
            v3 = apply(mpo, v2);
        end
        
        function testDerivatives(tc, spaces, depth, width, dualdepth, dualwidth, samebot)
            tc.assumeTrue(braidingstyle(spaces) ~= BraidingStyle.Fermionic, 'fZ2 test broken')

            mpo = tc.random_inf_mpo(spaces, depth, width, dualdepth, dualwidth, samebot);
            vspaces = repmat({spaces(end)}, 1, width);
            mps = mpo.initialize_mps(vspaces{:});
            
            [GL, GR] = environments(mpo, mps, mps);
            
            H_AC = AC_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_AC)
                AC_ = mps.AC(i);
                [AC_.var, lambda] = eigsolve(H_AC{i}, mps.AC(i).var, 1, 'largestabs');
                AC_2 = apply(H_AC{i}, AC_);
                tc.assertTrue(isapprox(AC_2, lambda * AC_.var, 'RelTol', 1e-6));
            end
            
            H_C = C_hamiltonian(mpo, mps, GL, GR);
            for i = 1:numel(H_C)
                [C_, lambda] = eigsolve(H_C{i}, mps.C(i), 1, 'largestabs');
                tc.assertTrue(isapprox(apply(H_C{i}, C_), lambda * C_, 'RelTol', 1e-6));
            end
        end
        
        function testMpoTensor(tc, spaces, dualdepth, dualwidth)
            tc.assumeTrue(braidingstyle(spaces) ~= BraidingStyle.Fermionic, 'fZ2 test broken')

            pspace = spaces(1); horzspace = spaces(2); vertspace = spaces(3);
            top = TestPeps.random_peps(pspace, horzspace, vertspace, 1, 1, dualdepth, dualwidth);
            bot = TestPeps.random_peps(pspace', vertspace, horzspace, 1, 1, dualdepth, dualwidth);
            O = PepsSandwich(top{1}, bot{1});
            t = MpoTensor(O);
        end
        
end
    
    methods (Static)
        function A = random_peps(pspace, horzspace, vertspace, depth, width, dualdepth, dualwidth)
            % spit out random peps unit cell
            A = cell(depth, width);
            if dualwidth, horzspace = horzspace'; end
            if dualdepth, vertspace = vertspace'; end

            if mod(depth, 2)
                vertspaces = [vertspace vertspace'];
            else
                vertspaces = [vertspace vertspace];
            end

            if mod(width, 2)
                horzspaces = [horzspace horzspace'];
            else
                horzspaces = [horzspace horzspace];
            end

            codomainspaces = reshape([vertspaces; horzspaces], 1, []);

            for d = 1:depth
                for w = 1:width
                    dmsp = pspace;
                    cdmsp = codomainspaces;
                    if ~mod(depth, 2) && ~mod(d, 2), cdmsp = conj(cdmsp); dmsp = conj(dmsp); end
                    if ~mod(width, 2) && ~mod(w, 2), cdmsp = conj(cdmsp); dmsp = conj(dmsp); end
                    A{d, w} = Tensor.randnc(dmsp, cdmsp);
                    A{d, w} = A{d, w} / norm(A{d, w});
                end
            end
        end
        
        function T = random_sandwich(pspace, horzspace, vertspace, depth, width, dualdepth, dualwidth, samebot)
            top = TestPeps.random_peps(pspace, horzspace, vertspace, depth, width, dualdepth, dualwidth);
            if samebot
                bot = cellfun(@conj, top, 'UniformOutput', false);
            else
                bot = TestPeps.random_peps(pspace', horzspace, vertspace, depth, width, dualdepth, dualwidth);
            end
            T = cellfun(@(t, b) PepsSandwich(t, b), top, bot, 'UniformOutput', false);
        end
        
        function mpo = random_mpo(spaces, depth, width, dualdepth, dualwidth, samebot)
            pspace = spaces(1); horzspace = spaces(2); vertspace = spaces(3); vspace = spaces(4);
            O = TestPeps.random_sandwich(pspace, horzspace, vertspace, depth, width, dualdepth, dualwidth, samebot);
            for d = depth:-1:1
                L = MpsTensor(Tensor.randnc([vspace leftvspace(O{d, 1})'], vspace));
                R = MpsTensor(Tensor.randnc([vspace rightvspace(O{d, end})'], vspace));
                mpo(d, 1) = FiniteMpo(L, O(d, :), R);
            end
        end
        
        function mpo = random_inf_mpo(spaces, depth, width, dualdepth, dualwidth, samebot)
            pspace = spaces(1); horzspace = spaces(2); vertspace = spaces(3);
            O = TestPeps.random_sandwich(pspace, horzspace, vertspace, depth, width, dualdepth, dualwidth, samebot);
            mpo = InfMpo(O);
        end
        
        function t = mpo_to_tensor(mpo)
            % contract FiniteMpo{PepsSandwich} to tensor...
            assert(depth(mpo) == 1, 'not implemented for 1 < depth');
            W = length(mpo);
            inds_top = arrayfun(@(x) [  3 + 3*(x-1), ...
                                        1 + 3*(x-1), ...
                                        -(2 + 2*(x-1)), ...
                                        1 + 3*x, ...
                                        -(4*W + 5 - 2*x)], ...
                                        1:W,  'UniformOutput', false);
            inds_bot = arrayfun(@(x) [  3 + 3*(x-1), ...
                                        2 + 3*(x-1), ...
                                        -(3 + 2*(x-1)), ...
                                        2 + 3*x, ...
                                        -(4*W + 4 - 2*x)], ...
                                        1:W,  'UniformOutput', false);
            tops = cellfun(@(x) x.top.var, mpo.O, 'UniformOutput', false);
            bots = cellfun(@(x) x.bot.var, mpo.O, 'UniformOutput', false);
            args = [tops; inds_top; bots; inds_bot];
            t = contract(mpo.L, [-1, 2, 1, -(4*W + 4)], args{:}, mpo.R, [-(2*W + 3), 1 + 3*W, 2 + 3*W, -(2 + 2*W)], ...
                    'Rank', [2*W+2, 2*W+2]);

        end
        
    end
end

