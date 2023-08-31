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
            tc.assumeTrue(depth == 1, 'Test ill-defined for multiline PEPS.')
            
            mpo = tc.random_mpo(spaces, depth, width, dualdepth, dualwidth, samebot);
            
            % test transpose (co)domain
            tc.assertTrue(isequal(domain(mpo.'), codomain(mpo)'));
            tc.assertTrue(isequal(codomain(mpo.'), domain(mpo)'));
            
            % test ctranspose(co)domain
            tc.assertTrue(isequal(domain(mpo'), codomain(mpo)));
            tc.assertTrue(isequal(codomain(mpo'), domain(mpo)));
            
            
            % test tensor behavior
            tc.assumeTrue(width < 2, 'FiniteMpo -> Tensor contraction not reasonable.')
            
            % test domain and codomain
            mpo_tensor = Tensor(mpo);
            tc.assertTrue(isequal(mpo.domain, mpo_tensor.domain), ...
                'domain should remain fixed after conversion.');
            tc.assertTrue(isequal(mpo.codomain, mpo_tensor.codomain), ...
                'codomain should remain fixed after conversion.');
            
            % test transpose
            tc.assertTrue(isapprox(Tensor(mpo).', Tensor(mpo.')), ...
                'transpose should not change mpo');
            
            % test ctranspose
            tc.assertTrue(isapprox(Tensor(mpo)', Tensor(mpo')), ...
                'ctranspose should not change mpo');
            
            % test apply
            v = initialize_fixedpoint(mpo);
            twistinds = find(isdual(space(v, 1:nspaces(v)))); % compensate for supertrace rules when using contract versus mtimes
            tc.assertTrue(isapprox(mpo.apply(v), mpo_tensor * twist(v, twistinds)));
        end
        
        function testFixedpoints(tc, spaces, depth, width, dualdepth, dualwidth, samebot)
            mpo = tc.random_mpo(spaces, depth, width, dualdepth, dualwidth, samebot);

            [V, D] = eigsolve(mpo, 'MaxIter', 1e3, 'KrylovDim', 32);
            tc.assertTrue(isapprox(mpo.apply(V), D * V));
            
            v2 = insert_onespace(V);
            v3 = apply(mpo, v2);
        end
        
        function testDerivatives(tc, spaces, depth, width, dualdepth, dualwidth, samebot)
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
            pspace = spaces(1); horzspace = spaces(2); vertspace = spaces(3);
            top = TestPeps.random_peps_unitcell(pspace, horzspace, vertspace, 1, 1, dualdepth, dualwidth);
            bot = TestPeps.random_peps_unitcell(pspace', vertspace, horzspace, 1, 1, dualdepth, dualwidth);
            O = PepsSandwich(top{1}, bot{1});
            t = MpoTensor(O);
        end
        
end
    
    methods (Static)
        function A = random_peps_unitcell(pspace, horzspace, vertspace, depth, width, dualdepth, dualwidth)
            % Spit out random peps unit cell with all kinds of arrows for testing
            
            if dualwidth, horzspace = horzspace'; end
            if dualdepth, vertspace = vertspace'; end
            pspaces = repmat(pspace, depth, width); pspaces(2:2:depth*width) = conj(pspaces(2:2:depth*width));
            westspaces = repmat(horzspace, depth, width);
            southspaces = repmat(vertspace, depth, width);
            
            % staggering
            if ~mod(width, 2)
                westspaces(2:2:depth*width) = conj(westspaces(2:2:depth*width));
                eastspaces = westspaces;
            else
                eastspaces = conj(westspaces);
            end
            if ~mod(depth, 2)
                southspaces(2:2:depth*width) = conj(southspaces(2:2:depth*width));
                northspaces = southspaces;
            else
                northspaces = conj(southspaces);
            end
            
            % fill up cell
            for d = depth:-1:1, for w = width:-1:1
                A{d, w} = PepsTensor.randnc(pspaces(d, w), ...
                    westspaces(d, w), southspaces(d, w), eastspaces(d, w), northspaces(d, w));
            end, end
        end
        
        function T = random_sandwich(pspace, horzspace, vertspace, depth, width, dualdepth, dualwidth, samebot)
            top = TestPeps.random_peps_unitcell(pspace, horzspace, vertspace, depth, width, dualdepth, dualwidth);
            if samebot
                bot = cellfun(@conj, top, 'UniformOutput', false);
            else
                bot = TestPeps.random_peps_unitcell(pspace', horzspace, vertspace, depth, width, dualdepth, dualwidth);
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
    end
end

