classdef TestCharge < matlab.unittest.TestCase
    % TestCharge - Unit tests for charges.
    
    properties
        tol = 1e-14
    end
    
    properties (TestParameter)
        smallset = struct( ...
            'Z1', Z1, ...
            'Z2', Z2([false, true]), ...
            'Z4', ZN(4, 0:3), ...
            'U1', U1(-2:2), ...
            'O2', O2([0 0 1 2], [0 1 2 2]), ...
            'SU2', SU2(1:4), ...
            'A4', A4(1:4), ...
            'Z2xU1', ProductCharge(Z2(0, 0, 0, 1, 1, 1), U1(-1, 0, 1, -1, 0, 1))...
            )
    end
    
    methods (Test)
        function fusionrules(testCase, smallset)
            for a = smallset, for b = smallset
                    cs = a * b;
                    N = Nsymbol(repmat(a, 1, length(cs)), repmat(b, 1, length(cs)), cs);
                    verifyEqual(testCase, sum(qdim(cs) .* N), qdim(a) * qdim(b), ...
                        'Fusion rules must preserve quantum dimensions.');
                    verifyTrue(testCase, all(N > 0), ...
                        'Nsymbol must be nonnegative for fusion product.');
                end, end
        end
        
        function conj(tc, smallset)
            for a = smallset
                abar = conj(a);
                
                tc.verifyTrue(conj(abar) == a, 'Conj should be an involution.');
                tc.verifyTrue(Nsymbol(a, abar, one(a)) == 1, ...
                    'a * abar should contain one.');
            end
        end
        
        function Fmove(testCase, smallset)
            % Fmove - Compare Fsymbol with fusiontensor.
            %   Fmove(testCase, smallset)
            %       verifies if an F-move on the fusion tensors is compatible with
            %       the Fsymbol.
            for a = smallset, for b = smallset, for c = smallset
                        for e = a * b
                            F3 = fusiontensor(a, b, e);
                            for f = b * c
                                F2 = conj(fusiontensor(b, c, f));
                                for d = intersect(e * c, a * f)
                                    F1 = conj(fusiontensor(a, f, d));
                                    F4 = fusiontensor(e, c, d);
                                    testCase.assertEqual(contract(F1, [1 6 4 -4], ...
                                        F2, [2 3 6 -3], F3, [1 2 5 -1], F4, [5 3 4 -2]), ...
                                        Fsymbol(a, b, c, d, e, f) * qdim(d), ...
                                        'AbsTol', testCase.tol, ...
                                        'Fsymbol incompatible with fusiontensor.');
                                end
                            end
                        end
                    end, end, end
        end
        
        function qdims(testCase, smallset)
            for a = smallset
                F = Fsymbol(a, conj(a), a, a, one(a), one(a));
                testCase.verifyTrue(isapprox(qdim(a), abs(1 / F(1))), ...
                    'Fsymbol and qdim incompatible');
            end
        end
        
        function Funitary(testCase, smallset)
            % Funitary - Test if the Fsymbols are unitary.
            %   Funitary(testCase, smallset)
            %       verifies whether the Fsymbol constitutes a unitary operation
            %       when interpreted as a matrix from (emn) to (fkl).
            for a = smallset, for b = smallset, for c = smallset
                        for d = prod([a, b, c])
                            Fmat = Fmatrix(a, b, c, d);
                            testCase.assertEqual(Fmat' * Fmat, eye(size(Fmat)), ...
                                'AbsTol', testCase.tol);
                            testCase.assertEqual(Fmat * Fmat', eye(size(Fmat)), ...
                                'AbsTol', testCase.tol);
                        end
                    end, end, end
        end
        
        function Fone(testCase, smallset)
            % Fone - Test if the Fsymbol is correct for the trivial charge.
            e = one(smallset);
            testCase.assertEqual(Fsymbol(e, e, e, e, e, e), 1, ...
                'AbsTol', testCase.tol);
        end
        
        function pentagon(testCase, smallset)
            % pentagon - Test if the Fsymbol fulfills the pentagon equation.
            for a = smallset, for b = smallset, for c = smallset, for d = smallset
                            for f = a * b, for h = c * d
                                    for g = f * c, for i = b * h
                                            for e = intersect(g * d, a * i)
                                                if hasmultiplicity(fusionstyle(smallset))
                                                    lhs = contract(...
                                                        Fsymbol(f, c, d, e, g, h), [-1 -2 -3 1], ...
                                                        Fsymbol(a, b, h, e, f, i), [-4 1 -5 -6]);
                                                    rhs = zeros(size(lhs));
                                                    for j = b * c
                                                        rhs = rhs + contract( ...
                                                            Fsymbol(a, b, c, g, f, j), [-4 -1 1 2], ...
                                                            Fsymbol(a, j, d, e, g, i), [2 -2 3 -6], ...
                                                            Fsymbol(b, c, d, i, j, h), [1 3 -3 -5]);
                                                    end
                                                else
                                                    lhs = Fsymbol(f, c, d, e, g, h) * ...
                                                        Fsymbol(a, b, h, e, f, i);
                                                    rhs = sum(arrayfun(@(j) ...
                                                        Fsymbol(a, b, c, g, f, j) * ...
                                                        Fsymbol(a, j, d, e, g, i) * ...
                                                        Fsymbol(b, c, d, i, j, h), b * c));
                                                end
                                                testCase.verifyTrue(isapprox(lhs, rhs, ...
                                                    'AbsTol', testCase.tol, 'RelTol', testCase.tol));
                                            end
                                        end, end
                                end, end
                        end, end, end, end
        end
        
        function Rmove(testCase, smallset)
            % Rmove - Compare Rsymbol with fusiontensor.
            %   Rmove(testCase, smallset)
            %       verifies if an R-move on the fusion tensors is compatible with
            %       the Fsymbol.
            for a = smallset, for b = smallset
                    for c = a * b
                        testCase.verifyEqual(fusiontensor(a, b, c), ...
                            contract(fusiontensor(b, a, c), [-2 -1 -3 1], ...
                            Rsymbol(a, b, c), [-4 1]), 'AbsTol', testCase.tol);
                    end
                end, end
            % TODO non-symmetric braiding tests
        end
        
        function hexagon(testCase, smallset)
            % hexagon - Test if the F- and Rsymbol fulfill the hexagon equation.
            %   hexagon(testCase, smallset)
            for a = smallset, for b = smallset, for c = smallset
                        for e = c * a, for g = c * b
                                for d = intersect(e * b, a * g)
                                    if hasmultiplicity(fusionstyle(smallset))
                                        lhs = contract(Rsymbol(c, a, e), [-1 1], ...
                                            Fsymbol(a, c, b, d, e, g), [1 -2 2 -4], ...
                                            Rsymbol(b, c, g), [2 -3]);
                                        rhs = zeros(size(lhs));
                                        for f = a * b
                                            rhs = rhs + contract(...
                                                Fsymbol(c, a, b, d, e, f), [-1 -2 1 2], ...
                                                Rsymbol(c, f, d), [2 3], ...
                                                Fsymbol(a, b, c, d, f, g), [1 3 -3 -4]);
                                        end
                                    else
                                        lhs = Rsymbol(c, a, e) * Fsymbol(a, c, b, d, e, g) ...
                                            * Rsymbol(b, c, g);
                                        rhs = sum(arrayfun(@(f) Fsymbol(c, a, b, d, e, f) * ...
                                            Rsymbol(c, f, d) * Fsymbol(a, b, c, d, f, g), ...
                                            a * b));
                                    end
                                    testCase.verifyEqual(lhs, rhs, 'AbsTol', testCase.tol);
                                end
                            end, end
                    end, end, end
        end
        
        function cumprod(testCase, smallset)
            % cumprod - Test if the cumprod functionality works.
            %   cumprod(testCase, smallset)
            
            % scalar case:
            for a = smallset, for b = smallset
                    [d, v] = cumprod([a; b]);
                    for c = smallset
                        [d, v] = cumprod([a; b; c]);
                    end
                end, end
            
            % matrix case:
            [d, v] = cumprod(combvec(smallset, smallset));
            [d, v] = cumprod(combvec(smallset, smallset, smallset));
            [d, v] = cumprod(combvec(smallset, smallset, smallset, smallset));
        end
    end
    
end
