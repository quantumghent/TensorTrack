function mpo = quantum1dHeisenberg(kwargs)
arguments
    kwargs.Spin = 1
    kwargs.J = 1
    kwargs.h = 0
    kwargs.L = Inf     % size of system
    kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'U1', 'SU2'})} = 'Z1'
end

J = kwargs.J;
h = kwargs.h;

Q = SU2(2 * kwargs.Spin + 1);

switch kwargs.Symmetry
    case 'SU2'
        assert(isscalar(J) || all(J == J(1)), ...
            'Different spin couplings not invariant under SU2');
        assert(h == 0, 'Magnetic field not invariant under SU2');
        
        pSpace = GradedSpace.new(Q, 1, false);
        vSpace = GradedSpace.new(SU2(3), 1, false);
        tSpace = one(vSpace);
        
        s = kwargs.Spin;
        L = Tensor.ones([tSpace pSpace], [pSpace vSpace]);
        L = L * (-J(1) * (s^2 + s));
        R = Tensor.ones([vSpace pSpace], [pSpace tSpace]);
        
        cod = SumSpace([one(vSpace) vSpace one(vSpace)], pSpace);
        dom = SumSpace(pSpace, [one(vSpace), vSpace, one(vSpace)]);
        O = MpoTensor.zeros(cod, dom);
        O(1, 1, 1, 1) = 1;
        O(3, 1, 3, 1) = 1;
        O(1, 1, 2, 1) = L;
        O(2, 1, 3, 1) = R;
        
    otherwise
        error('TBA');
end

mpo = InfJMpo(O);

if isfinite(kwargs.L)
    mpo = open_boundary_conditions(mpo, L);
end

end

