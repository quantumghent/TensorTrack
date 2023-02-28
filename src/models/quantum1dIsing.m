function mpo = quantum1dIsing(kwargs)
arguments
    kwargs.J = 1
    kwargs.h = 1
    kwargs.L = Inf     % size of system
    kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'Z2'})} = 'Z1'
end

J = kwargs.J;
h = kwargs.h;

sigma_x = [0 1; 1 0] / 2;
sigma_z = [1 0; 0 -1] / 2;

if strcmp(kwargs.Symmetry, 'Z1')
    pSpace = CartesianSpace.new(2);
    vSpace = one(pSpace);
    trivSpace = one(pSpace);
    
    S = Tensor([vSpace pSpace], [pSpace vSpace]);
    Sx = fill_matrix(S, sigma_x);
    Sz = fill_matrix(S, sigma_z);
    
    cod = SumSpace([vSpace vSpace vSpace], pSpace);
    dom = SumSpace(pSpace, [vSpace vSpace vSpace]);
    O = MpoTensor.zeros(cod, dom);
    O(1, 1, 1, 1) = 1;
    O(3, 1, 3, 1) = 1;
    O(1, 1, 2, 1) = -J * Sx;
    O(2, 1, 3, 1) = Sx;
    O(1, 1, 3, 1) = (-J * h) * Sz;
    
else
    pSpace = GradedSpace.new(Z2(0, 1), [1 1], false);
    vSpace = GradedSpace.new(Z2(1), 1, false);
    trivSpace = one(pSpace);
    
    Sx_l = fill_matrix(Tensor([trivSpace pSpace], [pSpace vSpace]), {1 1}) / 2;
    Sx_r = fill_matrix(Tensor([vSpace pSpace], [pSpace trivSpace]), {1 1}) / 2;
    Sz = fill_matrix(Tensor([trivSpace pSpace], [pSpace trivSpace]), {1 -1}) / 2;
    
    cod = SumSpace([one(vSpace) vSpace one(vSpace)], pSpace);
    dom = SumSpace(pSpace, [one(vSpace), vSpace, one(vSpace)]);
    O = MpoTensor.zeros(cod, dom);
    O(1, 1, 1, 1) = 1;
    O(3, 1, 3, 1) = 1;
    O(1, 1, 2, 1) = -J * Sx_l;
    O(2, 1, 3, 1) = Sx_r;
    O(1, 1, 3, 1) = (-J * h) * Sz;
end

mpo = InfJMpo(O);

if isfinite(kwargs.L)
    mpo = open_boundary_conditions(InfJMpo(O), kwargs.L);
end

end

