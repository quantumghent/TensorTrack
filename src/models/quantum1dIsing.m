function mpo = quantum1dIsing(kwargs)
% Hamiltonian for the 1D transverse-field Ising model.
%
% .. math::
%   H = -J \left(\sum_{\langle ij \rangle} S_i^x S_j^x  + h \sum_{i} S_i^z \right).
%
% Keyword arguments
% -----------------
% 'J' : :class:`double`
%   :math:`ZZ` coupling, defaults to :code:`1`.
%
% 'h' : :class:`double`
%   relative transverse field strength, defaults to :code:`1`.
%
% 'L' : :class:`int`
%   system size, defaults to :code:`Inf`.
%
% 'Symmetry' : :class:`char`
%   symmetry group ('Z1', 'Z2' or 'fZ2'), defaults to :code:`'Z1'`.
%
% Returns
% -------
% mpo : :class:`.InfJMpo`
%   Ising Hamiltonian as a Jordan block MPO.

arguments
    kwargs.J = 1
    kwargs.h = 1
    kwargs.L = Inf     % size of system
    kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'Z2', 'fZ2'})} = 'Z1'
end

J = kwargs.J;
h = kwargs.h;

sigma_x = [0 1; 1 0];
sigma_z = [1 0; 0 -1];

switch kwargs.Symmetry
    case 'Z1'
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
        
    case 'Z2'
        pSpace = GradedSpace.new(Z2(0, 1), [1 1], false);
        vSpace = GradedSpace.new(Z2(1), 1, false);
        trivSpace = one(pSpace);
        
        Sx_l = fill_matrix(Tensor([trivSpace pSpace], [pSpace vSpace]), {1 1});
        Sx_r = fill_matrix(Tensor([vSpace pSpace], [pSpace trivSpace]), {1 1});
        Sz = fill_matrix(Tensor([trivSpace pSpace], [pSpace trivSpace]), {1 -1});
        
        cod = SumSpace([one(vSpace) vSpace one(vSpace)], pSpace);
        dom = SumSpace(pSpace, [one(vSpace), vSpace, one(vSpace)]);
        O = MpoTensor.zeros(cod, dom);
        O(1, 1, 1, 1) = 1;
        O(3, 1, 3, 1) = 1;
        O(1, 1, 2, 1) = -J * Sx_l;
        O(2, 1, 3, 1) = Sx_r;
        O(1, 1, 3, 1) = (-J * h) * Sz;
        
    case 'fZ2'
        c_left = c_min('side', 'left');
        c_right = c_min('side', 'right');
        cdag_left = c_plus('side', 'left');
        cdag_right = c_plus('side', 'right');
        
        % twosite terms
        cdagc = contract(cdag_left, [-1 1 -4], c_right, [1 -2 -3], 'Rank', [2 2]);
        ccdag = contract(c_left, [-1 1 -4], cdag_right, [1 -2 -3], 'Rank', [2 2]);
        cc = contract(c_left, [-1 1 -4], c_right, [1 -2 -3], 'Rank', [2 2]);
        cdagcdag = contract(cdag_left, [-1 1 -4], cdag_right, [1 -2 -3], 'Rank', [2 2]);
        H_twosite = -J * (cdagc + ccdag + cdagcdag + cc);
        
        % onesite terms
        n = c_number();
        H_onesite = -J * 2 * h * (n - 1 / 2 * n.eye(n.codomain, n.domain));
        
        mpo = InfJMpo.twosite(H_twosite, H_onesite);
        return
        
    otherwise
        error('models:argerror', 'invalid symmetry');
end

mpo = InfJMpo(O);

if isfinite(kwargs.L)
    mpo = open_boundary_conditions(InfJMpo(O), kwargs.L);
end

end

