function mpo = quantum1dHeisenberg(kwargs)
% Hamiltonian for the 1D Heisenberg model.
%
% .. math::
%   H = J \sum_{\langle ij \rangle} \vec{S}_i \cdot \vec{S}_j + h \sum_{i} S_i^z
%
% Keyword arguments
% -----------------
% 'Spin' : :class:`double`
%   halfinteger or integer spin label, defaults to :code:`1`.
%
% 'J' : :class:`double`
%   exchange coupling, defaults to :code:`1`.
%
% 'h' : :class:`double`
%   magnetic field, defaults to :code:`0`.
%
% 'L' : :class:`int`
%   system size, defaults to :code:`Inf`.
%
% 'Symmetry' : :class:`char`
%   symmetry group ('U1' or 'SU2'), defaults to :code:`'SU2'`.
%
% Returns
% -------
% mpo : :class:`.InfJMpo`
%   Heisenberg Hamiltonian as a Jordan block MPO.

arguments
    kwargs.Spin = 1
    kwargs.J = 1
    kwargs.h = 0
    kwargs.L = Inf     % size of system
    kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'U1', 'SU2'})} = 'SU2'
end

J = kwargs.J;
h = kwargs.h;

Q = SU2(2 * kwargs.Spin + 1);

switch kwargs.Symmetry
    case 'SU2'
        assert(isscalar(J) || all(J == J(1)), ...
            'Different spin couplings not invariant under SU2');
        J = J(1);
        assert(h == 0, 'Magnetic field not invariant under SU2');
        
        H2 = J * sigma_exchange(kwargs.Spin, 'SU2');
        H1 = [];
        
    case 'U1'
        if isscalar(J)
            Jxy = J;
            Jz = J;
        else
            assert(length(J) == 2)
            Jxy = J(1);
            Jz = J(2);
        end
        assert(h == 0, 'TBA');
        
        Splus = sigma_plus(kwargs.Spin, 'U1');
        Smin = sigma_min(kwargs.Spin, 'U1');
        Sz = sigma_z(kwargs.Spin, 'U1');
        
        H2 = Jxy/2 * (contract(Splus, [-1 1 -3], conj(Splus), [-4 1 -2], 'Rank', [2 2]) + ...
            contract(Smin, [-1 1 -3], conj(Smin), [-4 1 -2], 'Rank', [2 2]));
        if Jz ~= 0
            H2 = H2 + Jz * contract(Sz, [-1 -3], Sz, [-2 -4], 'Rank', [2 2]);
        end
        
        if h == 0
            H1 = [];
        else
            H1 = h * Sz;
        end
        
    otherwise
        error('TBA');
end

mpo = InfJMpo.twosite(H2, H1);

if isfinite(kwargs.L)
    mpo = open_boundary_conditions(mpo, L);
end

end

