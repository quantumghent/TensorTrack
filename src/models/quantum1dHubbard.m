function mpo = quantum1dHubbard(u, mu, kwargs)
% Hamiltonian for the 1D Hubbard model.
%
% `math`:-\sum (c^+_i c_j + c^+_j c_i) + u\sum (1 - 2n_{up}) * (1-2n_{down}) - mu\sum (n_{up} + n_{down}):
%
% Arguments
% ---------

arguments
    u
    mu = 0
    kwargs.filling = 1
    kwargs.symmetry = 'fZ2xSU2xU1'
end

switch kwargs.symmetry
    case 'fZ2xSU2xU1'
        [p, q] = rat(kwargs.filling);
        if q > 30
            warning('filling %f is not a nice rational (%d // %d)', kwargs.filling, p, q)
        end
        
        pcharges = ProductCharge(fZ2(0, 1, 0), SU2(1, 2, 1), U1([0, 1, 2] * q - p));
        pspace = GradedSpace.new(pcharges, ones(size(pcharges)), false);
        
        acharge = ProductCharge(fZ2(1), SU2(2), U1(q));
        aspace = GradedSpace.new(acharge, 1, false);
        
        creation_L = fill_tensor(Tensor(pspace, [pspace aspace]), {sqrt(2) 1});
        annihilation_R = fill_tensor(Tensor([aspace pspace], pspace), {sqrt(2) 1});

        hopping = contract(creation_L, [-1 1 -4], annihilation_R, [1 -2 -3], 'Rank', [2 2]);
        hopping = hopping + hopping';
        
        interaction = fill_tensor(Tensor(pspace, pspace), {1 1 -1});
        
        chemical_potential = fill_tensor(Tensor(pspace, pspace), {0 2 1});
        
        mpo = repmat(InfJMpo.twosite(-hopping, ...
            u * interaction - mu * chemical_potential), 1, 2*q);
        
    otherwise
        error('tba:model', 'symmetry %s not implemented', kwargs.symmetry);
end




end
