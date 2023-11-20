function S = sigma_z(spin, symmetry)
arguments
    spin = 1/2
    symmetry = 'Z1'
end

switch symmetry
    case 'U1'
        charges = U1((-2 * spin):2:(2 * spin));
        degeneracies = ones(size(charges));
        pspace = GradedSpace.new(charges, degeneracies, false);
        
        S = Tensor.zeros(pspace, pspace);
        [mblocks, bcharges] = matrixblocks(S);
        for i = 1:length(charges)
            mblocks{charges(i) == bcharges} = spin + 1 - i;
        end
        S = S.fill_matrix(mblocks, bcharges);

    otherwise
        error('models:TBA', 'not implemented');
end

end

