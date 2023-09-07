function S = sigma_plus(spin, symmetry)
arguments
    spin = 1/2
    symmetry = 'Z1'
end

switch symmetry
%     case 'Z1'
        
    case 'U1'
        charges = U1((-2 * spin):2:(2 * spin));
        degeneracies = ones(size(charges));
        pspace = GradedSpace.new(charges, degeneracies, false);
        aspace = GradedSpace.new(U1(-2), 1, false);
        
        S = Tensor.zeros([pspace aspace], pspace);
        [mblocks, bcharges] = matrixblocks(S);
        for i = 1:length(charges)
            if charges(i) == U1(2 * spin), continue; end
            mblocks{charges(i) == bcharges} = 2 * pauliterm(spin, i, i + 1);
        end
        S = S.fill_matrix(mblocks, bcharges);
        
    otherwise
        error('models:TBA', 'not implemented');
end

end

