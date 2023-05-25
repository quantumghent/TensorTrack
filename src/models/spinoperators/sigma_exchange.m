function S = sigma_exchange(spin, symmetry)
arguments
    spin = 1/2
    symmetry = 'Z1'
end

switch symmetry
%     case 'Z1'
        
    case 'SU2'
        pspace = GradedSpace.new(SU2(2 * spin + 1), 1, false);
        aspace = GradedSpace.new(SU2(3), 1, false);
        
        Sleft = Tensor.ones(pspace, [pspace aspace]);
        Sright = -Tensor.ones([aspace pspace], pspace);
        
        S = (spin^2 + spin) * contract(Sleft, [-1 1 -4], Sright, [1 -2 -3], 'Rank', [2 2]);
        
    otherwise
        error('models:TBA', 'not implemented');
end
    
end

