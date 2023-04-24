function c_dagger = fcreation(kwargs)
arguments
    kwargs.side = 'left'
end


pspace = fZ2Space([0 1], [1 1], false);
vspace = fZ2Space(1, 1, false);

switch kwargs.side
    case 'right'
        c_dagger = Tensor.zeros([vspace pspace], pspace);
        c_dagger = fill_matrix(c_dagger, {1 0});
        
        
    case 'left'
        c_dagger = Tensor.zeros(pspace, [pspace vspace]);
        c_dagger = fill_matrix(c_dagger, {0 1});
        
    otherwise
        error('models:argerror', 'invalid side')
end

end

