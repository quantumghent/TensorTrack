function c = c_min(kwargs)
% Fermionic annihilation operator.
%
% Keyword arguments
% -----------------
% 'Side' : :class:`char`
%   side, 'left or 'right'.
%
% Returns
% c : :class:`.Tensor`
%   annihilation operator represented as a 3-leg tensor with :math:`fZ_2` symmetry.

arguments
    kwargs.side = 'left'
end

pspace = fZ2Space([0 1], [1 1], false);
vspace = fZ2Space(1, 1, false);

switch kwargs.side
    case 'right'
        c = Tensor.zeros([vspace pspace], pspace);
        c = fill_matrix(c, {0 1});
        
    case 'left'
        c = Tensor.zeros(pspace, [pspace vspace]);
        c = fill_matrix(c, {1 0});
        
    otherwise
        error('models:argerror', 'invalid side')
end

end

