function V = U1Space(charges, bonds, isdual)
% Convenience constructor for :math:`\mathrm{U}(1)`-graded spaces.
%
% Arguments
% --------
% charges : (1, :) :class:`int`-like
%   vector of :math:`\mathrm{U}(1)` charge labels.
%
% bonds : (1, :) :class:`int`
%   degeneracy of each charge.
%
% isdual : :class:`logical`
%   indicate if the space is dual.
arguments
    charges
    bonds
    isdual = false
end

V = GradedSpace.new(U1(charges), bonds, isdual);

end
