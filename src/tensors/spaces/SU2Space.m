function V = SU2Space(charges, bonds, isdual)
% Convenience constructor for :math:`\mathrm{SU}(2)`-graded spaces.
%
% Arguments
% --------
% charges : (1, :) :class:`int`-like
%   vector of :math:`\mathrm{SU}(2)` charge labels.
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

V = GradedSpace.new(SU2(charges), bonds, isdual);

end
