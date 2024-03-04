function V = fZ2Space(charges, bonds, isdual)
% Convenience constructor for :math:`fZ_2`-graded spaces.
%
% Arguments
% --------
% charges : (1, :) :class:`logical`-like
%   vector of :math:`fZ_2` charge labels.
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

V = GradedSpace.new(fZ2(charges), bonds, isdual);

end
