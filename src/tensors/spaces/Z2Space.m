function V = Z2Space(charges, bonds, isdual)
% Convenience constructor for :math:`Z_2`-graded spaces.
%
% Arguments
% --------
% charges : (1, :) :class:`logical`-like
%   vector of :math:`Z_2` charge labels.
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

V = GradedSpace.new(Z2(charges), bonds, isdual);

end
