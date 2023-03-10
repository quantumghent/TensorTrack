function V = SU2Space(charges, bonds, isdual)
% Convenience constructor for SU2-graded spaces.
arguments
    charges
    bonds
    isdual = false
end

V = GradedSpace.new(SU2(charges), bonds, isdual);

end
