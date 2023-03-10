function V = Z2Space(charges, bonds, isdual)
% Convenience constructor for Z2-graded spaces.
arguments
    charges
    bonds
    isdual = false
end

V = GradedSpace.new(Z2(charges), bonds, isdual);

end
