function V = U1Space(charges, bonds, isdual)
% Convenience constructor for U1-graded spaces.
arguments
    charges
    bonds
    isdual = false
end

V = GradedSpace.new(U1(charges), bonds, isdual);

end
