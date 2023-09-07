function V = fZ2Space(charges, bonds, isdual)
% Convenience constructor for fZ2-graded spaces.
arguments
    charges
    bonds
    isdual = false
end

V = GradedSpace.new(fZ2(charges), bonds, isdual);

end
