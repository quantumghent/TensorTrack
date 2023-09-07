function E = quantum1dHubbard_energy(u)
% computes the Hubbard groundstate energy at half-filling
f = @(x) x.^(-1) .* besselj(0, x) .* besselj(1, x) ./ (1 + exp(2 * u * x));
I = integral(f, 0, Inf);
E = -(u + 4 * I);

end
