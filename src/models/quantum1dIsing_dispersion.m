function e = quantum1dIsing_dispersion(k, kwargs)
% Compute the dispersion relation of the transverse field Ising model in 1d.
arguments
    k
    kwargs.J = 1.0
    kwargs.h = 0.5
end
g = 2 * kwargs.h;

e = kwargs.J * sqrt(g^2 + 1 - 2 * g * cos(k)) / 2;

end

