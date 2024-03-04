function E0 = quantum1dIsing_energy(J, h)
% Compute the groundstate energy of the transverse field Ising model in 1d.
%
%

E0 = -integral(@(k) quantum1dIsing_dispersion(k, 'J', J, 'h', h), 0, pi) / (2 * pi);

end

