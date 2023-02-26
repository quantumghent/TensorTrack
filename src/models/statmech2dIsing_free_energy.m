function f = statmech2dIsing_free_energy(beta)
% Compute the free energy of the classical 2d Ising model using Onsagers solution.

INTEGRAL_STEP_SIZE = 1e-6;

theta = 0:INTEGRAL_STEP_SIZE:(pi / 2);

sh = sinh(2 * beta);
ch = cosh(2 * beta);
x = 2 * sh / ch^2;

f = -1 / beta * (log(2 * cosh(2 * beta)) + ...
    1 / pi * trapz(theta, log(1/2 * (1 + sqrt(1 - x^2 * sin(theta).^2)))));

end
