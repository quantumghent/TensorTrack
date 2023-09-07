%% Groundstate search
alg = Vumps('which', 'smallestreal', 'maxiter', 100);

symmetry = 'U1';

mpo = quantum1dHeisenberg('Spin', 1, 'Symmetry', symmetry);

switch symmetry
    case 'SU2'
        vspace = GradedSpace.new(SU2(2:2:6), [5 5 1], false);
    case 'U1'
        vspace = GradedSpace.new(U1(-3:2:3), [2 5 5 2], false);
    otherwise
        error('invalid symmetry');
end

mps = initialize_mps(mpo, vspace);

if ~exist('gs_mps', 'var')
    [gs_mps] = fixedpoint(alg, mpo, mps);
end

lambda = expectation_value(gs_mps, mpo);
assert(isapprox(lambda, -1.401, 'RelTol', 1e-2));

%% Excitations
p = pi;
charge = U1(0);
qp = InfQP.randnc(gs_mps, gs_mps, p, charge);
tic;
[qp, mu] = excitations(QPAnsatz(), mpo, qp);
toc;
mu
