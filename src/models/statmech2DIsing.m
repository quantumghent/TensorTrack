function O = statmech2DIsing(kwargs)

arguments
    kwargs.beta = log(1 + sqrt(2)) / 2;
    kwargs.L = Inf     % size of system
    kwargs.Symmetry {mustBeMember(kwargs.Symmetry, {'Z1', 'Z2'})} = 'Z1'
end

if ~isfinite(kwargs.L)
    O = InfMpo({MpoTensor(bulk_mpo(kwargs.beta, kwargs.Symmetry))});
    
else
    r = boundary_mpo(kwargs.beta, kwargs.Symmetry);
    l = r';
    o = bulk_mpo(kwargs.beta, kwargs.Symmetry);
    
    O = FiniteMpo(MpsTensor(l), repmat({MpoTensor(o)}, 1, kwargs.L), MpsTensor(r));
end

end

function O = bulk_mpo(beta, symmetry)

if strcmp(symmetry, 'Z1')
    sz = [2 2 2 2];
    t = sqrtm([exp(beta) exp(-beta); exp(-beta) exp(beta)]);
    o = contract(diracdelta(sz), 1:4, t, [-1 1], t, [-2 2], t, [-3 3], t, [-4 4]);
    O = fill_tensor(Tensor.zeros(sz), o);
        
elseif strcmp(symmetry, 'Z2')
    s = GradedSpace.new(Z2(0, 1), [1 1], false);
    t = sqrtm(fill_matrix(Tensor(s, s), num2cell(2 * [cosh(beta) sinh(beta)])));
    o = Tensor.ones([s s], [s s]) / 2;
    
    O = contract(o, 1:4, t, [-1 1], t, [-2 2], t, [3 -3], t, [4 -4], 'Rank', [2 2]);
else
    error('invalid symmetry');
end

end

function O = boundary_mpo(beta, symmetry)

if strcmp(symmetry, 'Z1')
    sz = [2 2 2];
    t = sqrtm([exp(beta) exp(-beta); exp(-beta) exp(beta)]);
    o = contract(diracdelta(sz), 1:3, t, [-1 1], t, [-2 2], t, [-3 3]);
    O = fill_tensor(Tensor.zeros(sz), o);
        
elseif strcmp(symmetry, 'Z2')
    s = GradedSpace.new(Z2(0, 1), [1 1], false);
    t = sqrtm(fill_matrix(Tensor(s, s), num2cell(2 * [cosh(beta) sinh(beta)])));
    o = Tensor.ones([s s], s) / sqrt(2);
    
    O = contract(o, 1:3, t, [-1 1], t, [-2 2], t, [3 -3], 'Rank', [2 1]);
else
    error('invalid symmetry');
end

end