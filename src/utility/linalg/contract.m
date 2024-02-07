function C = contract(tensors, indices, kwargs)
% Compute a tensor network contraction.
%
% Usage
% -----
% :code:`C = contract(t1, idx1, t2, idx2, ...)`
%
% :code:`C = contract(..., 'Conj', conjlist)`
%
% :code:`C = contract(..., 'Rank', r)`
%
% Repeating Arguments
% -------------------
% tensors : :class:`.Tensor`
%   list of tensors that constitute the vertices of the network.
%
% indices : (1, :) :class:`int`
%   list of indices that define the links and contraction order, using `ncon-like syntax <https://arxiv.org/abs/1402.0939>`_.
%
% Keyword Arguments
% -----------------
% Conj : (1, :) :class:`logical`
%   optional list to flag that tensors should be conjugated.
%
% Rank : (1, 2) :class:`int`
%   optionally specify the rank of the resulting tensor.
%
% Returns
% -------
% C : :class:`.Tensor` or :class:`numeric`
%   result of the tensor network contraction.

% TODO contraction order checker, order specifier.

arguments (Repeating)
    tensors
    indices (1, :) {mustBeInteger}
end

arguments
    kwargs.Conj (1, :) logical = false(size(tensors))
    kwargs.Rank = []
    kwargs.Debug = false
    kwargs.CheckOptimal = false
end

assert(length(kwargs.Conj) == length(tensors));

if kwargs.CheckOptimal
    legcosts = zeros(2, 0);
    for i = 1:length(indices)
        legcosts = [legcosts [indices{i}; size(tensors{i}, 1:length(indices{i}))]];
    end
    legcosts = unique(legcosts.', 'rows');
    
    currentcost = contractcost(indices, legcosts);
    [sequence, cost] = netcon(indices, 1, 1, currentcost, 1, legcosts);
    
    if cost < currentcost
        warning('suboptimal contraction order.\n optimal: %s', ...
            num2str(sequence));
    end
end

for i = 1:length(tensors)
    [i1, i2] = traceinds(indices{i});
    tensors{i} = tensortrace(tensors{i}, i1, i2);
    indices{i}([i1 i2]) = [];
end

debug = kwargs.Debug;

% Special case for single input tensor
if nargin == 2
    C = tensors{1};
    if kwargs.Conj, C = conj(C); end
    
    if ~isempty(indices{1})
        [~, order] = sort(indices{1}, 'descend');
        C = permute(C, order);
    end
    
    return
end

% Generate trees
contractindices = cellfun(@(x) x(x > 0), indices, 'UniformOutput', false);
partialtrees = num2cell(1:length(tensors));
tree = generatetree(partialtrees, contractindices);

% contract all subtrees
[A, ia, ca] = contracttree(tensors, indices, kwargs.Conj, tree{1}, debug);
[B, ib, cb] = contracttree(tensors, indices, kwargs.Conj, tree{2}, debug);

% contract last pair
[dimA, dimB] = contractinds(ia, ib);

if debug, contractcheck(A, ia, ca, B, ib, cb); end

C = tensorprod(A, B, dimA, dimB, ca, cb, 'NumDimensionsA', length(ia));
ia(dimA) = [];  ib(dimB) = [];
ic = [ia ib];

% permute last tensor
if ~isempty(ic) && length(ic) > 1
    [~, order] = sort(ic, 'descend');
    C = permute(C, order);
end

end
