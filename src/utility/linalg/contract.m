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
% tensors : :class:`Tensor`
%   list of tensors that constitute the vertices of the network.
%
% indices : int
%   list of indices that define the links and contraction order, using ncon-like syntax.
%
% Keyword Arguments
% -----------------
% Conj : (1, :) logical
%   optional list to flag that tensors should be conjugated.
%
% Rank : (1, 2) int
%   optionally specify the rank of the resulting tensor.
%
% Returns
% -------
% C : :class:`Tensor` or numeric
%   result of the tensor network contraction.

% TODO contraction order checker, order specifier.

arguments (Repeating)
    tensors
    indices (1, :) {mustBeInteger}
end

arguments
    kwargs.Conj (1, :) logical = false(size(tensors))
    kwargs.Rank = []
end

assert(length(kwargs.Conj) == length(tensors));

for i = 1:length(tensors)
    if length(indices{i}) > 1
        assert(length(unique(indices{i})) == length(indices{i}), ...
            'Tensors:TBA', 'Traces not implemented.');
    end
end

% Special case for single input tensor
if nargin == 2
    [~, order] = sort(indices{1}, 'descend');
    C = tensors{1};
    if isnumeric(C)
        C = permute(C, order);
        if kwargs.Conj
            C = conj(C);
        end
    else
        if kwargs.Conj
            C = permute(C', order(length(order):-1:1), kwargs.Rank);
        else
            C = permute(C, order, kwargs.Rank);
        end
    end
    return
end

% Generate trees
contractindices = cellfun(@(x) x(x > 0), indices, 'UniformOutput', false);
partialtrees = num2cell(1:length(tensors));
tree = generatetree(partialtrees, contractindices);

% contract all subtrees
[A, ia, ca] = contracttree(tensors, indices, kwargs.Conj, tree{1});
[B, ib, cb] = contracttree(tensors, indices, kwargs.Conj, tree{2});

% contract last pair
[dimA, dimB] = contractinds(ia, ib);
C = tensorprod(A, B, dimA, dimB, ca, cb, 'NumDimensionsA', length(ia));
ia(dimA) = [];  ib(dimB) = [];
ic = [ia ib];

% permute last tensor
if ~isempty(ic) && length(ic) > 1
    [~, order] = sort(ic, 'descend');
    if isnumeric(C)
        C = permute(C, order);
    else
        if isempty(kwargs.Rank)
            kwargs.Rank = [length(order) 0];
        end
        C = permute(C, order, kwargs.Rank);
    end
end

end

function tree = generatetree(partialtrees, contractindices)
if length(partialtrees) == 1
    tree = partialtrees{1};
    return
end

if all(cellfun('isempty', contractindices)) % disconnected network
    partialtrees{end - 1} = partialtrees(end - 1:end);
    partialtrees(end) = [];
    contractindices(end) = [];
else
    tocontract = min(horzcat(contractindices{:}));
    tinds = find(cellfun(@(x) any(tocontract == x), contractindices));
    assert(length(tinds) == 2);
    partialtrees{tinds(1)} = partialtrees(tinds);
    partialtrees(tinds(2)) = [];
    contractindices{tinds(1)} = unique1(horzcat(contractindices{tinds}));
    contractindices(tinds(2)) = [];
end

tree = generatetree(partialtrees, contractindices);
end

function [C, ic, cc] = contracttree(tensors, indices, conjlist, tree)

if isnumeric(tree)
    C = tensors{tree};    
    ic = indices{tree};
    cc = conjlist(tree);
    return
end

[A, ia, ca] = contracttree(tensors, indices, conjlist, tree{1});
[B, ib, cb] = contracttree(tensors, indices, conjlist, tree{2});
[dimA, dimB] = contractinds(ia, ib);
C = tensorprod(A, B, dimA, dimB, ca, cb, 'NumDimensionsA', length(ia));

ia(dimA) = [];
ib(dimB) = [];
ic = [ia ib];
cc = false;
end
