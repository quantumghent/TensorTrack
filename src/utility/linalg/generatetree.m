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
