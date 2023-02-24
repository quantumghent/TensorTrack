function cost = contractcost(indices, legCosts)

cost = 0;

allInds = horzcat(indices{:});
numCont = max(allInds);

table = zeros(numCont, 2);
for ii = 1:length(indices)
    for jj = indices{ii}(indices{ii} > 0)
        if table(jj, 1) == 0
            table(jj, 1) = ii;
        else
            table(jj, 2) = ii;
        end
    end
end


%% Do contractions
ctr = 1;
contlist = 1:numCont;

while ~isempty(contlist)
    i1 = table(contlist(1), 1);
    i2 = table(contlist(1), 2);
    
    if i1 + i2 == 0
        contlist(1) = [];
        continue;
    else
        assert(i1 ~= i2, 'Tensor:isOptimalContract', 'Traces are not implemented.');
    end
    
    labels1 = indices{i1};
    labels2 = indices{i2};
    
    [pos1, pos2] = contractinds(labels1, labels2);
    unc1 = 1:length(labels1); unc1(pos1) = [];
    unc2 = 1:length(labels2); unc2(pos2) = [];
    contracting = labels1(pos1);
    
    % cost is dim(unc1) * dim(pos1) * dim(unc2)
    dims1 = zeros(1, length(unc1));
    for ii = 1:length(unc1)
        label = labels1(unc1(ii));
        dims1(ii) = legCosts(legCosts(:, 1) == label, 2);
    end
    dims2 = zeros(1, length(pos1));
    for ii = 1:length(pos1)
        label = labels1(pos1(ii));
        dims2(ii) = legCosts(legCosts(:, 1) == label, 2);
    end
    dims3 = zeros(1, length(unc2));
    for ii = 1:length(unc2)
        label = labels2(unc2(ii));
        dims3(ii) = legCosts(legCosts(:, 1) == label, 2);
    end
    
    cost = cost + prod([dims1 dims2 dims3]);
    
    % update remaining contractions
    indices{i1} = [indices{i1}(unc1) indices{i2}(unc2)];
    indices(i2) = [];
    table(table == i2) = i1;
    table(table > i2) = table(table > i2) - 1;
    contlist = contlist(~ismember(contlist, contracting));
    
    ctr = ctr + 1;
end


end
