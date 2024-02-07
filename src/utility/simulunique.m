function varargout = simulunique(varargin)
% Set unique over multiple arrays.


nlhs = max(1, nargout);
narginchk(1, 4);

nrhs = nargin;

% acceptable combinations, with optional inputs denoted in []:
% unique(A..., ['rows'], ['first'/'last'])
% unique(A..., ['rows'], ['sorted'/'stable'])
% where the optional arguments may have reversed order.

flagvals = ["rows" "first" "last" "sorted" "stable"];
flaginds = false(size(flagvals));

i = nrhs;
while i > 1
    flag = varargin{i};
    if ~isstring(flag) && ~ischar(flag), break; end
    
    foundflag = startsWith(flagvals, flag, 'IgnoreCase', true);
    if sum(foundflag) ~= 1
        error('utility:simulunique:argerror', 'Unknown flag.');
    end
    if flaginds(foundflag)
        error('utility:simulunique:argerror', 'Repeated flag.');
    end
    flaginds(foundflag) = true;
    i = i - 1;
end
if sum(flaginds(2:5)) == 0, flaginds(4) = true; end

first = flaginds(2);    last = flaginds(3);
assert(~(first && last), 'utility:simulunique:argerror', 'Flag combination not allowed.');
sorted = flaginds(4);   stable = flaginds(5);
assert(~(sorted && stable), ...
    'utility:simulunique:argerror', 'Flag combination not allowed,');
assert(~((sorted || stable) && (first || last)), ...
    'utility:simulunique:argerror', 'Flag combination not allowed.');

if flaginds(1)
    [varargout{1:nlhs}] = simuluniquerows(varargin(1:i), first, last, sorted, stable);
else
    [varargout{1:nlhs}] = simuluniqueels(varargin(1:i), first, last, sorted, stable);
end

end

function [C, indA, indC] = simuluniquerows(A, first, last, sorted, stable)

assert(all(cellfun(@ismatrix, A)), 'utility:simulunique:argerror', 'Input not matrices.');

numRows = size(A{1}, 1);
numCols = sum(cellfun(@(x) size(x, 2), A));

sortA = cell(size(A));
[indSortA, sortA{:}] = simulsortrows(A{:});

groupsSortA = any(sortA{1}(1:numRows-1,:) ~= sortA{1}(2:numRows, :), 2);
for i = 2:length(A)
    groupsSortA = groupsSortA | ...
        any(sortA{i}(1:numRows-1,:) ~= sortA{i}(2:numRows,:), 2);
end

if numRows ~= 0
    if last
        groupsSortA = [groupsSortA; true];
    else
        groupsSortA = [true; groupsSortA];
    end
end

% create output
C = cell(size(A));
if stable
    invIndSortA = invperm(indSortA);
    logIndA = groupsSortA(invIndSortA);
    for i = 1:length(A)
        C{i} = A{i}(logIndA, :);
    end
else
    for i = 1:length(A)
        C{i} = sortA{i}(groupsSortA, :);
    end
end

% find indA
if nargout > 1
    if stable
        indA = find(logIndA);
    else
        indA = indSortA(groupsSortA);
    end
end

% find indC
if nargout == 3
    if last
        if numRows == 0
            indC = cumsum(full(groupsSortA));
        else
            indC = cumsum([1; full(groupsSortA(1:end-1))]);
        end
        indC(indSortA) = indC;
    elseif sorted
        indC = cumsum(full(groupsSortA));
        indC(indSortA) = indC;
    else % stable
        if numCols == 0
            indC = ones(numRows, 1);
        elseif numRows == 0
            indC = zeros(0, 1);
        else
            indSortC = simulsortrows(C{:});
            
            lengthGroupSortA = diff(find([groupsSortA; true]));
            
            diffIndSortC = diff(indSortC);
            diffIndSortC = [indSortC(1); diffIndSortC];
            
            indLengthGroupSortA = cumsum([1; lengthGroupSortA]);
            indLengthGroupSortA(end) = [];
            
            indCOrderedBySortA(indLengthGroupSortA, 1) = diffIndSortC;
            
            if sum(lengthGroupSortA) ~= length(indCOrderedBySortA)
                indCOrderedBySortA(sum(lengthGroupSortA), 1) = 0;
            end
            
            indCOrderedBySortA = cumsum(indCOrderedBySortA);
            indC = indCOrderedBySortA(invIndSortA);
        end
    end
end
end


function simuluniqueels(arrays, first, last, sorted, stable)
error('TBA');
end