function [I, varargout] = simulsortrows(arrays, kwargs)
% Simultaneous sorting of several input arrays by row.
%   Sorts the rows such that equal rows of array{i} appear sorted by the rows of array{i+1}.
%
% This is achieved by sorting rows from end to front, making use of the fact that SORTROWS
% is stable and thus will not mess up the order of later rows when earlier rows are equal.
%
% Usage
% -----
% [I, array1, array2, ...] = simulsortrows(array1, array2, ..., kwargs)
%
% Arguments
% ---------
% array1, array2, ...
%   arrays of equal size that need to be sorted. These can be of any type that supports
%   SORTROWS.
%
% Keyword Arguments
% -----------------
% Col : int
%   vector of indices that specifies the columns used for sorting.
%
% Direction : 'ascend' or 'descend'
%   specify the sorting direction. You can also specify a different direction for each
%   column by using a cell array of 'ascend' and 'descend' of the same size as Col, such
%   that corresponding elements are sorted ascending or descending.
%
% Returns
% -------
% I : int
%   permutation vector that brings the input arrays into rowsorted order.
%
% array1, array2, ...
%   rowsorted arrays

arguments (Repeating)
    arrays
end

arguments
    kwargs.Col = []
    kwargs.Direction = 'ascend'
end

%% Input validation
sz = size(arrays{1});
assert(length(sz) == 2, 'Input should be matrices.');
for i = 2:length(arrays)
    assert(all(sz(1) == size(arrays{i}, 1)), 'Input arrays should have equal sizes.');
end


%% Sort with columns
if ~isempty(kwargs.Col)
    cols = cellfun(@(x) size(x, 2), arrays);
    arrayinds = zeros(1, length(cols));
    ctr = 0;
    for i = 1:length(arrays)
        arrayinds(ctr+(1:cols(i))) = i;
        ctr = ctr + cols(i);
    end
    
    col = kwargs.Col(end);
    [varargout{arrayinds(col)}, I] = sortrows(arrays{arrayinds(col)}, ...
        col - sum(cols(1:arrayinds(col)-1)));
    for k = [length(arrays):-1:arrayinds(col)+1 arrayinds(col)-1:-1:1]
        varargout{k} = arrays{k}(I, :);
    end
    
    for col = kwargs.Col(end-1:-1:1)
        [varargout{arrayinds(col)}, I_] = sortrows(varargout{arrayinds(col)}, ...
            col - sum(cols(1:arrayinds(col)-1)));
        for k = 1:length(arrays)
            if k == arrayinds(col), continue; end
            varargout{k} = varargout{k}(I_, :);
        end
        I = I(I_);
    end
    return;
end
    

%% Sort without columns
[varargout{length(arrays)}, I] = sortrows(arrays{end}, kwargs.Direction);
for k = length(arrays)-1:-1:1
    varargout{k} = arrays{k}(I, :);
end

for n = length(arrays)-1:-1:1
    [varargout{n}, I_] = sortrows(varargout{n}, kwargs.Direction);
    for k = 1:length(arrays)
        if k == n, continue; end
        varargout{k} = varargout{k}(I_, :);
    end
    I = I(I_);
end


end
