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
    kwargs.Col = 1:size(arrays{1}, 2)
    kwargs.Direction = 'ascend'
end

%% Input validation
sz = size(arrays{1});
assert(length(sz) == 2, 'Input should be matrices.');
for i = 2:length(arrays)
    assert(all(sz == size(arrays{i})), 'Input arrays should have equal sizes.');
end


%% Sort last array
[varargout{length(arrays)}, I] = sortrows(arrays{end}, kwargs.Col, kwargs.Direction);
for k = length(arrays)-1:-1:1
    varargout{k} = arrays{k}(I, :);
end


%% Sort other arrays
for n = length(arrays)-1:-1:1
    [varargout{n}, I_] = sortrows(varargout{n}, kwargs.Col, kwargs.Direction);
    for k = 1:length(arrays)
        if k == n, continue; end
        varargout{k} = varargout{k}(I_, :);
    end
    I = I(I_);
end

end
