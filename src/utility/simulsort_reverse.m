function [I, varargout] = simulsort_reverse(arrays, kwargs)
% Simultaneous sorting of several input arrays.
%
% Sorts the elements such that equal elements of array{i+1} appear sorted by
% array{i}.
%
% This is achieved by sorting from front to end, making use of the fact that SORT is stable
% and thus will not mess up the order of earlier elements when later elements are equal.
%
% Usage
% -----
% :code:`[I, array1, array2, ...] = simulsort(array1, array2, ..., kwargs)`
%
% Arguments
% ---------
% array1, array2, ...
%   arrays of equal size that need to be sorted. These can be of any type that supports
%   SORT.
%
% Keyword Arguments
% -----------------
% Dimension : :class:`int`
%   determine the dimension along which to sort. This behaves similarly to SORT, by default:
%   - for vectors, sorts the elements
%   - for matrices, sorts each column
%   - for N-D arrays, sorts along the first non-singleton dimension.
%
% Direction : :class:`char`, 'ascend' or 'descend'
%   specify the sorting direction, defaults to 'ascend'.
%
% Returns
% -------
% I : :class:`int`
%   permutation vector that brings the input arrays into sorted order.
%
% array1, array2, ...
%   sorted arrays.

arguments (Repeating)
    arrays
end

arguments
    kwargs.Dimension = find(size(arrays{1}) ~= 1, 1)
    kwargs.Direction {mustBeMember(kwargs.Direction, {'ascend', 'descend'})} = 'ascend'
end


%% Input validation
sz = size(arrays{1});
assert(length(sz) == 2, 'Not implemented for N-D arrays yet.');
for i = 2:length(arrays)
    assert(all(sz == size(arrays{i})), 'Input arrays should have equal size.');
end


%% Special cases
if all(sz == 1)
    varargout = arrays;
    I = 1;
    return
end
if any(sz == 0)
    varargout = arrays;
    I = [];
    return
end


%% Sort first array
[varargout{1}, I] = sort(arrays{1}, kwargs.Dimension, kwargs.Direction);
if kwargs.Dimension == 1
    for k = 2:length(arrays)
        for i = 1:size(I, 2)
            varargout{k}(:, i) = arrays{k}(I(:, i), i);
        end
    end
else
    for k = 2:length(arrays)
        for i = 1:size(I, 1)
            varargout{k}(i, :) = arrays{k}(i, I(i, :));
        end
    end
end

%% sort other arrays
for n = 2:length(arrays)
    [varargout{n}, I_] = sort(varargout{n}, kwargs.Dimension, kwargs.Direction);
    
    if kwargs.Dimension == 1
        % permute arrays
        for k = 1:length(arrays)
            if k == n, continue; end
            for i = 1:size(I, 2)
                varargout{k}(:, i) = varargout{k}(I_(:, i), i);
            end    
        end
        
        % permute index
        for i = 1:size(I, 2)
            I(:, i) = I(I_(:, i), i);
        end
        
    else
        % permute arrays
        for k = 1:length(arrays)
            if k == n, continue; end
            for i = 1:size(I, 1)
                varargout{k}(i, :) = varargout{k}(i, I_(i, :));
            end
        end
        
        % permute index
        for i = 1:size(I, 1)
            I(i, :) = I(i, I_(i, :));
        end
    end
end

end
