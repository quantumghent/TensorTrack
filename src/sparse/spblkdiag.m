function y = spblkdiag(varargin)
%SPBLKDIAG  Sparse block diagonal concatenation of matrix input arguments.
%
%                                     |A 0 .. 0|
%   Y = SPBLKDIAG(A,B,...)  produces  |0 B .. 0|
%                                     |0 0 ..  |

if nargin == 0
    y = [];
else
    checkAllAreMatrices(varargin); % Inputs must be all 2-dimensional
    X = varargin{1};
    clsX = class(X);
    if isobject(X) || ~isnumeric(X) || ~hasSameClass(varargin)
        y = X;
        if ~strcmp(clsX, 'logical')  %#ok<ISLOG>
           clsX = 'double';
        end
        for k = 2:nargin
            x = varargin{k};
            [p1,m1] = size(y);
            [p2,m2] = size(x);
            y = [y zeros(p1,m2,clsX); zeros(p2,m1,clsX) x]; %#ok
        end
    else      
        if hasSparse(varargin)
            y = matlab.internal.math.blkdiag(varargin{:}); % for sparse double
        else
            nz = sum(cellfun(@numel, varargin));
            [p, q] = getIndexVectors(varargin);
            y = spalloc(p(end), q(end), nz); %Preallocate
            for k = 1:nargin
                y(p(k)+1:p(k+1),q(k)+1:q(k+1)) = varargin{k};
            end
        end
    end
end

% Helper functions

% check if cellX contains all matrices.
function checkAllAreMatrices(cellX)
for i = 1:numel(cellX)
    if ~ismatrix(cellX{i})
        throwAsCaller(MException(message('MATLAB:blkdiag:inputMustBe2D')));
    end
end

% check if cellX contains all matrices of same class classA.
function tf = hasSameClass(cellX)
clsA = class(cellX{1});
for i = 2:numel(cellX)
    if ~strcmp(clsA, class(cellX{i}))
        tf = false;
        return
    end
end
tf = true;

% check if cellX contains all sparse matrices.
function tf = hasSparse(cellX)
for i = 1:numel(cellX)
    if issparse(cellX{i})
        tf = true;
        return
    end
end
tf = false;

% compute the index vector from each matrix.
function [p, q] = getIndexVectors(cellX)
numEl = numel(cellX);
p = zeros(1, numEl+1);
q = zeros(1, numEl+1);
for i = 1:numEl
    x = cellX{i};
    p(i+1) = size(x,1);
    q(i+1) = size(x,2);
end
p = cumsum(p);
q = cumsum(q);
