function mustBeSorted(A)
% mustBeSorted - Validate that value is sorted.
%   mustBeSorted(A) throws an error if A is not sorted.
%
%   Class support:
%       All classes that define these methods:
%           issorted

if ~issorted(A)
    throwAsCaller(createValidatorException('TensorTrack:validators:mustBeSorted'));
end

end
