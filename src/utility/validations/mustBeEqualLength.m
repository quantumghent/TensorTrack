function mustBeEqualLength(a, b)
% mustBeEqualLength - Validate that the inputs are of equal length.
%   mustBeEqualLength(a, b) throws an error if length(a) ~= length(b)
%
%   Class support:
%       All classes that define these methods:
%           length

if length(a) ~= length(b)
    throwAsCaller(createValidatorException('TensorTrack:validators:mustBeEqualLength'));
end

end
