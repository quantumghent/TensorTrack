function mustBeEqualLength(a, b)
% Validate that the inputs are of equal length.
%   
% :code:`mustBeEqualLength(a, b)` throws an error if :code:`length(a) ~= length(b)`
%
% Note
% ----
% Supported by all classes that define these methods: :code:`length`

if length(a) ~= length(b)
    throwAsCaller(createValidatorException('TensorTrack:validators:mustBeEqualLength'));
end

end
