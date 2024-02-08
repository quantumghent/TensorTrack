function mustBeSorted(A)
% Validate that input is sorted.
% 
% :code:`mustBeSorted(A)` throws an error if :code:`A` is not sorted.
%
% Note
% ----
% Supported by all classes that define these methods: :code:`issorted`

if ~issorted(A)
    throwAsCaller(createValidatorException('TensorTrack:validators:mustBeSorted'));
end

end
