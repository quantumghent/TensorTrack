function mustBeUnique(A)
% mustBeUnique - Validate that there are no repeated values.
%  mustBeUnique(A) throws an error if A contains repeated values.
%
%     Class support:
%        All classes that define these methods:
%           unique
%           length

if length(A) ~= length(unique(A))
   throwAsCaller(createValidatorException('validators:mustBeUnique'));
end

end
