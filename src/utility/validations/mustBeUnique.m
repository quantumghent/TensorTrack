function mustBeUnique(A)
% Validate that there are no repeated values.
%
% :code:`mustBeUnique(A)` throws an error if :code:`A` contains repeated values.
%
% Note
% ----
% Supported all classes that define these methods: :code:`unique`, :code:`length`

if length(A) ~= length(unique(A))
   throwAsCaller(createValidatorException('validators:mustBeUnique'));
end

end
