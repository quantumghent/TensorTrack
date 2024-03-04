function ind = sub2ind_(sz, I)
%sub2ind_ faster implementation of builtin sub2ind.

if isempty(I)
    ind = [];
    return;
elseif size(I,2) == 1
    ind = I;
    return;
end

numOfIndInput = size(I,2);

ind = I(:,1);
if numOfIndInput >= 2
    %Compute linear indices
    ind = ind + (I(:,2) - 1).*sz(1);
end 
    
if numOfIndInput > 2
    %Compute linear indices
    k = cumprod(sz);
    for i = 3:numOfIndInput
        ind = ind + (I(:,i)-1)*k(i-1);
    end
end

end