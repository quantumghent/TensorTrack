function sub2 = sub2sub(sz2, sz1, sub1)
%SUB2SUB Convert subscripts to subscripts corresponding to different array
%dimensions.

if isempty(sub1) || isequal(sz1, sz2)
    sub2 = sub1;
    return;
end
sub2 = zeros(size(sub1, 1), numel(sz2)); % preallocate new subs
nto1 = sum(cumsum(flip(sz1-1), 2) == 0, 2); % extract number of trailing ones in both sizes
nto2 = sum(cumsum(flip(sz2-1), 2) == 0, 2);
pos1_prev = 0;  pos2_prev = 0;
flag = true;
while flag
    [pos1, pos2] = find(cumprod(sz1(pos1_prev+1:end)).' == cumprod(sz2(pos2_prev+1:end)), 1);
    pos1 = pos1 + pos1_prev;
    pos2 = pos2 + pos2_prev;
    if prod(sz1(pos1_prev+1:pos1)) > 2^48-1
        error('Cannot map subscripts to new size as intermediate index exceeds MAXSIZE')
    end
    sub2(:, pos2_prev+1:pos2) = ind2sub_(sz2(pos2_prev+1:pos2), sub2ind_(sz1(pos1_prev+1:pos1), sub1(:, pos1_prev+1:pos1)));
    if      (isempty(pos2) && numel(sz2) - nto2 == 0) || ...
            (isempty(pos1) && numel(sz1) - nto1 == 0) || ...
            pos2 == numel(sz2) - nto2 || ...
            pos1 == numel(sz1) - nto1
        flag = false;
    else
        pos1_prev = pos1;
        pos2_prev = pos2;
    end
end
sub2(:, end-nto2+1:end) = 1;

end
