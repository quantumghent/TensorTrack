function d = diracdelta(sz)
% Construct delta tensor of given size. Note that all dimensions must have the same length.

assert(all(sz == sz(1)))
d = zeros(sz);

subs = repmat(1:sz(1), length(sz), 1).';
idx = sub2ind_(sz, subs);

d(idx) = 1;

end

