function x = between(x1,x,x2)

assert(x1 <= x2, 'range', 'x1 should be smaller than  or equal to x2');
if x < x1
    x = x1;
elseif x > x2
    x = x2;
end

end