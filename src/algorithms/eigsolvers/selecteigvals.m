function select = selecteigvals(lambda, howmany, which)

switch which
    case {'largestabs', 'lm'}
        [~, p] = sort(abs(lambda), 'descend');
    case {'largestreal', 'lr'}
        [~, p] = sort(real(lambda), 'descend');
    case {'largestimag', 'li'}
        [~, p] = sort(imag(lambda), 'descend');
    case {'smallestabs', 'sm'}
        [~, p] = sort(abs(lambda), 'ascend');
    case {'smallestreal', 'sr'}
        [~, p] = sort(real(lambda), 'ascend');
    case {'smallestimag', 'si'}
        [~, p] = sort(imag(lambda), 'ascend');
end

select = p(1:howmany);

end