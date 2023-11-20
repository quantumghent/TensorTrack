function p = pauliterm(spin, i, j)

p = sqrt((spin + 1) * (i + j - 1) - i * j) / 2;

end

