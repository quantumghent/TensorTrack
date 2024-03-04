function bool = iscircperm(p)

bool = true;
for i = 1:length(p)
    if isequal(circshift(p, i), 1:length(p))
        return
    end
end
bool = false;

end

