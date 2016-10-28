function [b] = isInCell(set,cellarray, lengthsarrays)
%Is the set contained in the cell array
b=false;
k=length(set);
[~,idx]=find(lengthsarrays==k);
for i=1:length(idx)
    if (sum(set==cellarray{idx(i)})==k)
        b=true;
        break
    end
end

end

