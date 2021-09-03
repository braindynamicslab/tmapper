function n = nodesize(nodemembers)
%NODESIZE count the number of members for each node
%   n = nodesize(nodemembers)

n = cell2mat(cellfun(@(x) length(x),nodemembers,'UniformOutput',0));
end

