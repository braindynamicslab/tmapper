function mu = nodemeasure(nodemembers)
%NODEMEASURE count the number of members for each node and normalize the
%sum to 1
%   mu = nodemeasure(nodemembers)

mu = nodesize(nodemembers);
mu = mu/sum(mu);
end

