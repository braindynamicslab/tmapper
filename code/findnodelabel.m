function m = findnodelabel(members,x_label)
%CELLMODE find the most frequent labels for members in the same node.
%   m = findnodelabel(members,x_label) 
% input:
%   members: a N-by-1 cell array, each cell contains a vector of indices
%   x_label: label for each member. 
% output:
%   m: a N-by-1 vector each element is the most frequent label for a
%   particular cell.
%{
created by MZ, 9-13-2019
%}

m = cell2mat(cellfun(@(x) mode(x_label(x)), members,'uniformoutput',0));
end

