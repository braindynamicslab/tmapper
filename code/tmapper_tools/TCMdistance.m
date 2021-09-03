function [TCM_pl] = TCMdistance(g,nodet,weighted)
%TCMDISTANCE calculate a time connectivity matrix where the edge is the
%minimal path length between two sample points.
%   [TCM_pl] = TCMdistance(g,nodet,weighted)
%   input:
%       g: matlab graph
%       nodet: members of each node as indices of time points, a cell array
%       weighted: whether to use the weight of the graph, default = false,
%       use as unweighted graph.
%   output:
%       TCM_pl: temporal connectivity using shortest path length on the
%       shape graph.
%{
created by MZ, 04-27-2019
modifications:
(8-23-2019) debug weight assignment -- no longer require input graph to be
weighted. also accomandate situation where t0~=0.
%}
if nargin<3 || isempty(weighted)
    weighted = false;
end
if ~weighted
    g.Edges.Weight=ones(height(g.Edges),1);
end
% -- extract time
t = unique(cell2mat(nodet));
Nt = length(t);
t_0 = min(t);

% -- construct TCM from graph distance
TCM_pl = NaN(Nt,Nt);

distmat = distances(g);

% -- assign diagonal elements to 0
for i=1:g.numnodes
        TCM_pl(nodet{i}-t_0+1,nodet{i}-t_0+1) = 0;
end

% -- assign off diagonal elements
for i=1:g.numnodes
    for j=i+1:g.numnodes
        TCM_pl(nodet{i}-t_0+1,nodet{j}-t_0+1) = nanmin(TCM_pl(nodet{i}-t_0+1,nodet{j}-t_0+1),distmat(i,j));
        TCM_pl(nodet{j}-t_0+1,nodet{i}-t_0+1) = nanmin(TCM_pl(nodet{j}-t_0+1,nodet{i}-t_0+1),distmat(j,i));
    end
end


end

