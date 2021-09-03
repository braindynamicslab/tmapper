function g = digraph2graph(dg)
%DIGRAPH2GRAPH convert a directed graph, to a undirected graph, where the
%weight of each edge of the latter is the average of the two edges between
%two corresponding nodes. 
%   g = digraph2graph(dg)
%{
created by MZ, Aug 2019.
Modifications:
(9/11/2019) modified computation of adjacency matrix to accomandate older
versions of matlab. added failsafe for graph without nodenames.

%}
if ~isfield(dg.Edges,'Weight')
    dg.Edges.Weight = ones(height(dg.Edges),1);
end

A = weightedAdj(dg);

if isfield(dg.Nodes, 'Name')
    g = graph((A + A')/2,dg.Nodes.Name);
else
    g = graph((A + A')/2);
end
end

