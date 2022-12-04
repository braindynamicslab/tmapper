function A = weightedAdj(g)
% WEIGHTEDAJD get the weighted adjacency matrix from a graph object "g".
% This is to accomandate for over versions of matlab. 
% A = weightedAdj(g) input a graph g and get weighted adjacency matrix A. 
%{
created by MZ, 9/11/2019
%}
    if ~isfield(g.Edges,'Weight')
        g.Edges.Weight = ones(height(g.Edges),1);
    end
    try
        A = adjacency(g,'weighted');
    catch me
        A = adjacency(g);
        nn = numnodes(g);
        [s,t] = findedge(g);
        A = sparse(s,t,g.Edges.Weight,nn,nn);
    end
end