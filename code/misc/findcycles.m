function [cycles, lenCycles] = findcycles(G)
% FINDCYCLES find cycles in a directed graph based on adjacency matrix G.
% [cycles, lenCycles] = findcycles(G)
%   input:
%       G: adjacency matrix. 
%   output:
%       cycles: N-by-1 cell array, where N is the number of cycles, each
%       cell contains a 1-by-M_n vector, where M_n is the number of nodes
%       in the n-th cycle.
%       lenCycles: N-by-1 vector, the n-th element is the the number of
%       nodes in the n-th cycle.
%{
Author: Mengsen Zhang ~ <mengsenzhang@gmail.com> ~ 3/20/2020 
modifications:

%}


% -- adjacency matrix must be sparse
if ~issparse(G)
    G = sparse(G);
end

% -- list of nodes
numNodes = size(G,1); 
nodes = 1:numNodes;
% -- list of cycles
cycles = cell(numNodes,1);
N_cycles = 0;

% -- find cycles starting from each node
while ~isempty(nodes)
    origin = nodes(1);
    % -- spanning tree from the first node on the list 
   [D,P]=graphtraverse(G,origin);
   for d = D
       if G(d,origin)% if loops back to the node
           N_cycles = N_cycles + 1;
           % -- add path to list of cycles
           cycles{N_cycles} = graphpred2path(P,d);
           % -- remove nodes in the cycle from the list of nodes
           nodes = setdiff(nodes, cycles{N_cycles});
       end
   end
end

cycles(N_cycles+1:end)=[];

% -- length of the cycles
lenCycles = cell2mat(cellfun(@length, cycles, 'UniformOutput',0));
end