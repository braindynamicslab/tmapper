function [allbd,pclusters_nodes,pclusters_interior,pclusters_boundary,...
    pcluster_conn_dir,pclusters_intcrtpts, allupath, Tp] ...
    = CyclePathDecomp(dg,varargin)
%CYCLEPATHDECOMP Cycle-Path decomposition of a directed graph. Cycles on
%the digraph is first computed. Cycles are then clustered such that cycles
%with `clusterthres` amount of overlap belong to the same cluster.
%Boundaries between clusters of cycles all calculated, i.e. nodes where one
%cycle enters into another. Boundary nodes are further used to cut cycles
%into paths. The function returns the final results of path clustering. 
% 
%   [allbd,pclusters_nodes,pclusters_interior,pclusters_boundary,...
%     pcluster_conn_dir,pclusters_intcrtpts, allupath, Tp] ...
%     = CyclePathDecomp(dg,varargin)
% input: 
%   dg: a digraph
% output:
%   allbd: n-by-1 vector, containing all boundary points between cycle
%   clusters and resulted path clusters.
%   pclusters_nodes: N-by-1 cell array for N path clusters. Each cell
%   contains the index of nodes in each path cluster.
%   pclusters_interior: N-by-1 cell array. Each cell contains the index of
%   nodes in the interior of each path cluster (not boundary).
%   pclusters_boundary: N-by-1 cell array. Each cell contains the index of
%   nodes that are the boundary of the path cluster.
%   pcluster_conn_dir: N-by-N cell array. Cell (ii,jj) contains the
%   boundary node connecting path cluster ii to path cluster jj. 
%   pclusters_intcrtpts: N-by-1 cell array. Each cell contains the index of
%   critical points (nodes with more than 1 sources or targets) interior to
%   the cluster. 
%   allupath: M-by-1 cell array. A list of all the unique path.
%   Tp: M-by-1 vector of cluster assignments of `allupath`.
% parameters:
%   clusterthres: threshold of overlap between cycles in the same cluster.
%   plotmat: plot overlap matrix, and cluster assignment for diagnostics.
%   Default true.
%   plotmds: plot 2D projection of cycles/paths using MDS. Default false.
%   plothist: plot linkage histograms. Default false.
%   reordermat: sort overlap matrix according to cluster assigments.
%   Default true.
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 9-3-2020 ~
%}

p=inputParser;
p.addParameter('clusterthres',0.5) % threshold for cycle/path clustering
p.addParameter('plotmat',true)
p.addParameter('plotmds',false)
p.addParameter('plothist',false)
p.addParameter('reordermat',true)
p.parse(varargin{:})
par=p.Results;

% -- find all the cycles
[~,~,~,allcycles]=CycleCount2p(dg.adjacency);
% -- cluster cycles and find boundaries
T=CycleCluster(allcycles,par.clusterthres,...
    'plotmat',par.plotmat, 'plotmds',par.plotmds,...
    'plothist',par.plothist, 'reordermat', par.reordermat);
% :: find the boundaries between cycle clusters
[~, ~, ...
    ~, clusters_boundary, ~,...
    ~, ~] ...
    = CycleClusterConn(dg, allcycles, T);
% :: plot blocks for visual check
if par.plotmat
    addDiagBlock(gca,sort(T),zeros(max(T),3));
end
disp(['# cycle clusters = ' num2str(max(T))])

% -- break cycles into path
allbd = unique(cell2mat(clusters_boundary));
allupath=Cycles2Paths(allcycles,allbd);
% :: clustering path
Tp=CycleCluster(allupath,par.clusterthres,...
    'plotmat',par.plotmat, 'plotmds',par.plotmds,...
    'plothist',par.plothist, 'reordermat', par.reordermat);
[~, pcluster_conn_dir, ...
    pclusters_nodes, pclusters_boundary, pclusters_interior,...
    ~, pclusters_intcrtpts] ...
    = CycleClusterConn(dg, allupath, Tp);
% :: plot blocks for visual check
if par.plotmat
    addDiagBlock(gca,sort(Tp),zeros(max(Tp),3));
end
disp(['# path clusters = ' num2str(max(Tp))])


end

