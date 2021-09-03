function [cluster_conn, cluster_conn_dir, ...
    clusters_nodes, clusters_boundary, clusters_interior,...
    clusters_crtpts, clusters_intcrtpts] = CycleClusterConn(dg, allcycles, clusterIdx)
% CYCLECLUSTERCONN connectivity between M clusters of N cycles.
%   [cluster_conn, cluster_conn_dir] = CycleClusterConn(dg, allcycles, clusterIdx)
% input:
%   dg: a digraph
%   allcycles: N-by-1 cell array, each cell contains the path of one cycle
%   clusterIdx: N-by-1 integer array, each integer is the index of a
%   loop-cluster. Integers are from 1 to M. 
% output:
%   cluster_conn: M-by-M cell array. cell (i,j) = cell (j,i) contain nodes
%   that are the boundary between loop clusters.
%   cluster_conn_dir: M-by-M cell array. cell (i,j) contain nodes that are
%   in the boundary of cluster j that receives a link from at least one
%   node from cluster i (excluding shared nodes). 
%   clusters_nodes: M-by-1 cell array. Each cell contains the indices of all
%   nodes in each cluster. 
%   clusters_boundary: M-by-1 cell array. Each cell the boundary nodes of
%   a cluster. 
%   clusters_interior: M-by-1 cell array. Each cell contains the nodes that
%   are not the boundary of the cluster. 
%   clusters_crtpts: M-by-1 cell array. Each cell contains the crtical
%   points in each cluster. Critical points are the nodes that have more
%   than one targets or sources. 
%   clusters_intcrtpts: M-by-1 cell array. Each cell contains the critial
%   points in each cluster that are not in the boundary.
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 9-1-2020 ~
modifications:
(9-3-2020) perserve order of nodes of clusters_nodes and clusters_interior
%}
clusters = unique(clusterIdx);
N_clusters = length(clusters);

% -- find all critical points (all nodes with more than 1 sources or
% targets)
outdgr = outdegree(dg);
indgr = indegree(dg);
crtpts = find(outdgr>1 | indgr>1);

% -- find and classify critical points in each loop-cluster
[clusters_nodes,clusters_crtpts,clusters_intcrtpts,...
    clusters_inneighbors, clusters_outneighbors,...
    clusters_boundary, clusters_interior]...
    = deal(cell(N_clusters,1));

for ii=1:N_clusters
    clusters_nodes{ii}=unique(cell2mat(allcycles(clusterIdx==ii)'),'stable');
    clusters_crtpts{ii}=intersect(unique(cell2mat(allcycles(clusterIdx==ii)')),crtpts);
    
    % -- find the neighbors of critical points
    N_cp = length(clusters_crtpts{ii});
    innbg = [];
    outnbg = [];
    for n=1:N_cp
        [~,innbgidx]=dg.inedges(clusters_crtpts{ii}(n));
        [~, outnbgidx]=dg.outedges(clusters_crtpts{ii}(n));
        innbg=[innbg; innbgidx];
        outnbg=[outnbg; outnbgidx];
    end
    % :: neighbors of the cluster
    innbg = setdiff(innbg,clusters_nodes{ii});
    outnbg= setdiff(outnbg,clusters_nodes{ii});
    
    % -- find boundaries
    bd =[];
    if ~isempty(innbg)
        for n=innbg'
            [~, bdidx]=dg.outedges(n);
            bd = [bd; bdidx];
        end
    end
    if ~isempty(outnbg)
        for n=outnbg'
            [~, bdidx]=dg.inedges(n);
            bd = [bd; bdidx];
        end
    end
    bd=intersect(unique(bd),clusters_nodes{ii});
    
    % -- store results
    clusters_inneighbors{ii}=innbg;
    clusters_outneighbors{ii}=outnbg;
    clusters_boundary{ii}=bd;
    clusters_interior{ii}=setdiff(clusters_nodes{ii},bd,'stable');
    clusters_intcrtpts{ii}=setdiff(clusters_crtpts{ii},bd);
end

% -- find critical points shared between loop-clusters
[cluster_conn,cluster_conn_dir,cluster_conn_in,cluster_conn_out] ...
    = deal(cell(N_clusters));

for ii=1:N_clusters
    for jj=1:N_clusters
        % :: undirected connectivity (boundary nodes shared between ii and jj)
        cluster_conn{ii,jj} = intersect(clusters_boundary{ii},clusters_nodes{jj});
        % :: nodes in ii going into jj
        cluster_conn_in{ii,jj} = intersect(clusters_inneighbors{jj},clusters_nodes{ii});
        % :: nodes in jj going from ii
        cluster_conn_out{ii,jj} = intersect(clusters_outneighbors{ii},clusters_nodes{jj});
    end
end

% :: directed connectivity (nodes where ii enters jj)
for ii=1:N_clusters
    for jj=1:N_clusters
        % :: boundary from cluster ii to cluster jj
        bd=[];
        if ~isempty(cluster_conn{ii,jj})
            for nn = cluster_conn{ii,jj}'
                [~, bdidx] = dg.inedges(nn);
                if ~isempty(intersect(bdidx, cluster_conn_in{ii,jj}))
                    bd = [bd; nn];
                end
            end
        end
        
        cluster_conn_dir{ii,jj} = intersect(unique(bd),cluster_conn{ii,jj});
    end
end
end