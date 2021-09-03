function [h1, h2, cb, cb_, hg] = plotgraphtcm(g,x_label,t,nodemembers, varargin)
%PLOTGRAPHTCM plot a shape graph and its correspomding geodesic recurrence
%plot (a kind of temporal connectivity matrix). 
%   [h1,h2] = plotgraphtcm(g,x_label,t,nodemembers, ...)
% input:
%   g: a graph or digraph (MATLAB object). 
%   x_label: a label for each member of each node of the graph, assumed to
%   be a N-by-1 vector of integer indices, where N is the number of unique
%   members of all nodes.
%   t: time associated with each member of each node. a N-by-1 vector.
%   nodemembers: a numnodes-by-1 cell array. Each cell contains a vector of
%   integer indices, indicating which time points belong to this node.
% parameters:
%   nodesizerange: [s_min, s_max], where s_min is the smallest markersize
%   when plotting the nodes of the graph, and s_max is the largest
%   markersize.
%   colorlabel: label of the color axis (meaning of the color-coding of
%   nodes), default "attractor index".
% output:
%   h1: axis handle of the first subplot (the graph).
%   h2: axis handle of the second subplot (the geodesic recurrence plot).
%   cb: handle of the colorbar of graph
%   cb_: handle of the colorbar of recurrence plot
%   hg: handle of the graph
%{
created by MZ, 9/13/2019
modifications:
(9/30/2019) adjust color range, node size for special cases, and
calculation of distances (change to unweighted)
(2/11/2020) add colormap options
(11/11/2020) add more handle output
%}

p = inputParser;
p.addParameter('nodesizerange',[1 20]);
p.addParameter('colorlabel','attractor index')
p.addParameter('cmap','jet') % colormap
p.addParameter('normalize',false) % whether or not 
p.parse(varargin{:});

par = p.Results;

% -- check nodes
if nargin<4 || isempty(nodemembers)
    nodemembers = num2cell((1:g.numnodes)');
end

% -- check labels for members
if nargin<2 || isempty(x_label)
    x_label = ones(length(unique(cell2mat(nodemembers(:)))),1);
end

% -- define node size
nodesize = cell2mat(cellfun(@(x) length(x), nodemembers, 'UniformOutput',0));
bsinglemember = all(nodesize==1);% there is only one member associated with each node.
buniform = length(unique(nodesize))==1; % if all nodes are of the same size
if ~buniform%adjust nodesize with rank
    nodesize = rankval(nodesize);% the marker size reflects the rank of the node size
    nodesize = rescale(nodesize, min(par.nodesizerange), max(par.nodesizerange));
else%adjust nodesize with number of nodes
    nodesize = ones(g.numnodes,1)*(par.nodesizerange(1)+range(par.nodesizerange)/g.numnodes);
end

% -- define node labels
nodelabel = findnodelabel(nodemembers,x_label);

% -- plotting
figure('position',[10,10,1000,400]);
% plot graph
subplot(1,2,1)
hg=plot(g,'EdgeAlpha',0.3,'EdgeColor','k','NodeCData',nodelabel,'NodeLabel','',...
    'Layout','force','ArrowSize',5,'MarkerSize',nodesize); ;
axis equal
axis off
cb=colorbar;
cb.Label.String = par.colorlabel;
colormap(gca, par.cmap)
caxis([min(x_label) max(x_label)])
h1 = gca;

% plot geodesic recurrence plot (aka TCM)
if bsinglemember
    D_geo = distances(g,'Method','unweighted');
else
    D_geo = TCMdistance(g,nodemembers);
end
subplot(1,2,2)
imagesc(t,t,D_geo);
cb_ = colorbar;
cb_.Label.String = 'path length';
colormap(gca, 'hot')
axis square
xlabel('time (s)')
ylabel('time (s)')
title('geodesic recurrence plot')
h2 = gca;
end

