function [clusterIdx] = CycleCluster(allcycles,thres,varargin)
%CYCLECLUSTER cluster cycles based on the overlap between them. The
%fraction of overlap between two cycles is the number of shared nodes
%divided by the number of the union of nodes from two cycles. 
%   [clusterIdx] = CycleCluster(allcycles,thres,...)
% input:
%   allcycles: N-by-1 cell array, each cell contains the path of one cycle
%   thres: a number between 0 and 1, a cut-off threshold for single-linkage
%   clustering. Such that if overlap>thres, the two cycles belong to
%   the same cluster. 
% output: 
%   clusterIdx: N-by-1 integer array, each integer is the index of a
%   loop-cluster. Integers are from 1 to M. 
% parameters:
%   plotmat: plot the overlap matrix (N-by-N). Default true.
%   plotmds: plot the 2D mds projection. Default false.
%   plothist: plot the histogram of linkage distance. Default true.
%   reordermat: reorder rows/cols of overlap matrix according to cluster
%   assignment. Default true.
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 9-1-2020 ~
%}

p=inputParser;
p.addParameter('plotmat',true)
p.addParameter('plotmds',false)
p.addParameter('plothist',true)
p.addParameter('reordermat',true)
p.parse(varargin{:})
par=p.Results;


% -- compute overlap between cycles
Nc=length(allcycles);
prct_overlap = nan(Nc,Nc);

for ii=1:Nc
    for jj=1:Nc 
        prct_overlap(ii,jj) = length(intersect(allcycles{ii},allcycles{jj}))...
            /length(union(allcycles{ii},allcycles{jj}));
    end
end


% -- clustering
% :: visualize clusters with MDS
if par.plotmds
    Y=cmdscale(1-prct_overlap,2);
    figure
    scatter(Y(:,1),Y(:,2))
end

% :: single linkage
Z=linkage(squareform(1-prct_overlap),'complete');
% :: visualize linkage
if par.plothist
    % figure
    % dendrogram(Z)
    figure;
    histogram(1-Z(:,3),linspace(0,1,21))
    hold on
    plot(thres*[1 1],ylim)
    xlabel('overlap')
    ylabel('count')
end

clusterIdx = cluster(Z,'cutoff',1-thres,'criterion','distance');

% -- visualize overlaps between cycles
if par.plotmat
    figure
    if par.reordermat
        [~,cycorder]=sort(clusterIdx);
        imagesc(prct_overlap(cycorder,cycorder));
    else
        imagesc(prct_overlap);
    end
    set(gca,'ydir','normal')
    hcb=colorbar;
    hcb.Label.String='fraction of overlap';
    caxis([0 1])
    axis square
    if par.reordermat
        xlabel('cycle/path index (reordered)')
        ylabel('cycle/path index (reordered)')
    else
        xlabel('cycle/path index')
        ylabel('cycle/path index')
    end
end

end

