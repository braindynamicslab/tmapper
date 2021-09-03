% tmapper_CME_cycles.m
% --------------------
%{
Here we compare the tmapper shape graphs using measures that are
dynamically relevant. In particular, we study the cycles and decision
points in the graph.

%}

addpath(genpath('../'))
addpath(genpath('../../GWTLBgpu/'))
addpath('plot_utils','tmapper_tools','other_utils')
%% ===== configure paths, import data ===== %%
clear all
close all
clc

k=5;% default 5
d=2;% deftaul 2

savepath = '../results/CME/'; % path for behavioral data and shape graphs
figpath = [savepath 'cycleDistr/'];
% load([ savepath 'sumdat_tmapper_' par2filename({'k',k,'d',d}) '.mat'])
load([ savepath 'sumdat_tmapper_cyc_' par2filename({'k',k,'d',d}) '.mat'])
N_ss = height(sumdat); % number of subjects
%% ===== path decomposition ====== %%
ss_id=1;
dg = sumdat.g_simp{ss_id};
[~,~,~,allcycles]=CycleCount2p(dg.adjacency);
T=CycleCluster(allcycles,0.5);
addDiagBlock(gca,sort(T),zeros(max(T),3));

tic
[cluster_conn, cluster_conn_dir, ...
    clusters_nodes, clusters_boundary, clusters_interior,...
    clusters_crtpts, clusters_intcrtpts] ...
    = CycleClusterConn(dg, allcycles, T);
toc

allbd = unique(cell2mat(clusters_boundary));
tic
allupath=Cycles2Paths(allcycles,allbd);
toc
Tp=CycleCluster(allupath,0.5);
addDiagBlock(gca,sort(Tp),zeros(max(Tp),3));

[pcluster_conn, pcluster_conn_dir, ...
    pclusters_nodes, pclusters_boundary, pclusters_interior,...
    pclusters_crtpts, pclusters_intcrtpts] ...
    = CycleClusterConn(dg, allupath, Tp);

[allbd,pclusters_nodes,pclusters_interior,pclusters_boundary,...
    pcluster_conn_dir,pclusters_intcrtpts, allupath, Tp] ...
    = CyclePathDecomp(dg);


% -- plot transition map with labelled paths
figure('position',[810 10 800 800])
ho=plot(dg,'layout','force','NodeCData',zeros(dg.numnodes,1),...
    'MarkerSize',rescale(tiedrank(cellfun(@length,sumdat.members{ss_id})),5,25),'EdgeColor','k','NodeLabel',1:height(dg.Nodes),'linewidth',1);
colormap jet
hbp=colorbar;
hbp.Label.String = 'path cluster index';
axis square 
axis off
% title({'theoretically predicted shape graph', sprintf('for w_{EE}=%g, w_{EI}=%g',...
%     a.w_EE,a.w_EI)})

for n=max(Tp)+1-unique(Tp')
    ho.NodeCData(pclusters_nodes{n}) = n;
end
ho.NodeCData(unique(cell2mat(clusters_boundary))) = 0;
%% ===== distribution of cycles ===== %%
% for n = 1:N_ss
%     [sumdat.cycCount{n}, sumdat.cycLen{n}, sumdat.cycPath{n}]=...
%         CycleCount2p(sumdat.g_simp{n}.adjacency);
% end

% save([ savepath 'sumdat_tmapper_cyc_' par2filename({'k',k,'d',d}) '.mat'],...
%     'sumdat','t','tidx','wee','wei','x_label')
%% ===== visualize distribution of cycles ==== %%
behscore = 'MemPC';
cmap = jet(N_ss);
[~,idx]=sort(sumdat{:,behscore});
figure('position',[10 20 1000 500])
hold on
for n = 1:N_ss
    plot(sumdat.cycLen{idx(n)},sumdat.cycCount{idx(n)},'color',cmap(n,:))
end
legend(sumdat.id(idx))
xlabel('cycle length')
ylabel('count')
title(['colored by ' behscore])
xlim([1 50])

% print([ figpath 'cycLenbySS_nonorm_' behscore '_' ...
%     par2filename({'k',k,'d',d}) '.png'],'-dpng')
%% ===== visualize distribution of cycles (percentage) ==== %%
% -- normalize distributions: convert to percentage
sumdat.cycPrct = cellfun(@(y) y./sum(y),sumdat.cycCount,...
    'UniformOutput',0);
% -- plot ss-level distribution
% behscore = 'AveragePC';
cmap = jet(N_ss);
[~,idx]=sort(sumdat{:,behscore});
figure('position',[10 20 1000 500])
hold on
for n = 1:N_ss
    plot(sumdat.cycLen{idx(n)},100*sumdat.cycPrct{idx(n)},'color',cmap(n,:))
end
legend(sumdat.id(idx))
xlabel('cycle length')
ylabel('Percentage (%)')
title(['colored by ' behscore])
xlim([1 50])

% print([ figpath 'cycLenbySS_prct_' behscore '_' ...
%     par2filename({'k',k,'d',d}) '.png'],'-dpng')
%% ===== visualize distribution of cycles (max length) ==== %%
% -- normalize distributions: convert to percentage
sumdat.cycLenNorm = cellfun(@(y) y./max(y),sumdat.cycLen,...
    'UniformOutput',0);
% -- plot ss-level distribution
% behscore = 'AveragePC';
cmap = jet(N_ss);
[~,idx]=sort(sumdat{:,behscore});
figure('position',[10 20 1000 500])
hold on
for n = 1:N_ss
    plot(sumdat.cycLenNorm{idx(n)},100*sumdat.cycPrct{idx(n)},'color',cmap(n,:))
end
legend(sumdat.id(idx))
xlabel('relative cycle length')
ylabel('Percentage (%)')
title(['colored by ' behscore])
xlim([0 1])

% print([ figpath 'cycLenNormBySS_prct_' behscore '_' ...
%     par2filename({'k',k,'d',d}) '.png'],'-dpng')
%% ===== average cycle length ===== %%
mapperscore = 'STDCycLenNorm';

sumdat.AvgCycLen = cellfun(@(c,l) (c'*l)/sum(c),sumdat.cycCount, sumdat.cycLen);
sumdat.AvgCycLenNorm = cellfun(@(c,l) (c'*l)/sum(c),sumdat.cycCount, sumdat.cycLenNorm);
sumdat.STDCycLen = cellfun(@(p) ...
    std(cell2mat(cellfun(@(x) size(x,2),p,'uniformoutput',0))),...
    sumdat.cycPath);
sumdat.STDCycLenNorm = cellfun(@(p,l) ...
    std(cell2mat(cellfun(@(x) size(x,2),p,'uniformoutput',0))/max(l)),...
    sumdat.cycPath,sumdat.cycLen);

figure;
scatter(sumdat{:,mapperscore},sumdat{:,behscore})
xlabel(mapperscore)
ylabel(behscore)
[RHO,P]=corr(sumdat{:,mapperscore},sumdat{:,behscore},'type','Pearson');
disp(['rho=' num2str(RHO) ', p=' num2str(P)])
%% ===== cycle distribution by task performance ==== %%
behscore = 'AveragePC';
normMethod = 'count';% 'count' or 'pdf'

% -- median split by task performance
[~,idx]=sort(sumdat{:,behscore});
grpvar = nan(N_ss,1);
grpvar(idx(1:N_ss/2)) = 1;% lower group
grpvar(idx(N_ss/2+1:end)) = 2;% higher group

% -- compute a histogram
binctrs = 2:2:50;
binedges = 1:2:51;
N_bins = length(binctrs);
sumdat.cycHist=cellfun(@(l,c) histcounts(repelem(l,c,1),binedges,'Normalization',normMethod),...
    sumdat.cycLen, sumdat.cycCount,'UniformOutput',0);
[cycHsMean,cycHsSE]=grpstats(cell2mat(sumdat.cycHist),grpvar,{'mean','sem'});
% -- average and plot
figure
plotSpectraSEM(binctrs,cycHsMean,cycHsSE,'cmap',[0,0,1;1,0,0]);
xlabel('cycle length (TR)')
ylabel(normMethod)
title(behscore)
legend('mean(-)','SE(-)','mean(+)','SE(+)')
xlim([min(binctrs) max(binctrs)])

% print([ figpath 'cycLen_' normMethod '_splitBy_' behscore '_' ...
%     par2filename({'k',k,'d',d}) '.png'],'-dpng')

% -- statistical tests
Y=cell2mat(sumdat.cycHist);
% :: ANOVA
[grp_bin, grp_beh]=meshgrid(1:N_bins,grpvar);
figure
[p,tbl,stats,terms]=anovan(Y(:),[grp_beh(:), grp_bin(:)],...
    'model','full','varnames',{behscore,'bin index'});%,'nested',[0 1; 0 0]);
[c,m,h,gnames]=multcompare(stats,'Dimension',[1 2]);

% -- ttests
[h,p,ci,stats]=ttest2(Y(grpvar==1,:),Y(grpvar==2,:));
%% ===== task distribution on cycles ===== %%
% -- find the task labels for each TR on cycles
findCycleTask = @(x,members) histcounts(x_label(cell2mat(members(unique(x(:))))),0.5:5.5); % find task label for all nodes on cycles of a particular length.
sumdat.cycleTaskCount = cellfun(@(X, members) cell2mat(cellfun(@(x) findCycleTask(x,members), X, 'uniformoutput',0)),...
        sumdat.cycPath, sumdat.members, 'uniformoutput',0);
% -- fraction of different tasks for each given cycle length
sumdat.cycleTaskPrc = cellfun(@(X) X./sum(X,2), ...
                                sumdat.cycleTaskCount,'UniformOutput',0);
% -- fraction of different cycle length for each given task
sumdat.cyclePrcByTask = cellfun(@(X) X./sum(X,1), ...
                                sumdat.cycleTaskCount,'UniformOutput',0);
%% -- plot the two fractions
ssIdx = 18;

figure('position',[50 50 1120 420]);
xl = [min(sumdat.cycLen{ssIdx})-1 max(sumdat.cycLen{ssIdx})+1];
subplot(1,2,1)
hp = plot(sumdat.cycLen{ssIdx},sumdat.cycleTaskPrc{ssIdx});
colorlines(hp,CMEcmap)
xlim(xl)
ylim([0 1])
xlabel('cycle length (TR)')
ylabel('fraction of task')
title(sumdat.id{ssIdx})
legend(CMEtasknames,'Color','none','EdgeColor','None')

subplot(1,2,2)
hp = plot(sumdat.cycLen{ssIdx},sumdat.cyclePrcByTask{ssIdx});
colorlines(hp,CMEcmap)
xlim(xl)
xlabel('cycle length (TR)')
ylabel('fraction of cycles')
title(sumdat.id{ssIdx})

print([figpath 'ss_level/subj' sumdat.id{ssIdx}(4:5) '_' par2filename({'k',k,'d',d}) '_distrTaskCyc.png'],'-dpng')