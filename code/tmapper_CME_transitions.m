% tmapper_CME_transitions.m
%{

plot geodesics averaged across subjects 
%}

addpath(genpath('../'))
addpath(genpath('../../GWTLBgpu/'))
addpath('plot_utils','tmapper_tools','other_utils')
%% ===== configure paths, import data ===== %%
clear all
% close all
clc

k=5;% default 5
d=2;% deftaul 2

savepath = '../results/CME/'; % path for behavioral data and shape graphs
figpath = [savepath 'cycleDistr/'];
% load([ savepath 'sumdat_tmapper_' par2filename({'k',k,'d',d}) '.mat'])
load([ savepath 'sumdat_tmapper_cyc_' par2filename({'k',k,'d',d}) '.mat'])
N_ss = height(sumdat); % number of subjects
N_tasks = 5;
N_t = length(x_label);
%% ===== compute source/sink distance for each subject ===== %%
[sourceAvgDis,sinkAvgDis,sourcesinkDiff]=deal(nan(N_ss, N_t));

for n=1:N_ss
    sourceAvgDis(n,:) = mean(normtcm(TCMdistance(sumdat.g_simp{n},sumdat.members{n})),2); % average distance from a node (source) to other nodes
    sinkAvgDis(n,:) = mean(normtcm(TCMdistance(sumdat.g_simp{n},sumdat.members{n})),1)'; % average distance to a node (sink) from other nodes  
    sourcesinkDiff(n,:) = (sourceAvgDis(n,:)-sinkAvgDis(n,:));
%     sourcesinkDiff(n,:) = abs(sourceAvgDis(n,:)-sinkAvgDis(n,:));
end
%% ===== plot source-sink difference for a single subject ===== %%
ssIdx = 10;%15;
node_sourcesinkDiff=cellfun(@(x) nanmean(sourcesinkDiff(ssIdx,x)), sumdat.members{ssIdx,1});

% -- node colormap
cmap = vals2colormap(node_sourcesinkDiff,'jet');
% -- node size
nodesize = cell2mat(cellfun(@(x) length(x), sumdat.members{ssIdx,1}, 'UniformOutput',0));
nodesize = rankval(nodesize);% the marker size reflects the rank of the node size
nodesize = rescale(nodesize, 1, 10);
% -- plotting
figure('position',[3418 177 728 571]);
hg=plot(sumdat.g_simp{ssIdx},'EdgeAlpha',0.3,'EdgeColor','k','NodeCData',node_sourcesinkDiff,'NodeLabel','',...
    'Layout','force','ArrowSize',5,'MarkerSize',nodesize); 
axis equal
axis off
cb=colorbar;
cb.Label.String = "source-sink difference (a.u.)";
% colormap(gca, par.cmap)
caxis([-1 1]*max(abs(node_sourcesinkDiff)))
colormap bluewhitered
h1 = gca;

% -- additional visual details
h.ArrowPosition=0.7;
h.ArrowSize=6;
h.EdgeAlpha=0.5;
title(sumdat.id{ssIdx})
%% ===== calculating mean and se across subjects ===== %%
[srdisMean,srdisSE]=grpstats(sourceAvgDis,ones(N_ss,1),{'mean','sem'});
[skdisMean,skdisSE]=grpstats(sinkAvgDis,ones(N_ss,1),{'mean','sem'});
[ssdiffMean,ssdiffSE]=grpstats(sourcesinkDiff,ones(N_ss,1),{'mean','sem'});

%% ===== plot ====== %%
xl = [min(t),max(t)];
yl1 = [0.1 0.5];
yl2 = [-0.26 0.26];
% yl2 = [0 0.35];
figure('position',[10,10,1000,500]); 
% - average distance as a source vs a sink
subplot(2,1,1)
hold on
xlim(xl)
ylim(yl1)
addBackBlock(gca,x_label,CMEcmap,'xlim',xl,'ylim',yl1,...
    'ignore',1,'ulabelnames',CMEtasknames);
[hm,hse]=plotSpectraSEM(repmat(t',2,1),[srdisMean;skdisMean],[srdisSE;skdisSE],'cmap',[0,0,1;1,0,0]);

legend(hm,'as source', 'as sink');
ylabel({'normalized','distance','(a.u.)'});

% - difference in average distance as a source vs a sink
subplot(2,1,2)
hold on
xlim(xl)
ylim(yl2)
addBackBlock(gca,x_label,CMEcmap,'xlim',xl,'ylim',yl2,...
    'ignore',1,'ulabelnames',CMEtasknames);
plot(xl,[0 0],':','Color',0.3*ones(1,3) )
[hm_,hse_]=plotSpectraSEM(t',ssdiffMean,ssdiffSE,'cmap',[0,0,0]);
legend(hm_,'source-sink');
xlabel('time (s)');
ylabel({'normalized','distance','(a.u.)'});

%% ===== extract on the source-sink separation in different tasks ===== %%
win = 20; % (TR) window of initial and final segment of each task block 

% -- store statistics for each tast
srsksep = nan(N_ss, N_tasks);
srskavg = nan(N_ss, N_tasks);
srskavg_mid = nan(N_ss, N_tasks);

for taskIdx = 2:N_tasks
    [block_start, block_end, block_size]=findtaskn(x_label==taskIdx);
    
    % -- average diff in the initial segment of the blocks
    seginit= mean([sourcesinkDiff(:,block_start(1):block_start(1)+win-1),...
                    sourcesinkDiff(:,block_start(2):block_start(2)+win-1)],2);
                
    % -- average diff in the final segment of the blocks
    segend= mean([sourcesinkDiff(:,block_end(1)-win+1:block_end(1)),...
                    sourcesinkDiff(:,block_end(2)-win+1:block_end(2))],2);
    
    % -- average diff in the middle segment of the blocks
    srskavg_mid(:,taskIdx) = mean([sourcesinkDiff(:,block_start(1)+win:block_end(1)-win),...
                        sourcesinkDiff(:,block_start(2)+win:block_end(2)-win)],2);
        
    % -- diff between initial and final seg
    srsksep(:,taskIdx) = seginit - segend;
    srskavg(:,taskIdx) = (abs(seginit)+abs(segend))/2;
end
%% ===== extract on the abs source-sink diff in different tasks ===== %%
win = 20; % (TR) window of initial and final segment of each task block 

% -- store statistics for each tast
srsksep = nan(N_ss, N_tasks);
srskavg = nan(N_ss, N_tasks);
srskavg_mid = nan(N_ss, N_tasks);

for taskIdx = 2:N_tasks
    [block_start, block_end, block_size]=findtaskn(x_label==taskIdx);
    
    % -- average diff in the initial segment of the blocks
    seginit= mean(abs([sourcesinkDiff(:,block_start(1):block_start(1)+win-1),...
                    sourcesinkDiff(:,block_start(2):block_start(2)+win-1)]),2);
                
    % -- average diff in the final segment of the blocks
    segend= mean(abs([sourcesinkDiff(:,block_end(1)-win+1:block_end(1)),...
                    sourcesinkDiff(:,block_end(2)-win+1:block_end(2))]),2);
    
    % -- average diff in the middle segment of the blocks
    srskavg_mid(:,taskIdx) = mean(abs([sourcesinkDiff(:,block_start(1)+win:block_end(1)-win),...
                        sourcesinkDiff(:,block_start(2)+win:block_end(2)-win)]),2);
        
    % -- diff between initial and final seg
    srsksep(:,taskIdx) = seginit - segend;
    srskavg(:,taskIdx) = (seginit+segend)/2;
end
%% ===== compare source-sink separation across tasks ===== %%
edges=-0.5:0.05:0.5;
ctrs=edge2ctr(edges);
Nbins= length(ctrs);

srsksepdist = nan(Nbins,N_tasks);

for n = 1:N_tasks
    srsksepdist(:,n)=histcounts(srsksep(:,n),edges,'Normalization','probability');
end

figure
colorlines(plot(ctrs',srsksepdist), CMEcmap)
legend(CMEtasknames)

% -- stats
Y = srsksep(:,2:N_tasks);
grp = repmat(2:N_tasks,N_ss,1);
% grp(grp==5)=3;
% grp(grp==4)=2;
[p,tbl,stats]=anovan(Y(:),grp(:),'varname','taskidx');%,'nested',[0 1; 0 0]);
[c,m,h,gnames]=multcompare(stats);

[grpmean,grpse]=grpstats(Y(:),grp(:),{'mean','sem'});
% -- plotting
figure;
hb=bar(grpmean,0.5,'facecolor','w');
set(gca,'XTickLabel',CMEtasknames(2:N_tasks))
hold on
errorbar(1:4, grpmean,grpse,grpse,'LineStyle','none','color',[0 0 0])
ylabel('source-sink separation')
scatter(grp(:)-1+(repmat((1:18)',4,1)-9)*0.01,Y(:),'filled','MarkerFaceColor',lines(1))

[h,p,ci,stats]=ttest(Y(grp==5))
%% ===== correlation with performance ===== %%
taskIdx=2;[2,4];[3,5];
behscore=sumdat.AverageRT;
figure
scatter(mean(srsksep(:,taskIdx),2),behscore,'filled')
lsline
[r,p]=corr(mean(srsksep(:,taskIdx),2),behscore);
titlestr = par2title({'r',r,'p',p});
title(CMEtasknames(taskIdx) + " " + titlestr);

%% ===== compare overall source/sinkness at entrance and exit across tasks ===== %%
edges=0:0.05:1;
ctrs=edge2ctr(edges);
Nbins= length(ctrs);

srskavgdist = nan(Nbins,N_tasks);

for n = 1:N_tasks
    srskavgdist(:,n)=histcounts(srskavg(:,n),edges,'Normalization','probability');
end

figure
colorlines(plot(ctrs',srskavgdist), CMEcmap)
legend(CMEtasknames)

% -- stats
Y = (srskavg(:,2:N_tasks));
grp = repmat(2:N_tasks,N_ss,1);
% grp(grp==5)=3;
% grp(grp==4)=2;
[p,tbl,stats]=anovan(Y(:),grp(:),'varname','taskidx');%,'nested',[0 1; 0 0]);
[c,m,h,gnames]=multcompare(stats);

[grpmean,grpse]=grpstats(Y(:),grp(:),{'mean','sem'});
% -- plotting
figure;
hb=bar(grpmean,0.5,'facecolor','w');
set(gca,'XTickLabel',CMEtasknames(2:N_tasks))
hold on
errorbar(1:4, grpmean,grpse,grpse,'LineStyle','none','color',[0 0 0])
ylabel('absolute source-sink difference')
scatter(grp(:)-1+(repmat((1:18)',4,1)-9)*0.01,Y(:),'filled','MarkerFaceColor',lines(1))
ylim([0 0.3]) % -- for begin and

[h,p,ci,stats]=ttest(Y(grp==2))
%% ===== compare overall source/sinkness at different seg tasks ===== %%
% -- stats
Y = [srskavg(:,2:N_tasks); (srskavg_mid(:,2:N_tasks))];
grp1 = repmat(2:N_tasks,N_ss,1);% task group
grp1 = [grp1; grp1];
grp2 = [ones(N_ss,N_tasks-1);2*ones(N_ss,N_tasks-1)]; % segment group
% grp(grp==5)=3;
% grp(grp==4)=2;
[p,tbl,stats]=anovan(Y(:),[grp1(:), grp2(:)],...
    'varname',{'taskidx','segidx'},'model','full');%,'nested',[0 1; 0 0]);
[c,m,h,gnames]=multcompare(stats,'dimension',[1 2]);
[c,m,h,gnames]=multcompare(stats,'dimension',[1]);
[c,m,h,gnames]=multcompare(stats,'dimension',[2]);


[grpmean,grpse]=grpstats(Y(:),[grp1(:)],{'mean','sem'});
% -- plotting
figure;
hb=bar(grpmean,0.5,'facecolor','w');
set(gca,'XTickLabel',CMEtasknames(2:N_tasks))
% set(gca,'XTickLabel',{'entry/exit','middle'})
hold on
errorbar(1:length(unique(grp1(:))), grpmean,grpse,grpse,'LineStyle','none','color',[0 0 0])
ylabel('absolute source-sink difference')
scatter(grp1(:)-1+(repmat((1:18)',8,1)-9)*0.01,Y(:),'filled','MarkerFaceColor',lines(1))
% scatter(grp2(:)+(repmat((1:18)',8,1)-9)*0.01,Y(:),'filled','MarkerFaceColor',lines(1))

ylim([0 0.4]) 
