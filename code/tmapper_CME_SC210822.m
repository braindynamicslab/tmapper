% tmapper_CME.m
% --------------
%{
analyses of temporal mapper construction based on the CME data. 
(created by MZ, 2-11-2020)
%}

addpath(genpath('../'))
addpath(genpath('../../GWTLBgpu/'))
addpath('plot_utils','tmapper_tools','other_utils')

%% ===== load data ===== %%
clear all
close all
clc


ss_id = 1; % subject id
load(sprintf('../data/CME/SBJ%02d_Shine_375.mat',ss_id))
savepath = '../results/CME/ss_level/';
%% ===== construct tmapper ===== %%
k=7;% default 5
d=3;% deftaul 2
fprefix = ['subj' num2str(ss_id,'%02d') '_k' num2str(k) '_d' num2str(d)];

tic
g = tknndigraph (x,k,tidx,'reciprocal',true,'timeExcludeSpace', true); % tknn
[g_simp, members, ~] = filtergraph(g,d,'reciprocal',true);% simplified graph
toc

% -- plotting
[a1, a2, cb,~,h] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,15],'cmap',CMEcmap);
cb=CMEcbar(cb);
% [a1, a2] = plotgraphtcm(g,x_label,t,[],'nodesizerange',[1,10]);
title(a1,['subj-' num2str(ss_id,'%02d') ', k=' num2str(k) ', d=' num2str(d)]);
title(a2,'geodesic recurrence plot')

% -- additional annotations
colormap(a2,'gray')
% colormap(a2,'hot')
a2=addDiagBlock(a2,x_label,CMEcmap);

% -- additional visual details
set(gcf,'position', [700, 490, 1332, 491])
h.ArrowPosition=0.7;
h.ArrowSize=6;
h.EdgeAlpha=0.5;
% -- save
% print([savepath fprefix '_tgraph.png'],'-dpng')
%% ===== inter- and intra-block distance ===== %%
% --- compute TCM (normalized)
tcm = normtcm(TCMdistance(g_simp,members));
% --- plot block statistics
plotBlockMatStats(tcm,x_label,'x',t,'y',t,'xlabel','time (s)', ...
    'ylabel','time (s)', 'blockcmap',CMEcmap, 'ptnames', CMEtasknames([],'short',1),...
    'asymOpt','afterStats');

% -- save
print([savepath fprefix '_blockAvg.png'],'-dpng')
%% ===== inter- and intra-task distance ===== %%
% --- compute TCM (normalized)
tcm = normtcm(TCMdistance(g_simp,members));
% --- plot task statistics
[sortedlab, sortedlabIdx] = sort(x_label);
plotBlockMatStats(tcm(sortedlabIdx,sortedlabIdx),sortedlab,'x',t,'y',t,'xlabel','', ...
    'ylabel','', 'blockcmap',CMEcmap, 'ptnames', CMEtasknames([],'short',1),...
    'asymOpt','afterStats','tickOpt','block');

% -- save
print([savepath fprefix '_taskAvg.png'],'-dpng')
return
%% ===== plot distance distribution ===== %%
% --- compute TCM
tcm = normtcm(TCMdistance(g_simp,members));
N_t = length(tcm);
tcm_row = mat2cell(tcm,ones(1,N_t),N_t);
tcm_col = mat2cell(tcm',ones(1,N_t),N_t);

% --- compute histograms
nbins = 20;
binedges = linspace(0,1,nbins+1);
histfun = @(x) histcounts(x,binedges,'Normalization','pdf');
source_hist = cell2mat(cellfun(histfun,tcm_row,'UniformOutput',0))';
sink_hist = cell2mat(cellfun(histfun,tcm_col,'UniformOutput',0))';

% --- plot distance distribution as a function of time
xl = [min(t),max(t)];
yl_annot = [-0.05 1.05];
cbar_xpos = 0.92;
figure('position',[10,10,1000,500]); 
% - source distribution
subplot(2,1,1)
hold on
ylim(yl_annot)
xlim(xl)
imagesc(xl,[0 1],source_hist)
colormap gray
set(gca,'ydir','normal')
addBackBlock(gca,x_label,CMEcmap,'xlim',xl,'ylim',yl_annot,'ignore',1,'ulabelnames',CMEtasknames,'alpha',0.4);
ylabel({'distance','as source','(a.u.)'})
cb = colorbar;
cb.Position(1) = cbar_xpos;
xlabel(cb,'probability density')

% - sink distribution
subplot(2,1,2)
hold on
ylim(yl_annot)
xlim(xl)
imagesc(xl,[0 1],sink_hist)
colormap gray
set(gca,'ydir','normal')
addBackBlock(gca,x_label,CMEcmap,'xlim',xl,'ylim',yl_annot,'ignore',1,'ulabelnames',CMEtasknames,'alpha',0.4);
ylabel({'distance','as sink','(a.u.)'})
xlabel('time (s)')
cb = colorbar;
cb.Position(1) = cbar_xpos;
xlabel(cb,'probability density')

% -- save
% print([savepath fprefix '_histdistdyn.png'],'-dpng')
%% ===== plot average distance ===== %%
sourceAvgDis = mean(normtcm(TCMdistance(g_simp,members)),2); % average distance from a node (source) to other nodes
sinkAvgDis = mean(normtcm(TCMdistance(g_simp,members)),1)'; % average distance to a node (sink) from other nodes 

xl = [min(t),max(t)];

figure('position',[10,10,1000,500]); 
% - average distance as a source vs a sink
subplot(2,1,1)
hold on
xlim(xl)
ylim([0 1])
addBackBlock(gca,x_label,CMEcmap,'xlim',xl,'ylim',[0,1],...
    'ignore',1,'ulabelnames',CMEtasknames);
plot(t,sourceAvgDis);
plot(t,sinkAvgDis);
legend('as source', 'as sink');
ylabel({'normalized','distance','(a.u.)'});

% - difference in average distance as a source vs a sink
subplot(2,1,2)
hold on
xlim(xl)
ylim([-1 1])
addBackBlock(gca,x_label,CMEcmap,'xlim',xl,'ylim',[-1,1],...
    'ignore',1,'ulabelnames',CMEtasknames);
plot(xl,[0 0],'w')
lh=plot(t,sourceAvgDis-sinkAvgDis,'k');
legend(lh,'source-sink');
xlabel('time (s)');
ylabel({'normalized','distance','(a.u.)'});

% -- save
% print([savepath fprefix '_avgdistdyn.png'],'-dpng')
%% ===== null model for parameter selection ===== %%
