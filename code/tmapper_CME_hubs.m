% tmapper_CME_hubs.m
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
N_tasks = 5;
%% ===== degrees of each node of each subject ===== %%
for n=1:N_ss
    sumdat.indeg{n}=full(sum(sumdat.g_simp{n}.adjacency,1)');
    sumdat.outdeg{n}=full(sum(sumdat.g_simp{n}.adjacency,2));
    sumdat.deg{n}= sumdat.indeg{n} + sumdat.outdeg{n};
end
%% ===== split nodes to high and low deg ===== %%
% :: split by rank
% rank_thres = 5; % top N as high-degree nodes
% highdegidx=cellfun(@(x) max(tiedrank(x))-tiedrank(x)<rank_thres,sumdat.deg,'UniformOutput',0);
% :: split by percentile
pct_thres = 2; % top x% to nodes
highdegidx=cellfun(@(x) x>=prctile(x,100-pct_thres),sumdat.deg,'UniformOutput',0);

N_highdeg=cellfun(@sum,(highdegidx));
highdegAvgTR=cellfun(@(idx,mem) length(unique(vertcat(mem{idx}))),...
    highdegidx,sumdat.members)./N_highdeg;
highdegAvgTR_bytask=nan(N_ss,N_tasks);
for n=1:N_tasks
    highdegAvgTR_bytask(:,n)=cellfun(@(idx,mem) sum(x_label(unique(vertcat(mem{idx})))==n),...
        highdegidx,sumdat.members)./N_highdeg;
end

rankAvgPC = tiedrank(sumdat.AveragePC);
[rankAvgPC,idx]=sort(rankAvgPC);
figure;
colorlines(plot(rankAvgPC,highdegAvgTR_bytask(idx,:)),CMEcmap)
legend(CMEtasknames,'location','best')
ylabel('#TR / node')
xlabel('Rank of Average PC')
%% ===== distr of #TR node by degree rank thres ===== %%
% -- define task
taskIdx=5;
% -- threshold
N_thres = 50;
rank_thres_set = logspace(log10(5),log10(400),N_thres);

% -- histogram parameters
maxNTR = 50;
N_bins = maxNTR;
edges=linspace(0.5,maxNTR+0.5,N_bins+1);
ctrs=edge2ctr(edges);
% -- preallocate memory to store distribution
distrXthres=nan(N_bins,N_thres);
for n = 1:N_thres
    % :: find the index of high-degree nodes
    highdegidx=cellfun(@(x) max(tiedrank(x))-tiedrank(x)<rank_thres_set(n),...
        sumdat.deg,'UniformOutput',0);
    % :: count task-TRs in nodes
    N_highdeg=cellfun(@sum,(highdegidx));
    % :: average # TR/node
    highdegAvgTR=cellfun(@(idx,mem) sum(x_label(unique(vertcat(mem{idx})))==taskIdx),...
        highdegidx,sumdat.members)./N_highdeg;
    % :: calculate distributions of #TR/node
    distrXthres(:,n) = histcounts(highdegAvgTR,edges,'Normalization','pdf');
end

% -- plot
figure
imagesc(rank_thres_set, ctrs, distrXthres)
set(gca,'XScale','log','ydir','normal')
xlabel('#nodes included')
ylabel('#TR / node')
title(CMEtasknames(taskIdx))
cb=colorbar;
caxis([0 1])
xlim([min(rank_thres_set) max(rank_thres_set)])
cb.Label.String = 'probability density';
%% ===== distr of TR fraction by degree rank thres ===== %%
% -- define task
taskIdx=3;
% -- threshold
N_thres = 50;
rank_thres_set = logspace(log10(5),log10(400),N_thres);

% -- histogram parameters
N_bins = 50;
edges=linspace(0,1,N_bins+1);
ctrs=edge2ctr(edges);
% -- preallocate memory to store distribution
distrXthres=nan(N_bins,N_thres);
for n = 1:N_thres
    % :: find the index of high-degree nodes
    highdegidx=cellfun(@(x) max(tiedrank(x))-tiedrank(x)<rank_thres_set(n),...
        sumdat.deg,'UniformOutput',0);
    % :: fraction of TRs in high-degree nodes
    highdegAvgTR=cellfun(...
        @(idx,mem) sum(x_label(unique(vertcat(mem{idx})))==taskIdx)...
        ./length(unique(vertcat(mem{idx}))),...
        highdegidx,sumdat.members);
    % :: calculate distributions of #TR/node
    distrXthres(:,n) = histcounts(highdegAvgTR,edges,'Normalization','pdf');
end

% -- plot
figure
imagesc(rank_thres_set, ctrs, distrXthres)
set(gca,'XScale','log','ydir','normal')
xlabel('#nodes included')
ylabel('fraction of TRs')
title(CMEtasknames(taskIdx))
cb=colorbar;
caxis([0 N_bins])
cb.Label.String = 'probability density';
xlim([min(rank_thres_set) max(rank_thres_set)])

%% ===== distr of TR fraction by degree percentile thres ===== %%
% -- define task
taskIdx=5;
% -- threshold
N_thres = 50;
prt_thres_set = logspace(log10(1.25),log10(100),N_thres);

% -- histogram parameters
N_bins = 50;
edges=linspace(0,1,N_bins+1);
ctrs=edge2ctr(edges);
% -- preallocate memory to store distribution
distrXthres=nan(N_bins,N_thres);
for n = 1:N_thres
    % :: find the index of high-degree nodes
    highdegidx=cellfun(@(x) x>=prctile(x,100-prt_thres_set(n)),...
        sumdat.deg,'UniformOutput',0);
    % :: fraction of TRs in high-degree nodes
    highdegAvgTR=cellfun(...
        @(idx,mem) sum(x_label(unique(vertcat(mem{idx})))==taskIdx)...
        ./length(unique(vertcat(mem{idx}))),...
        highdegidx,sumdat.members);
    % :: calculate distributions of #TR/node
    distrXthres(:,n) = histcounts(highdegAvgTR,edges,'Normalization','pdf');
end

% -- plot
figure
imagesc(prt_thres_set, ctrs, distrXthres)
set(gca,'XScale','log','ydir','normal')
xlabel('% nodes included')
ylabel('fraction of TRs')
title(CMEtasknames(taskIdx))
cb=colorbar;
caxis([0 N_bins])
cb.Label.String = 'probability density';
xlim([min(prt_thres_set) max(prt_thres_set)])
%% ===== average TR fraction by degree percentile thres ===== %%

% -- threshold
N_thres = 10;
prt_thres_set = [100./2.^[N_thres-1:-1:0]];%logspace(log10(1.25),log10(100),N_thres);

% -- preallocate memory to store distribution
sumTRXthres=nan(N_thres,N_tasks);

for taskIdx=1:N_tasks
    for n = 1:N_thres
        % :: find the index of high-degree nodes
        highdegidx=cellfun(@(x) x>=prctile(x,100-prt_thres_set(n)),...
            sumdat.deg,'UniformOutput',0);
        % :: fraction of TRs in high-degree nodes
        highdegAvgTR=cellfun(...
            @(idx,mem) sum(x_label(unique(vertcat(mem{idx})))==taskIdx)...
            ./length(unique(vertcat(mem{idx}))),...
            highdegidx,sumdat.members);
        % :: calculate distributions of #TR/node
        sumTRXthres(n,taskIdx) = nansum(highdegAvgTR);
    end
end

fracTRXthres = sumTRXthres ./ sum(sumTRXthres,2);

% -- plot group average
figure
b=bar(fracTRXthres,'stacked');
set(gca,'xtick',[1,5,10],'xticklabel',prt_thres_set([1,5,10]))
for n=1:N_tasks
    b(n).FaceColor=CMEcmap(n);
end
ylim([0 1])
legend(CMEtasknames,'location','northeastoutside')
xlabel('% nodes included')
ylabel('fraction of TRs')

%% ===== TR fraction for each subj and each task ===== %%
N_tasks = 5;
prt_thres = 2; % top 2 percent nodes (of high degree)
sumdat.fTR_highdeg = nan(N_ss,N_tasks);

% :: find the index of high-degree nodes
sumdat.highdegidx=cellfun(@(x) find(x>=prctile(x,100-prt_thres)),...
    sumdat.deg,'UniformOutput',0);
for n =1:N_tasks 
    % :: fraction of TRs in high-degree nodes
    sumdat.fTR_highdeg(:,n)=cellfun(...
        @(idx,mem) sum(x_label(unique(vertcat(mem{idx})))==n)...
        ./length(unique(vertcat(mem{idx}))),...
        sumdat.highdegidx,sumdat.members);
end
%% ===== correlation between TR fraction and task performance ===== %%
summemmat = sum(sumdat.fTR_highdeg(:,[3,5]), 2);
diffmemmat = diff(sumdat.fTR_highdeg(:,[3,5]),[], 2);
performvaridx = 11:13;
perform = sumdat{:,performvaridx};
[r,p]=corr(perform,summemmat);
YLIMS=[65, 100; 1 2.4; 0 45];
for n=1:size(perform,2)
    mdl=fitlm(summemmat,perform(:,n));
    % -- plot
    figure;
    scatter(summemmat,perform(:,n),'filled')
    lsline
    xlabel('TR fraction of memory + math')
    ylabel(sumdat.Properties.VariableNames{performvaridx(n)})
    axis square
    xlim([0.3 1])
    ylim(YLIMS(n,:))
    par2title({'r',r(n),'p',p(n), 'R^2',mdl.Rsquared.Ordinary})
end

% legend(sumdat.Properties.VariableNames{performvaridx})
%% ===== null models for correlation ===== %%
N_perm = 1000;
randseed=7;
rng(randseed)

[r_null, p_null, rsq_null] = deal(nan(N_perm,3));

% -- permutations
tic
for n = 1:N_perm
    disp(n)
    randmemmat = summemmat(randperm(N_ss));
    [r_null(n,:),p_null(n,:)]=corr(perform,randmemmat);
    for m=1:3
        mdl=fitlm(randmemmat,perform(:,m));
        rsq_null(n,m) = mdl.Rsquared.Ordinary;
    end
    toc
end

% -- plotting null distributions
N_bins = 20;
edges = linspace(-1,1,N_bins+1);

for n=1:3
    disp(confint(r_null(:,n)))
    figure
    histogram(r_null(:,n),edges,'Normalization','pdf','FaceColor','k','EdgeColor','w')
    hold on
    plot(repmat(confint(r_null(:,n)),2,1),repmat(ylim',1,2),':r')
    plot(ones(2,1)*r(n),ylim,'g')
    xlabel('r')
    ylabel('probability density')
    legend('null distribution','95% CI','','real')
    title(sumdat.Properties.VariableNames{performvaridx(n)})
end

figure
histogram(p_null)
figure
histogram(rsq_null)
%% ===== CCA ===== %%
% [A,B,r,U,V,stats]=canoncorr(sumdat{:,11:13},sumdat.fTR_highdeg(:,[3,5]));
% stats
% r
% figure
% scatter(U(:,1),V(:,1),'filled')
% hold on
% scatter(U(:,2),V(:,2),'filled')
%% ===== Recurrent processes through the hubs ====== %%
sumdat.allcycles = cellfun(@reorgCycles,sumdat.cycPath,'UniformOutput',0);
sumdat.allcyclelen = cellfun(@(x) cellfun(@length,x), sumdat.allcycles,'UniformOutput',0);

% -- find index of cycles that pass through the hubs
for n=1:N_ss
    sumdat.highdegCycIdx{n,1}=...
        find(cellfun(@(cyc) any(ismember(sumdat.highdegidx{n},cyc)), sumdat.allcycles{n}));
end
% -- find length of cycles
sumdat.highdegCycLen = cellfun(@(idx,cyclen) cyclen(idx), ...
    sumdat.highdegCycIdx, sumdat.allcyclelen,'UniformOutput',0);

% -- plot histograms of distributions over cycle length
N_bins = 21;
edges = linspace(1.5,85.5,N_bins+1);
ctrs = edge2ctr(edges);
% :: color by performance
behscore = 'MemRT';
cmap=rankcolor(jet(N_ss),sumdat{:,behscore});
figure
hold on
for n=1:N_ss
    distr_cyc=histcounts(sumdat.highdegCycLen{n},edges,'Normalization','pdf');
    plot(ctrs,distr_cyc,'Color',cmap(n,:))
end

% -- plot histrogram of distributions over cycle length percentage
% N_bins = 20;
% edges = linspace(0,100,N_bins+1);
% ctrs = edge2ctr(edges);
% figure
% hold on
% for n=1:N_ss
%     distr_cyc=histcounts(sumdat.highdegCycLen{n},prctile(sumdat.highdegCycLen{n},edges),...
%         'Normalization','cdf');
%     plot(ctrs,distr_cyc,'Color',cmap(n,:))
% end
%% ===== cycle distribution by task performance ==== %%
behscore = 'VidPM';
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
sumdat.cycHist=cellfun(@(x) histcounts(x,binedges,'Normalization',normMethod),...
    sumdat.highdegCycLen,'UniformOutput',0);% allcyclelen %highdegCycLen
[cycHsMean,cycHsSE]=grpstats(cell2mat(sumdat.cycHist),grpvar,{'mean','sem'});
% -- average and plot
figure
plotSpectraSEM(binctrs,cycHsMean,cycHsSE,'cmap',[0,0,1;1,0,0]);
xlabel('cycle length')
ylabel(normMethod)
title(behscore)
xlim([min(binctrs) max(binctrs)])
ylim([0 60])
ha=gca;
% print([ figpath 'cycLen_' normMethod '_splitBy_' behscore '_' ...
%     par2filename({'k',k,'d',d}) '.png'],'-dpng')

% -- statistical tests
Y=cell2mat(sumdat.cycHist);
% :: ANOVA
[grp_bin, grp_beh]=meshgrid(1:N_bins,grpvar);
figure
[p,tbl,stats,terms]=anovan(Y(:),[grp_beh(:), grp_bin(:)],...
    'model','full','varnames',{behscore,'bin index'});%,'nested',[0 1; 0 0]);
[c,m,h,gnames]=multcompare(stats,'Dimension',1,'Display','off');
disp(['p=' num2str(c(end))])
% [c,m,h,gnames]=multcompare(stats,'Dimension',2);
[grpmean,grpse]=grpstats(Y(:),grp_beh(:),{'mean','sem'})
[c,m,h,gnames]=multcompare(stats,'Dimension',[1 2]);
% :: index of within-bin comparisons
compidx=ismember(c(:,1:2), [(1:2:N_bins*2-1)',(2:2:N_bins*2)'],'row');
% :: plot significance levels
addStars(ha,binctrs, max(cycHsMean) + 12, c(compidx,end))
title(legend(ha,'SE(-)','mean(-)','SE(+)','mean(+)','location','best'),behscore)

% -- ttests
[h,p,ci,stats]=ttest2(Y(grpvar==1,:),Y(grpvar==2,:));
%% ===== cycle overlap by task performance ===== %%
binwidth = 5;
maxCycLen = 55;
edges = 1.5:binwidth:maxCycLen+0.5;
N_bins = length(edges)-1;
ctrs = edge2ctr(edges);

overlap_cyc = nan(N_ss,N_bins);

for nss=1:N_ss
    for nbin=1:N_bins
        disp([nss nbin])
        tic
        % -- extract cycles within a range of cycle length
        cycidx = (sumdat.highdegCycLen{nss}>edges(nbin)) ...
            & (sumdat.highdegCycLen{nss}<edges(nbin+1));
        % -- average overlap between cycles
        overlaptmp = CyclePathOverlap(sumdat.allcycles{nss}(sumdat.highdegCycIdx{nss}(cycidx))');
        overlap_cyc(nss,nbin) = nanmean(belowdiag(overlaptmp));
        toc
    end
end



% :: color by performance
behscore = 'AveragePM';
cmap=rankcolor(jet(N_ss),sumdat{:,behscore});

figure;
colorlines(plot(overlap_cyc'),cmap)

% -- average
normMethod = 'count';% 'count' or 'pdf'

% -- median split by task performance
[~,idx]=sort(sumdat{:,behscore});
grpvar = nan(N_ss,1);
grpvar(idx(1:N_ss/2)) = 1;% lower group
grpvar(idx(N_ss/2+1:end)) = 2;% higher group

[overlapMean,overlapSE]=grpstats(overlap_cyc,grpvar,{'mean','sem'});

% -- plot comparison
figure
plotSpectraSEM(ctrs,overlapMean,overlapSE,'cmap',[0,0,1;1,0,0]);
% title(CMEtasknames(taskIdx))
xlabel('cycle length')
ylabel('average % overlap')
%% ===== task TR distribution on cycles ===== %%
% N_bins = 21;
% edges = linspace(1.5,85.5,N_bins+1);
binwidth = 5;
maxCycLen = 55;
edges = 1.5:binwidth:maxCycLen+0.5;
N_bins = length(edges)-1;
ctrs = edge2ctr(edges);
% :: whether or not to exclude hubs themselves
excludeHubs = true; 

fTR_cyc = nan(N_ss,N_bins,N_tasks);

for nss=1:N_ss
    for nbin=1:N_bins
        % -- extract cycles within a range of cycle length
        cycidx = (sumdat.highdegCycLen{nss}>edges(nbin)) ...
            & (sumdat.highdegCycLen{nss}<edges(nbin+1));
        % -- find unique nodes in the set of cycles
        nodesInBin = unique(cell2mat(sumdat.allcycles{nss}(sumdat.highdegCycIdx{nss}(cycidx))'));
        if excludeHubs
            nodesInBin = setdiff(nodesInBin, sumdat.highdegidx{nss});
        end
        % -- find unique TRs in the nodes
        TRsInBin = unique(cell2mat(sumdat.members{nss}(nodesInBin)));
        % -- find TR fraction for each task
        fTR_cyc(nss,nbin,:) = mean(x_label(TRsInBin) == 1:N_tasks);
    end
end

m = squeeze(nanmean(fTR_cyc,1));
% figure
% colorlines(plot(ctrs,m'),CMEcmap)
% legend(CMEtasknames)

% -- plot group average
figure
b=bar(ctrs,m,'stacked');
for n=1:N_tasks
    b(n).FaceColor=CMEcmap(n);
end
legend(CMEtasknames,'location','northeastoutside')
xlabel('cycle length')
ylabel('TR fraction')
if excludeHubs
    title('hubs exlucded')
else
    title('hubs included')
end

% -- plot individual 
taskIdx = 3;
figure
behscore = 'AveragePC';
cmap=rankcolor(jet(N_ss),sumdat{:,behscore});
hold on
for n=1:N_ss
    plot(ctrs,squeeze(fTR_cyc(n,:,taskIdx)),'color',cmap(n,:))
end

%% ===== fraction TR by performance ===== %%
taskIdx = 5;
behscore = 'AveragePM';
% -- create grouping variable
grpvar=nan(N_ss,1);
grpvar(rankval(sumdat{:,behscore})<=N_ss/2)=1;
grpvar(grpvar~=1)=2;
% -- group stats
Nss_thres = 5;% a valid bin must have a minimal # subjects
validBinIdx = sum(~isnan(fTR_cyc(:,:,taskIdx)),1)>Nss_thres;
N_validBin = sum(validBinIdx);
Y = squeeze(fTR_cyc(:,validBinIdx,taskIdx));
Y(isnan(Y))=0;% optional: set "nan" to "0"
[fTRMean,fTRSE]=grpstats(Y,grpvar,{'mean','sem'});

% -- plot comparison
figure
plotSpectraSEM(ctrs(validBinIdx),fTRMean,fTRSE,'cmap',[0,0,1;1,0,0]);
% title(CMEtasknames(taskIdx))
xlabel('cycle length')
ylabel([CMEtasknames(taskIdx) 'TR fraction'])
ha=gca;
% -- statistical testing (ANOVA)
[grp_bin, grp_beh]=meshgrid(1:N_validBin,grpvar);
figure
[p,tbl,stats,terms]=anovan(Y(:),[grp_beh(:), grp_bin(:)],...
    'model','full','varnames',{behscore,'bin index'});%,'nested',[0 1; 0 0]);
[c,m,h,gnames]=multcompare(stats,'Dimension',1,'Display','off');
disp(['p=' num2str(c(end))])
% [c,m,h,gnames]=multcompare(stats,'Dimension',2);
[grpmean,grpse]=grpstats(Y(:),grp_beh(:),{'mean','sem'})
[c,m,h,gnames]=multcompare(stats,'Dimension',[1 2]);
% :: index of within-bin comparisons
compidx=ismember(c(:,1:2), [(1:2:N_validBin*2-1)',(2:2:N_validBin*2)'],'row');
% :: plot significance levels
addStars(ha,ctrs, max(fTRMean) + 0.04, c(compidx,end))
title(legend(ha,'SE(-)','mean(-)','SE(+)','mean(+)','location','best'),behscore)

% -- ttests
[h,p,ci,stats]=ttest2(Y(grpvar==1,:),Y(grpvar==2,:));