% tmapper_parameter.m
%{
 null model-based parameter selection
%}

addpath(genpath('../'))
addpath(genpath('../../GWTLBgpu/'))
addpath('plot_utils','tmapper_tools','other_utils')
%% ===== load simulated data
clear all
close all
clc

% -- define local parameters
wee = 2.8;
wei = 1;

% -- load data
examplef = dir(['../results/varG_examples/' par2filename({'w_EE',wee,'w_EI',wei}) '*.mat']);
[~,examplefname] = fileparts(examplef.name);
load(fullfile(examplef.folder,examplef.name))

% -- subsample data
Nds = 720; % step of down sampling
x = a.X(1:Nds:end,1:a.N); % subsampled data (S_E)
% x = B.B(1:Nds:end,:); % subsampled BOLD
% x = phaseRand(x);
t = a.t(1:Nds:end); % subsampled time
Nt = length(t);
tidx = (1:Nt)'; % index of t, s.t. could be used to en
x_label = state_regime(1:Nds:end); % regime of attractor
%% ===== load human data
clear all
close all
clc

ss_id = 4; % subject id
load(sprintf('../data/CME/SBJ%02d_Shine_375.mat',ss_id))
%% ===== compute ground truth
[g0,dwelltime,nodemembers] = symDyn2digraph(x_label); % a priori digraph
% -- geodesics and TCM
geod0 = distances(g0);
tcm0 = TCMdistance(g0,nodemembers); % geodesic distance bw time points on the directed graph

% -- normalize distance and node measures
[geod0_n, mNode0] = normgeo(geod0,nodesize(nodemembers));
tcm0_n = normtcm(tcm0,'normtype','max');
%% ===== construct simplified graph
k = 19;
d = 2;
% g = tknndigraph (x,k,tidx,'reciprocal',true);
g = tknndigraph (x,k,tidx,'reciprocal',true,'timeExcludeSpace', true);

tic
[g_simp, members, ~] = filtergraph(g,d,'reciprocal',true);
toc
% -- plotting
[a1, a2] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,10]);
% [a1, a2] = plotgraphtcm(g,x_label,t,[],'nodesizerange',[1,10]);
title(a1,['k=' num2str(k) ', d=' num2str(d)]);
title(a2,'geodesic recurrence plot')
%% ===== null model ===== %%
% xnull = phaseRand(x,8);

k = 30;
d = 2;

g = tknndigraph (xnull,k,tidx,'reciprocal',true,'timeExcludeSpace', true);

tic
[g_simp, members, ~] = filtergraph(g,d,'reciprocal',true);
toc
% -- plotting
[a1, a2] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,10]);
% [a1, a2] = plotgraphtcm(g,x_label,t,[],'nodesizerange',[1,10]);
title(a1,['k=' num2str(k) ', d=' num2str(d)]);
title(a2,'geodesic recurrence plot')

%% ===== difference between real and null model ===== %%
met = 'L_2';
null_name = 'PR';
nullseed = 7;
kopts = 3:1:30;
d = 2;

N_kopts = length(kopts);
[diff2null, diff2truth,nulldiff2truth, hubcount,hubcountnull] = deal(nan(1,N_kopts));
[ADJ, ADJnull, TCM,TCMnull] = deal(cell(1,N_kopts));

rng(nullseed)
switch null_name
        case 'perm'
            % - permute time points - %
            rorder = randperm(Nt);
            xnull = x(rorder,:);
        case 'PR'
            % - or phase randomize - %
            xnull = phaseRand(x);
end
    
figure
hold on
for nk = N_kopts:-1:1
    disp(kopts(nk))
    % -- real model
    tic
    g = tknndigraph (x,kopts(nk),tidx,'reciprocal',true,'timeExcludeSpace', true);
    [g_simp, members, nsize] = filtergraph(g,d,'reciprocal',true);
	toc
    % -- null model
    tic
    g_null = tknndigraph (xnull,kopts(nk),tidx,'reciprocal',true,'timeExcludeSpace', true);
    [g_simp_null, members_null, nsize_null] = filtergraph(g_null,d,'reciprocal',true);
    toc
    % -- prepare distance matrices
    tic
    switch met
        case 'L_2'
            geod_n = normtcm(TCMdistance(g_simp, members),'normtype','max');
            geod_n_null = normtcm(TCMdistance(g_simp_null, members_null),'normtype','max');
            nanidx = isnan(geod_n(:)) | isnan(geod_n_null(:));
            diff2null(nk) = norm(geod_n(~nanidx)-geod_n_null(~nanidx));
            nanidx = isnan(geod_n(:)) | isnan(tcm0_n(:));
            diff2truth(nk) = norm(geod_n(~nanidx)-tcm0_n(~nanidx));
            nulldiff2truth(nk) = norm(geod_n_null(~nanidx)-tcm0_n(~nanidx));
        case 'GW'
            [geod_n, mNode] = normgeo(distances(g_simp,'Method','unweighted'),nsize);
            [geod_n_null, mNode_null] = normgeo(distances(g_simp,'Method','unweighted'),nsize);
            diff2null(nk) = emd2RTLB_hetero(geod_n, geod_n_null, mNode, mNode_null);
            diff2truth(nk) = emd2RTLB_hetero(geod_n, geod0_n, mNode, mNode0);
    end
    toc
       
    scatter(kopts(nk), diff2null(nk),'k','filled')
    scatter(kopts(nk), diff2truth(nk),'g','filled')
    scatter(kopts(nk), norm(geod_n(:)), 'b','filled')
    scatter(kopts(nk), norm(geod_n_null(:)), 'r','filled')
    drawnow
    
    % -- store graphs
    ADJ{nk} = g_simp.adjacency;
    ADJnull{nk} = g_simp_null.adjacency;
    TCM{nk}=geod_n;
    TCMnull{nk} = geod_n_null;
    % -- graph properties
%     degdistr{nk} = sort(sum(g_simp.adjacency,1)' + sum(g_simp.adjacency,2));
%     degdistrnull{nk} = sort(sum(g_simp_null.adjacency,1)' + sum(g_simp_null.adjacency,2));
    hubcount(nk) = sum((sum(g_simp.adjacency,1)' + sum(g_simp.adjacency,2))>2)./g_simp.numnodes;
    hubcountnull(nk) = sum((sum(g_simp_null.adjacency,1)' + sum(g_simp_null.adjacency,2))>2)./g_simp_null.numnodes;
end

legend('diff2null','diff2truth','norm','normnull')
xlabel('k')

% --- other heauristic metric
figure
plot(kopts,hubcount,'-o')
hold on
plot(kopts,hubcountnull,'-o')
plot(kopts,hubcount-hubcountnull,'-o')
legend('data','null','data-null')
xlabel('k')
%% ===== analyzing change of degree distribution ===== %%
degdistr = cellfun(@full,degdistr,'UniformOutput',0);
degdistrnull = cellfun(@full,degdistrnull,'UniformOutput',0);

figure;
plot(kopts,cellfun(@std, degdistr))
hold on
plot(kopts,cellfun(@std, degdistrnull))
%% ===== analyzing distribution entropy ===== %%
% -- degree distribution entropy of tmapper
DEG = cellfun(@(x) [sum(x), sum(x')]',ADJ,'UniformOutput',0);
DEGnull = cellfun(@(x) [sum(x), sum(x')]',ADJnull,'UniformOutput',0);
H_DEG = cellfun(@entropyDiscrete,DEG);
H_DEGnull = cellfun(@entropyDiscrete,DEGnull);
% -- average distance distribution entropy of TCM
AVGD = cellfun(@(x) [sum(x), sum(x')]',TCM,'UniformOutput',0);
AVGDnull = cellfun(@(x) [sum(x), sum(x')]',TCMnull,'UniformOutput',0);
% AVGD = cellfun(@(x) x(:),TCM,'UniformOutput',0);
% AVGDnull = cellfun(@(x) x(:),TCMnull,'UniformOutput',0);
H_AVGD = cellfun(@entropyDiscrete,AVGD);
H_AVGDnull = cellfun(@entropyDiscrete,AVGDnull);
% -- cycle length distribution
% tic
% [~,~,~,CYC]=cellfun(@CycleCount2p,ADJ,'UniformOutput',0);
% toc
% CYCLEN = cellfun(@(x) cellfun(@length,x),CYC,'UniformOutput',0);
% H_CYCLEN = cellfun(@entropyDiscrete,CYCLEN);

figure
plot(kopts, [H_DEG; H_DEGnull; H_AVGD; H_AVGDnull]')
legend('H(degree)','H(degree,null)','H(avgdist)','H(avgdist,null)')

% plot(kopts, [H_DEG; H_DEGnull; H_AVGD; H_AVGDnull; H_CYCLEN]')
% legend('H(degree)','H(degree,null)','H(avgdist)','H(avgdist,null)','H(cyclen)')

% -- compare degree-distribution entropy to diff to ground truth
figure
% hold on
yyaxis left
plot(kopts, diff2truth);
hold on
plot(kopts, nulldiff2truth,'--');
xlabel('k')
ylabel('disimilarity to ground truth (a.u.)')
title('simulated BOLD')
yyaxis right
plot(kopts, H_DEG);
hold on
plot(kopts, H_DEGnull);
ylabel('degree-distribution entropy')
legend('original',null_name,'original',null_name)


%% ===== degree distribution entropy for all subjects ===== %% 
clear all
close all
clc

datapath = '../data/CME/'; % path of fMRI data
sumdatpath = [datapath 'behavior_cme.xlsx']; % path of summary data
savepath = '../results/CME/'; % path for saving results

sumdat = readtable(sumdatpath); % behavioral and meta-data
N_ss = height(sumdat); % number of subjects

% -- compute degree entropy for a range of k
kopts = 3:1:30;
d = 2;
N_kopts = length(kopts);

H_DEG_all = nan(N_ss,N_kopts);

figure
hold on
for n = 1:N_ss
    % --- load data
    disp(sumdat.id(n));
    load(sprintf('../data/CME/%s_Shine_375.mat',sumdat.id{n})); 
    tic
    for nk = N_kopts:-1:1
        disp(kopts(nk))
        % -- construct shape graph
        g = tknndigraph (x,kopts(nk),tidx,'reciprocal',true,'timeExcludeSpace', true);
        [g_simp, members, nsize] = filtergraph(g,d,'reciprocal',true);
        
        % -- compute entropy
        H_DEG_all(n,nk) = entropyDiscrete([sum(g_simp.adjacency), sum(g_simp.adjacency)]');
    end
    toc
    
    % -- plot
    plot(kopts,H_DEG_all(n,:))
    drawnow
end

% -- save
save('results/parselection/H_deg.mat','H_DEG_all','sumdat','kopts','d','N_ss','N_kopts')

% -- a prettier plot
cmap = rankcolor(jet(N_ss),sum(H_DEG_all'));
figure
colorlines(plot(kopts,H_DEG_all),cmap)
xlabel('k')
ylabel('degree distribution entropy')
cb=colorbar;
caxis([1 N_ss])
colormap(jet(N_ss))
set(cb,'ytick',1:N_ss,'yticklabel',1:N_ss)
cb.Label.String = 'rank of mean entropy';
title('effect of k by subject')

% -- an average plot
[m,se]=grpstats(H_DEG_all,ones(N_ss,1),{'mean','sem'});
figure
plotSpectraSEM(kopts,m,se)
legend('SE','mean')
xlabel('k')
ylabel('degree distribution entropy')
title('effect of k averaged across subjects')
%% ===== degree distribution entropy for all subjects' null models ===== %% 
clear all
close all
clc

datapath = '../data/CME/'; % path of fMRI data
sumdatpath = [datapath 'behavior_cme.xlsx']; % path of summary data
savepath = '../results/CME/'; % path for saving results

sumdat = readtable(sumdatpath); % behavioral and meta-data
N_ss = height(sumdat); % number of subjects

% -- null model parameters
null_name = 'perm';
nullseed = 7;
rng(nullseed)

% -- compute degree entropy for a range of k
kopts = 3:1:30;
d = 2;
N_kopts = length(kopts);

H_DEG_all = nan(N_ss,N_kopts);

figure
hold on
for n = 1:N_ss
    % --- load data
    disp(sumdat.id(n));
    load(sprintf('../data/CME/%s_Shine_375.mat',sumdat.id{n})); 
    
    % --- create null model
    switch null_name
        case 'perm'
            % - permute time points - %
            rorder = randperm(Nt);
            xnull = x(rorder,:);
        case 'PR'
            % - or phase randomize - %
            xnull = phaseRand(x);
    end
    tic
    
    for nk = N_kopts:-1:1
        disp(kopts(nk))
        % -- construct shape graph
        g = tknndigraph (xnull,kopts(nk),tidx,'reciprocal',true,'timeExcludeSpace', true);
        [g_simp, members, nsize] = filtergraph(g,d,'reciprocal',true);
        
        % -- compute entropy
        H_DEG_all(n,nk) = entropyDiscrete([sum(g_simp.adjacency), sum(g_simp.adjacency)]');
    end
    toc
    
    % -- plot
    plot(kopts,H_DEG_all(n,:))
    drawnow
end

% -- save
save('results/parselection/H_deg.mat','H_DEG_all','sumdat','kopts','d','N_ss','N_kopts')

% -- a prettier plot
cmap = rankcolor(jet(N_ss),sum(H_DEG_all'));
figure
colorlines(plot(kopts,H_DEG_all),cmap)
xlabel('k')
ylabel('degree distribution entropy')
cb=colorbar;
caxis([1 N_ss])
colormap(jet(N_ss))
set(cb,'ytick',1:N_ss,'yticklabel',1:N_ss)
cb.Label.String = 'rank of mean entropy';
title(['effect of k by subject (' null_name ')'])

% -- an average plot
[m,se]=grpstats(H_DEG_all,ones(N_ss,1),{'mean','sem'});
figure
plotSpectraSEM(kopts,m,se)
legend('SE','mean')
xlabel('k')
ylabel('degree distribution entropy')
title(['effect of k averaged across subjects (' null_name ')'])

% -- correlation of total entropy with behavior
[r,p]=corr(sum(H_DEG_all')',sumdat.AveragePC);
disp([r p])
[r,p]=corr(sum(H_DEG_all')',sumdat.AverageRT);
disp([r p])
[r,p]=corr(sum(H_DEG_all')',sumdat.AveragePM);
disp([r p])