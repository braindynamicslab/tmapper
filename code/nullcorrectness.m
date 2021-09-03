% nullcorrectness.m
% -----------------
% calculate null distributions of shape graph correctness score (i.e.
% dissimilarity to ground truth measured by L2 or GW).
% Here we try different types of null models.

addpath(genpath('../'))
addpath(genpath('../../GWTLBgpu/'))
%% ===== load data
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
% x = a.X(1:Nds:end,1:a.N); % subsampled data (S_E)
x = B.B(1:Nds:end,:); % subsampled BOLD
t = a.t(1:Nds:end); % subsampled time
Nt = length(t);
tidx = (1:Nt)'; % index of t, s.t. could be used to en
x_label = state_regime(1:Nds:end); % regime of attractor
%% ===== compute ground truth
[g0,dwelltime,nodemembers] = symDyn2digraph(x_label); % a priori digraph
% -- geodesics and TCM
geod0 = distances(g0);
tcm0 = TCMdistance(g0,nodemembers); % geodesic distance bw time points on the directed graph

% -- normalize distance and node measures
[geod0_n, mNode0] = normgeo(geod0,nodesize(nodemembers));
tcm0_n = normtcm(tcm0);
%% ===== permutation of time (variable k,d) ===== %%
% -- options
null_name = 'PR';
% -- parameters that work for the original time series
k = 14;
d = 10;

N_perm = 200;% number of random permutations

% -- range of k, d
k_opts = 2:2:k;
d_opts = 2:2:d;

N_ks = length(k_opts);
N_ds = length(d_opts);

% -- preallocate storage for correctness scores
[diff2geod0_GW,diff2geod0_L2] = deal(nan(N_perm,N_ks,N_ds));

% -- permute, construct graph, and scoring
for n = 1:N_perm
    disp(n)
    tic
    
    switch null_name
        case 'perm'
            % - permute time points - %
            rorder = randperm(Nt);
            x_ = x(rorder,:);
        case 'PR'
            % - or phase randomize - %
            x_ = phaseRand(x);
    end
    
    for nk = 1:N_ks
        % -- construct knn graph
        g = tknndigraph (x_,k_opts(nk),tidx,'reciprocal',true);
        for nd = 1:N_ds
            % -- simplify graph
            [g_simp, members, nsize] = filtergraph(g,d_opts(nd),'reciprocal',true);
            
            % ----- compute correctness using L2 ----- %
            geod = TCMdistance(g_simp, members);
            % -- normalize the matrix
            geod_n = geod/max(geod(:));%normgeo(geod);
            % -- compare to ground truth
            diff2geod0_L2(n,nk,nd) = norm(geod_n(:)-tcm0_n(:));
            
            % ----- computing correctness using GW ----- %
            geod = distances(g_simp,'Method','unweighted');
            % -- normalize the matrix
            [geod_n, mNode] = normgeo(geod,nsize);
            % -- compare to ground truth
            diff2geod0_GW(n,nk,nd) = emd2RTLB_hetero(geod_n, geod0_n, mNode, mNode0);
            
            % -- skip over larger d's if there's only one node
            if g_simp.numnodes < 2
                break;
            end
        end
    end
    toc
end

% -- find min dissimilarity to ground truth 
diff2geod0_bestL2 = nanmin(nanmin(diff2geod0_L2,[],3),[],2);
diff2geod0_bestGW = nanmin(nanmin(diff2geod0_GW,[],3),[],2);

% -- plotting
% [a1, a2] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,25]);
% title(a1,['k=' num2str(k) ', d=' num2str(d)]);

% -- saving data 
cName = 'HCP11avg';
mkdir(['../data/' cName '/'])
save(['../data/' cName '/null_' null_name '_' par2filename({'w_EE',wee,'w_EI',wei}) '_BOLD.mat'],...
    'diff2geod0_bestL2','diff2geod0_bestGW','diff2geod0_L2','diff2geod0_GW',...
    'k_opts','d_opts','N_ks','N_ds')
%% ===== plotting null distribution and actual performance
% ----- dissimilarity to ground truth of the original time series
k = 14;
d = 10;
g = tknndigraph (x,k,tidx,'reciprocal',true);
[g_simp, members, nsize] = filtergraph(g,d,'reciprocal',true);

% ----- compare to null distribution
nbins = 20;
% -- using L2
CI = confint(diff2geod0_bestL2);
% compare to ground truth
geod = TCMdistance(g_simp, members);
geod_n = normtcm(geod);%normgeo(geod);
diff2geod0 = norm(geod_n(:)-tcm0_n(:));
% plot
figure
histogram(diff2geod0_bestL2,nbins,'Normalization','pdf'); hold on
plotevent(CI)
plotevent(diff2geod0,'g')
xlabel('dissimilarity to ground truth')
ylabel('pdf')
title('L2')
% -- using GW
CI = confint(diff2geod0_bestGW);
% compare to ground truth
geod = distances(g_simp,'Method','unweighted');
[geod_n, mNode] = normgeo(geod,nsize);
diff2geod0 = emd2RTLB_hetero(geod_n, geod0_n, mNode, mNode0);
% plot
figure
histogram(diff2geod0_bestGW,nbins); 
plotevent(CI)
plotevent(diff2geod0,'g')
xlabel('dissimilarity to ground truth')
ylabel('pdf')
title('GW')
%% ===== permutation of time (matched k,d) ===== %%
% -- options
null_name = 'perm';
% -- parameters that work for the original time series
k = 30;
d = 2;

N_perm = 1000;% number of random permutations

% -- range of k, d
k_opts = 3:1:k;
d_opts = 2:2:d;

N_ks = length(k_opts);
N_ds = length(d_opts);

% -- preallocate storage for correctness scores
[diff2geod0_GW,diff2geod0_L2] = deal(nan(N_perm,N_ks,N_ds));

% -- permute, construct graph, and scoring
for n = 1:N_perm
    tic
    
    switch null_name
        case 'perm'
            % - permute time points - %
            rorder = randperm(Nt);
            x_ = x(rorder,:);
        case 'PR'
            % - or phase randomize - %
            x_ = phaseRand(x);
    end
    
    for nk = 1:N_ks
        disp([n k_opts(nk)])
        % -- construct knn graph
        g = tknndigraph (x_,k_opts(nk),tidx,'reciprocal',true);
        for nd = 1:N_ds
            % -- simplify graph
            [g_simp, members, nsize] = filtergraph(g,d_opts(nd),'reciprocal',true);
            
            % ----- compute correctness using L2 ----- %
            geod = TCMdistance(g_simp, members);
            % -- normalize the matrix
            geod_n = normtcm(geod);
            % -- compare to ground truth
            diff2geod0_L2(n,nk,nd) = norm(geod_n(:)-tcm0_n(:));
            
            % ----- computing correctness using GW ----- %
            geod = distances(g_simp,'Method','unweighted');
            % -- normalize the matrix
            [geod_n, mNode] = normgeo(geod,nsize);
            % -- compare to ground truth
            diff2geod0_GW(n,nk,nd) = emd2RTLB_hetero(geod_n, geod0_n, mNode, mNode0);
            
            disp([diff2geod0_L2(n,nk,nd) diff2geod0_GW(n,nk,nd)])
        end
    end
    toc
end

% -- find min dissimilarity to ground truth 
diff2geod0_bestL2 = nanmin(nanmin(diff2geod0_L2,[],3),[],2);
diff2geod0_bestGW = nanmin(nanmin(diff2geod0_GW,[],3),[],2);

% -- plotting
% [a1, a2] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,25]);
% title(a1,['k=' num2str(k) ', d=' num2str(d)]);

% -- saving data 
cName = 'null';
mkdir(['data/' cName '/'])
save(['data/' cName '/null_' null_name '_' par2filename({'w_EE',wee,'w_EI',wei,'k',k,'d',d}) '_BOLD.mat'],...
    'diff2geod0_bestL2','diff2geod0_bestGW','diff2geod0_L2','diff2geod0_GW',...
    'k_opts','d_opts','N_ks','N_ds')
%% ===== original time (varying k,d) ===== %%
% -- parameters that work for the original time series
k = 30;
d = 2;

% -- range of k, d
k_opts = 3:1:k;
d_opts = 2:2:d;

N_ks = length(k_opts);
N_ds = length(d_opts);

% -- preallocate storage for correctness scores
[diff2geod0_GW,diff2geod0_L2] = deal(nan(N_ks,N_ds));

% -- construct graph, and scoring
tic
for nk = 1:N_ks
    disp([k_opts(nk)])
    % -- construct knn graph
    g = tknndigraph (x,k_opts(nk),tidx,'reciprocal',true);
    for nd = 1:N_ds
        % -- simplify graph
        [g_simp, members, nsize] = filtergraph(g,d_opts(nd),'reciprocal',true);

        % ----- compute correctness using L2 ----- %
        geod = TCMdistance(g_simp, members);
        % -- normalize the matrix
        geod_n = normtcm(geod);
        % -- compare to ground truth
        diff2geod0_L2(nk,nd) = norm(geod_n(:)-tcm0_n(:));

        % ----- computing correctness using GW ----- %
        geod = distances(g_simp,'Method','unweighted');
        % -- normalize the matrix
        [geod_n, mNode] = normgeo(geod,nsize);
        % -- compare to ground truth
        diff2geod0_GW(nk,nd) = emd2RTLB_hetero(geod_n, geod0_n, mNode, mNode0);

        disp([diff2geod0_L2(nk,nd) diff2geod0_GW(nk,nd)])
    end
end
toc

% -- find min dissimilarity to ground truth 
diff2geod0_bestL2 = nanmin(nanmin(diff2geod0_L2,[],2),[],1);
diff2geod0_bestGW = nanmin(nanmin(diff2geod0_GW,[],2),[],1);

% -- plotting
% [a1, a2] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,25]);
% title(a1,['k=' num2str(k) ', d=' num2str(d)]);

% -- saving data 
cName = 'real';
mkdir(['data/' cName '/'])
% save(['data/' cName '/real_' par2filename({'w_EE',wee,'w_EI',wei,'k',k,'d',d}) '_SE.mat'],...
%     'diff2geod0_bestL2','diff2geod0_bestGW','diff2geod0_L2','diff2geod0_GW',...
%     'k_opts','d_opts','N_ks','N_ds')
%% ===== plotting the results (multiple k) ===== %%
nullname = 'PR';
datatype = 'SE';
measure = 'GW';

% -- max k, d explored
k=30;
d=2;

nulldat=load(['data/null/null_' nullname '_' par2filename({'w_EE',wee,'w_EI',wei,'k',k,'d',d}) '_' datatype '.mat']);
realdat=load(['data/real/real_' par2filename({'w_EE',wee,'w_EI',wei,'k',k,'d',d}) '_' datatype '.mat']);

siglevel = 0.05;
CI=confint(eval("nulldat.diff2geod0_" + measure),siglevel);
realdis = eval("realdat.diff2geod0_" + measure);

figure;
hold on
% -- plot CI
nullcolor = [1 0 0];
nullalpha = 0.3;
realcolor = [0 1 0];
plotColorBand(nulldat.k_opts,CI(1,:),CI(2,:),nullcolor,nullalpha)
plot(nulldat.k_opts, mean(eval("nulldat.diff2geod0_" + measure)),'color',nullcolor)
plot(realdat.k_opts,realdis,'color',realcolor);
legend(sprintf('null CI (p=%g)',siglevel),'null mean', 'real','location','best')
xlabel('k')
ylabel(['dissimilarity to ground truth (' measure ')'])
title([nullname ' of ' datatype])

outdir='results/nullCI/';
mkdir(outdir);

fname = [outdir datatype '_' nullname '_' measure '_' par2filename({'w_EE',wee,'w_EI',wei,'k',k,'d',d})];
% print([fname '.png'],'-dpng')
% print([fname '.eps'],'-depsc','-painters')
% savefig([fname '.fig'])
%% ===== plotting the results (single k) ===== %%
nullname = 'PR';
datatype = 'SE';
measure = 'L2';

% -- max k, d explored
k=30;
d=2;

nulldat=load(['data/null/null_' nullname '_' par2filename({'w_EE',wee,'w_EI',wei,'k',k,'d',d}) '_' datatype '.mat']);
realdat=load(['data/real/real_' par2filename({'w_EE',wee,'w_EI',wei,'k',k,'d',d}) '_' datatype '.mat']);

% -- k chosen as optimal
switch datatype
    case 'SE'
        k_choice = 16;
    case 'BOLD'
        k_choice = 14;
end
k_choice_idx = find(nulldat.k_opts==k_choice);

% -- define confidence interval
siglevel = 0.05;
CI=confint(eval("nulldat.diff2geod0_" + measure),siglevel);
realdis = eval("realdat.diff2geod0_" + measure);

% -- define histogram
nbins = 200;
switch nullname
    case 'perm'
        switch measure
            case 'L2'
                xl = [150 550]; 
            case 'GW'
                xl = [0.08 0.28];
        end
        yl = [0 0.65];
    case 'PR'
        switch measure
            case 'L2'
                xl = [150 700]; 
            case 'GW'
                xl = [0.08 0.7];
        end
        yl = [0 1];
end
edges=linspace(xl(1),xl(2),nbins+1);
ctrs = edge2ctr(edges);

% -- other plot parameters

nullcolor = [1 0 0];
nullalpha = 0.3;
realcolor = [0 1 0];

% -- plotting
figure('position', [1221,635,560,160]);
hold on
% :: confidence interval
patch([CI(:,k_choice_idx); flipud(CI(:,k_choice_idx))]', repelem(yl,1,2),nullcolor,...
    'edgecolor','none','facealpha',nullalpha)
% :: null distribution
histogram(eval("nulldat.diff2geod0_" + measure + "(:,k_choice_idx)"),edges,...
    'edgecolor','none','facecolor','k','facealpha',0.5,'Normalization','probability')
% :: real value
line(repmat(eval("realdat.diff2geod0_" + measure + "(k_choice_idx,:)"),2,1), yl,...
    'color',realcolor,'linewidth',2)
xlim(xl)
ylim(yl)
ylabel('probability')
xlabel("dissimilarity (" + measure + ")")
title(nullname + " " + datatype)
% -- plot CI
nullcolor = [1 0 0];
nullalpha = 0.3;
realcolor = [0 1 0];


outdir='results/nullCI/';
mkdir(outdir);

fname = [outdir datatype '_' nullname '_' measure '_' par2filename({'w_EE',wee,'w_EI',wei,'k',k,'d',d,'kc',k_choice})];
print([fname '.png'],'-dpng')
print([fname '.eps'],'-depsc','-painters')
savefig([fname '.fig'])