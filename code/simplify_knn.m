% simplify_knn.m
% ---------------
%{
Simplify knn graphs to have fewer nodes.
The objective is to use as few nodes as possible will representing the
correct topology.

%}
addpath(genpath('../'))
addpath(genpath('../../GWTLBgpu/'))
addpath('plot_utils','tmapper_tools','other_utils')
%% ===== load data
% clear all
% close all
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
distMat = pdist2(x,x);

% % -- delay embedding
% x = [x(1:end-1,:),x(2:end,:)];% delay embedding
% t = t(1:end-1);
% Nt = Nt - 1;
% tidx = (1:Nt)';
% x_label = x_label(1:end-1);
%% ===== construct simplified graph
k = 16;
d = 10;
% g = tknndigraph (x,k,tidx,'reciprocal',true);
g = tknndigraph (x,k,tidx,'reciprocal',true,'timeExcludeSpace', true);

tic
[g_simp, members, ~] = filtergraph(g,d,'reciprocal',true);
toc
% -- plotting
[a1, a2] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,15]);
% [a1, a2] = plotgraphtcm(g,x_label,t,[],'nodesizerange',[1,10]);
title(a1,['k=' num2str(k) ', d=' num2str(d)]);
title(a2,'geodesic recurrence plot')
%% ===== plot mds
x_proj=cmdscale(pdist(x),2);
figure;scatter(x_proj(:,1), x_proj(:,2),50,x_label,'filled')
colormap jet
axis equal
%% ===== compute ground truth
[g0,dwelltime,nodemembers] = symDyn2digraph(x_label); % a priori digraph
% -- geodesics and TCM
geod0 = distances(g0);
tcm0 = TCMdistance(g0,nodemembers); % geodesic distance bw time points on the directed graph
% -- plotting
[a1, a2,~,~,h] = plotgraphtcm(g0,x_label,t,nodemembers,'nodesizerange',[1,15]);
title(a1,'ground truth');
% -- normalize distance and node measures
[geod0_n, mNode0] = normgeo(geod0,nodesize(nodemembers));
tcm0_n = normtcm(tcm0);
%% ===== compute similarity to ground truth ===== %%
% -- parameter ranges
k_opts = 3:2:100;
d_opts = 2:50;

N_ks = length(k_opts);
N_ds = length(d_opts);

% -- store comparison
diff2geod0 = nan(N_ds, N_ks);
% -- metric for comparison
met = 'GW';'L_2'; 
for nk = 1:N_ks
    tic
    % knn graph
    g = tknndigraph (x,k_opts(nk),tidx,'reciprocal',true);
    for nd = 1:N_ds
        disp([nk nd])
        % -- simplified graph
        [g_simp, members, nsize] = filtergraph(g,d_opts(nd),'reciprocal',true);
        % -- calculated and normalized TCM or geodesic distance
        switch met
            case 'L_2'
                geod = TCMdistance(g_simp, members);
                % -- normalize the matrix
                geod_n = normtcm(geod);%normgeo(geod);
            case 'GW'
                geod = distances(g_simp,'Method','unweighted');
                % -- normalize the matrix
                [geod_n, mNode] = normgeo(geod,nsize);
        end
        
        % -- compare to ground truth
        switch met 
            case 'L_2'
                diff2geod0(nd,nk) = norm(geod_n(:)-tcm0_n(:));
            case 'GW'
                diff2geod0(nd,nk) = emd2RTLB_hetero(geod_n, geod0_n, mNode, mNode0);
        end
    end
    toc
end

figure;
imagesc(k_opts, d_opts,diff2geod0)
cb = colorbar;
cb.Label.String = [met ' distance'];
xlabel('k')
ylabel('d')
colormap hot
title('dissimilarity to ground truth')

%% ===== compute similarity to previous k and d ===== %%
% -- parameter ranges
k_opts = 3:2:100;
d_opts = 2:50;

N_ks = length(k_opts);
N_ds = length(d_opts);

% -- store data
[GEOD_N,MNODES] = deal(cell(N_ds, N_ks));
% -- store comparison
[diff2lastk,diff2lastd] = deal(nan(N_ds, N_ks));
% -- metric for comparison
met = 'GW';'L_2'; 
for nk = 1:N_ks
    disp(nk)
    tic
    % knn graph
    g = tknndigraph (x,k_opts(nk),tidx,'reciprocal',true);
    for nd = 1:N_ds
        % -- simplified graph
        [g_simp, members, nsize] = filtergraph(g,d_opts(nd),'reciprocal',true);
        % -- normalize geodesics
        switch met
            case 'L_2'
                geod = TCMdistance(g_simp, members);
                % -- normalize the matrix
                geod_n = geod/max(geod(:));%normgeo(geod);
            case 'GW'
                geod = distances(g_simp,'Method','unweighted');
                % -- normalize the matrix
                [geod_n, mNode] = normgeo(geod,nsize);
        end
        % -- store geodesics
        GEOD_N{nd,nk} = geod_n;
        MNODES{nd,nk} = mNode;
        % -- compare to previous k
        if nk>1
            switch met 
                case 'L_2'
                    diff2lastk(nd,nk) = norm(geod_n(:)-GEOD_N{nd,nk-1}(:));
                case 'GW'
                    diff2lastk(nd,nk) = emd2RTLB_hetero(geod_n, GEOD_N{nd,nk-1}, mNode, MNODES{nd,nk-1});
            end
        end
        % -- compare to previous d
        if nd>1
            switch met 
                case 'L_2'
                    diff2lastd(nd,nk) = norm(geod_n(:)-GEOD_N{nd-1,nk}(:));
                case 'GW'
                    diff2lastd(nd,nk) = emd2RTLB_hetero(geod_n, GEOD_N{nd-1,nk}, mNode, MNODES{nd-1,nk});
            end
        end
    end
    toc
end

figure;
imagesc(k_opts, d_opts,diff2lastk)
cb = colorbar;
cb.Label.String = [met ' distance'];
xlabel('k')
ylabel('d')
colormap hot
title('dissimilarity to (k-1)')

figure;
imagesc(k_opts, d_opts,diff2lastd)
cb = colorbar;
cb.Label.String = [met ' distance'];
xlabel('k')
ylabel('d')
colormap hot
title('dissimilarity to (d-1)')