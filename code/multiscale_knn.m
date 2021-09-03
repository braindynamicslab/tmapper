% multiscale_knn.m
% ----------------
% compute knn-graph, temporal knn-graph and temporal knn-digraph for
% different values of k, to see 
%  1. how the topology recovered by TDA approaches that of the ground truth.
%  2. how the topology of knn-graphs/digraphs change with respect to itself
%  across k's.
%  3. whether 1 and 2 are related, so that we can infer an appropriate k
%  without prior knowledge of the true topology of the dynamics.
%{
created by MZ, 8-30-2019

%}

addpath('../')
addpath(genpath('../../GWnets/'))
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
x = a.X(1:Nds:end,1:a.N); % subsampled data (S_E)
% x = B.B(1:Nds:end,:); % subsampled BOLD
t = a.t(1:Nds:end); % subsampled time
Nt = length(t);
tidx = (1:Nt)'; % index of t, s.t. could be used to en
x_label = state_regime(1:Nds:end); % regime of attractor
distMat = pdist2(x,x);

% -- compute ground truth
[g0,dwelltime,nodemembers] = symDyn2digraph(x_label); % a priori digraph
g0_undir = digraph2graph(g0); % undirected version of a priori graph
% direct geodesics
geod0 = distances(g0);
geod0_undir = distances(g0_undir);
mNode0 = nodesize(nodemembers);
tcm0 = TCMdistance(g0,nodemembers); % geodesic distance bw time points on the directed graph
tcm0_undir = TCMdistance(g0_undir,nodemembers); % geodesic distance bw time points on the undirected graph

% -- normalize distance and node measures
geod0_n = geod0/max(geod0(:));
geod0_undir_n = geod0_undir/max(geod0_undir(:));
mNode0 = mNode0/sum(mNode0); % prob distribtuion must sum to 1
tcm0_n = tcm0/max(tcm0(:));
tcm0_undir_n = tcm0_undir/max(tcm0_undir(:));
%% ===== compute knn graphs and compare to ground truth
% -- what type of graph  (knn, tknn,tknndi)
graphtype = 'tknndi';'knn';'tknn';
% -- whether to enforce reciprocity
rec = 1;

% -- define k - spatial resolutions
k_opts = 3:2:100;
% k_opts = 3:20:1000;
Nk = length(k_opts);

% -- determine metric
met = 'GW weighted';'TCM';'GW'; 'GW bw TCM'; 

% -- preallocate storage for graphs and geodesics
gs = cell(Nk,1); % graphs
[geods, geods_n] = deal(nan(Nt,Nt,Nk)); % raw & normalized geodesics
[diff2geod0, diff2geod0_undir,diff2lastk] = deal(nan(Nk,1)); % comparisons

for n = 1:Nk
    disp(n)
    tic
    
    % -- compute graph & geodesic distances
    switch graphtype
        case 'knn'
            gs{n} = eval(sprintf('%sgraph(x,k_opts(n),''reciprocal'',%g);',graphtype,rec));
        otherwise
            gs{n} = eval(sprintf('%sgraph(x,k_opts(n),tidx,''reciprocal'',%g);',graphtype, rec));
    end
    geod = distances(gs{n}); 
    geods(:,:,n) = geod;
    
    % -- check for inf (set to max)
    geod(geod==Inf) = max(geod(geod<Inf));
    
    % -- normalize geodesics
    geod_n = geod/max(geod(:));     
    geods_n(:,:,n) = geod_n;
    
    % -- compare to ground truth
    switch met
        case 'TCM'
            diff2geod0(n) = norm(geod_n(:)-tcm0_n(:));
            diff2geod0_undir(n) = norm(geod_n(:)-tcm0_undir_n(:));
        case 'GW bw TCM'
            diff2geod0(n) = emd2RTLB_unih(geod_n,tcm0_n);
            diff2geod0_undir(n) = emd2RTLB_unih(geod_n,tcm0_undir_n);
        case 'GW'
            diff2geod0(n) = emd2RTLB_unih(geod_n,geod0_n);
            diff2geod0_undir(n) = emd2RTLB_unih(geod_n,geod0_undir_n);
        case 'GW weighted'
            diff2geod0(n) = emd2RTLB_hetero(geod_n,geod0_n,ones(Nt,1)/Nt,mNode0);
            diff2geod0_undir(n) = emd2RTLB_hetero(geod_n,geod0_undir_n,ones(Nt,1)/Nt,mNode0);
    end
    
    % -- compare to previous k
    if n>1
        lastgeod_n = squeeze(geods_n(:,:,n-1));
        switch met
            case 'TCM'
                diff2lastk(n) = norm(geod_n(:)-lastgeod_n(:));
            case 'GW bw TCM'
                diff2lastk(n) = emd2RTLB_unih(geod_n, lastgeod_n);
            case 'GW'
                diff2lastk(n) = emd2RTLB_unih(geod_n, lastgeod_n);
            case 'GW weighted'%tknn graph is always unweighted
                diff2lastk(n) = emd2RTLB_unih(geod_n, lastgeod_n);
        end
    end
    
    toc
end

% -- plotting results
figure
hold on
% compre to graph truth
plot(k_opts,diff2geod0)
plot(k_opts,diff2geod0_undir)
% compare to previous k
plot(k_opts,diff2lastk)
xlabel('k')
ylabel('dissimilarity')
legend('vs. true directed graph','vs. true undirected graph','vs. (k-1)-nn graph')
switch graphtype
    case 'knn'
        title([met ': knn-graph'])
    case 'tknn'
        title([met ': temporal knn-graph (undirected)'])
    case 'tknndi'
        title([met ': temporal knn-graph (directed)'])
end
