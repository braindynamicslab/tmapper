% filterConstruct.m
% -----------------
% construct filters for mapper

addpath('../')

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
%% ===== prepare data
Nds = 720; % step of down sampling
x = a.X(1:Nds:end,1:a.N); % subsampled data
% x = B.B(1:Nds:end,:); % BOLD
t = a.t(1:Nds:end); % subsampled time
tidx = (1:length(t))'; % index of t, s.t. could be used to en
x_label = state_regime(1:Nds:end); % regime of attractor
distMat = pdist2(x,x);
%% ===== construct nearest neighbor graph
% [~, mapping] = bdl_isomap(x,2,num_k);
k = 6;700;5;
delta = [];5;
g = tknndigraph (x,k,tidx,'reciprocal',true);
g = tknngraph (x,k,tidx,'reciprocal',true);
D_geo = distances(g);
D_geo_norm = D_geo/max(D_geo(:));

[Div,Conv] = localDivConv(g.adjacency,x,tidx);
c = zeros(length(t),1);
% c(landmk_idx) = 1;
figure;
plot(g,'EdgeAlpha',0.3,'EdgeColor','k','NodeCData',x_label,'Layout','force','MarkerSize',5); 
axis equal
axis off
cb=colorbar;
cb.Label.String = 'attractor index';
colormap(gca, 'jet')
title(sprintf('$w_{EE}=%g, w_{EI}=%g, k=%g, \\delta=%g$',...
    a.w_EE,a.w_EI,k,delta),'Interpreter','latex')

% -- Ball Mapper
resolution = 100;
gain = 26;
[landmk_idx, dist_to_landmk, epsilon] = px_maxmin(D_geo_norm,'metric',resolution,'n');

%Add points to bigBins

pts_in_bigBin = cell(resolution,1);
for ii = 1:resolution
    tmp = dist_to_landmk(ii,:);
    pt_idx = tmp < max((gain/100)*(4*epsilon), eps);
    % Gain should be at least 25 to make sure we're covering the dataset
    pts_in_bigBin{ii} = find(pt_idx);
end

idx_bigBin = landmk_idx;
Adja = bdl_create_cluster_graph_v0(pts_in_bigBin);
[~, idx] = sort(cell2mat(cellfun(@numel,pts_in_bigBin,'UniformOutput',0)));
attractor_size_rank=[];
attractor_size_rank(idx,1) = (1:length(Adja));
figure;plot(graph(Adja),'EdgeAlpha',0.3,'EdgeColor','k',...
    'NodeCData',cell2mat(cellfun(@(x) mode(x_label(x)),pts_in_bigBin,'UniformOutput',0)),...
    'Layout','force','MarkerSize',attractor_size_rank/max(attractor_size_rank)*20)
colormap jet

% -- Mapper
[y,E] = cmdscale((D_geo_norm+D_geo_norm')/2);
sl_histo_bins = 20;
pts_in_bigBin = bdl_binning_v0([normrange(Div) normrange(Conv)],10,50);
% pts_in_bigBin = bdl_binning_v0(y(:,1:2),20,70);
% pts_in_smallBin = bdl_pc_sl_hist_v0(pts_in_bigBin,distMat, ...
%                             sl_histo_bins);
pts_in_smallBin = bdl_pc_sl_hist_v0(pts_in_bigBin,D_geo_norm, ...
                            sl_histo_bins);
pts_in_smallBin_pruned = bdl_pruning_v0(pts_in_smallBin);     

adja        = bdl_create_cluster_graph_v0(pts_in_smallBin);
adja_pruned = bdl_create_cluster_graph_v0(pts_in_smallBin_pruned);

[~, idx] = sort(cell2mat(cellfun(@numel,pts_in_smallBin_pruned,'UniformOutput',0)));
attractor_size_rank=[];
attractor_size_rank(idx,1) = (1:length(adja_pruned));
figure;plot(graph(adja_pruned),'EdgeAlpha',0.3,'EdgeColor','k',...
    'NodeCData',cell2mat(cellfun(@(x) mode(x_label(x)),pts_in_smallBin_pruned,'UniformOutput',0)),...
    'Layout','force','MarkerSize',attractor_size_rank/max(attractor_size_rank)*20)
colormap jet

