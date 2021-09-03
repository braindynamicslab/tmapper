% construction.m
%{
explain the procedures and results of construction of tmapper from
simulated data. 

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
%% ===== display raw data using MDS ===== %%
x_proj=cmdscale(pdist(x),2);
figure
scatter(x_proj(:,1),x_proj(:,2),50,x_label,'filled')
colormap jet
axis equal
axis off
%% ===== construct simplified graph
k = 14;
d = 10;
% g = tknndigraph (x,k,tidx,'reciprocal',true);
g = tknndigraph (x,k,tidx,'reciprocal',true,'timeExcludeSpace', true);

tic
[g_simp, members, ~] = filtergraph(g,d,'reciprocal',true);
toc
% -- plotting
[a1, a2,~,~,h] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,15]);
% [a1, a2, ~,~,h] = plotgraphtcm(g,x_label,t,[],'nodesizerange',[3,10]);
title(a1,['k=' num2str(k) ', d=' num2str(d)]);
title(a2,'geodesic recurrence plot')

layout(h,'force','Iterations',200)
h.ArrowPosition=0.7;
h.ArrowSize=6;
h.EdgeAlpha=0.5;