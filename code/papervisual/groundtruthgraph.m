% groundtruthgraph.m
%{
plot ground truth
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
%% ===== adjust edge weights
g0_=g0;
msize =h.MarkerSize;
sweight= msize(findnode(g0_,g0_.Edges.EndNodes(:,1)))'; % weight of source node
tweight= msize(findnode(g0_,g0_.Edges.EndNodes(:,2)))'; % weight of target node
g0_.Edges.Weight = sweight+tweight;
[a1, a2,~,~,h_] = plotgraphtcm(g0_,x_label,t,nodemembers,'nodesizerange',[1,20]);
h_.ArrowPosition=0.7;
h_.ArrowSize=6;
h_.EdgeAlpha=0.5;
% h_.ArrowPosition=sweight./(sweight+tweight);