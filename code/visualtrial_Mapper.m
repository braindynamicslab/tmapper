% visualtrial_Mapper.m
% --------------------
% visualize simulated dynamics using mapper.
% (8-14-2019)

addpath('../')
addpath('2019-neuromapper/code/mapperCode/bdl_mapper_v0')
addpath(genpath('2019-neuromapper/code/dimReducMethods/'));
addpath(genpath('2019-neuromapper/code/tools/'));
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
x = a.X(1:Nds:end,1:a.N);
xlabel = state_regime(1:Nds:end);
distMat = pdist2(x,x);
%% ===== mapper options
resdir              = '../results/';
mkdir(resdir)
options.filename    = examplefname; 
options.save_to     = [resdir,'res_',options.filename];

options.dimreduce   = 'bdl_isomap';
options.dim_embed   = 2;
options.filter      = 'identity'; %'time', 'identity' 

% -- binning parameters
options.binning     = 'ball';
options.resolution  = 30; 
options.gain        = 50;
options.clustering  = 'sl_histo'; %single linkage with histogram
options.sl_histo_bins = 1;

% -- For methods using knn, specify the k value
% Used for 'bdl_isomap', 'Isomap', 'LLE', among others
% Also used for ballMapper
options.knnparam    = 30; % for trefoil, if winding is w, try 3*w 

% For iterative methods, specify max number of iterations
options.maxiter     = 200;

% Color parameter
options.colors = xlabel;

%-- Run mapper on data

res = bdl_mapper_v0(x,distMat,options);
%% ===== alternative plotting
% -- node assignment
nodecolor = cell2mat(cellfun(@(x) mode(xlabel(x)),res.smallBinPruned,'UniformOutput',0))';
figure
plot(res.graph,'layout','force','NodeCData',nodecolor)
colormap jet