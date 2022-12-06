%% Generate simulated data
% Mengsen Zhang (MZ, 10/27/2022

%% ===== simulate time series =====
clear all
close all
clc

% -- keep simulation parameters of interest
sigma = 0.01;
dt = 0.0005;
% -- create data output folder
outfolder = "./network_" + par2filename({'sigma',sigma,'dt',dt});
mkdir(outfolder);
% -- number of trials to simulate
N_trials = 50;
SEEDS = 1:2:2*N_trials;
% -- tmapper parameters
SI = 0.72; % (s) sampling interval
k_se = 16;
d_se = 10;
k_bold = 14;
d_bold = 10;

% -- connectome type
Cname = 'HCP11avg';% name of structural connectivity
% -- define trial duration
T = 1200; % (s)

% -- model design
wee = 2.8;
wei = 1;
G_range = [1.7 5];
G = @(t) cosTentTrain(t,'yrange',G_range,'period',T/2);

% -- define connectome
Nnodes = 66;
switch Cname
    case "HCP11avg"
        load ../data/connectome/HCP11avg_66.mat
    case "HCP100avg"
        load ../data/connectome100/HCP100avg_66.mat
end
C = connAvg66;
indg = mean(C,2); % in-degree

% -- load attractor labels
load(['../data/labelAttractors/XssCluster_' par2filename({'w_EE',wee,'w_EI',wei}) '.mat'],'fps','clusterBySE','clustermean')

% -- create model set parameters
a = WWWC(C,G);
a.T = T;
a.dt = dt;
a.w_EE = wee;
a.w_IE = wee;
a.w_EI = wei;
a.sigma = sigma;

% -- choose initial condition
load([ '../data/' Cname '/Gsweep_Xss0-05_' par2filename({'w_EE',wee,'w_EI',wei}) '.mat'],'fixpts','Gopts');
fixpts0 = fixpts{knnsearch(Gopts',a.G(0))}; % options of steady states
[~,fixptsRank] = sort(mean(fixpts0(:,1:a.N),2));
X0=fixpts0(fixptsRank(2),:)'; % choose the second smallest fixed points (the smallest is zero)

%  -- generate a trial (loop over N_trials if desired)

n = 1;
disp("trial#" + n)
% -- simulate gating variables
disp("simulating SE...")
seed_tr = SEEDS(n);
rng(seed_tr);
tic
a = a.HeunSolver(X0);
toc
% -- simulate BOLD
disp("simulating BOLD...")
b = BallWindBOLD(a.X(:,1:a.N));
b.T = a.T;
b.dt = a.dt;
b.t_transient = 30;
tic
b = b.EularSolver;
toc

% -- subsampling SE
Nds = round(SI/a.dt); % step of down sampling
spidx = 1:Nds:a.Nt; % sample indices

% -- label states by attractor index
disp("attractor assignment...")
states = [a.G(a.t(spidx)) a.X(spidx,:)];
tic
[idx,D] = knnsearch(fps,states);
toc
x_label = clusterBySE(idx);

% -- create ground truth transition map
disp("constructing graphs...")
tic
[dg0,dwelltime0,nodemembers0] = symDyn2digraph(x_label);
geod0 = distances(dg0);
nsize0 = nodesize(nodemembers0);
[geod0_n,mNode0] = normgeo(geod0,nsize0);
% -- create ground truth recurrence plot
tcm0 = TCMdistance(dg0,nodemembers0); 
tcm0_n = normtcm(tcm0,'normtype','max');

% ----- create network from SE ----- %
x_se = a.X(spidx,1:a.N); % subsampled data
figure;plot(x_se);

N_t = length(spidx);
tidx = (1:N_t)';
g_se = tknndigraph (x_se,k_se,tidx,'reciprocal',true,'timeExcludeSpace', true);
[g_simp_se, members_se, nsize_se] = filtergraph(g_se,d_se,'reciprocal',true);
% -- compute geodesic
geod_se = distances(g_simp_se,'Method','unweighted');
[geod_n_se, mNode_se] = normgeo(geod_se,nsize_se);
% -- compute recurrence plot
tcm_se = TCMdistance(g_simp_se, members_se);
tcm_n_se = normtcm(tcm_se,'normtype','max');

% ----- create network from BOLD ----- %
x_bold = b.B(spidx,:); % subsampled data
figure; plot(x_bold);

g_bold = tknndigraph (x_bold,k_bold,tidx,'reciprocal',true,'timeExcludeSpace', true);
[g_simp_bold, members_bold, nsize_bold] = filtergraph(g_bold,d_bold,'reciprocal',true);
% -- compute geodesic
geod_bold = distances(g_simp_bold,'Method','unweighted');
[geod_n_bold, mNode_bold] = normgeo(geod_bold,nsize_bold);
% -- compute recurrence plot
tcm_bold= TCMdistance(g_simp_bold, members_bold);
tcm_n_bold = normtcm(tcm_bold,'normtype','max');
