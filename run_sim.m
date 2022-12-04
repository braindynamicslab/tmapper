%% Run Temporal Mapper on simulated data

% Load code and data
addpath(genpath('code/'));
load('data/sim.mat');
tidx = (1:length(data.t))'; % indexing variable

% Set parameters
k = 14; % nearest neighbor parameter
d = 2;

% Compute temporal mapper
g = tknndigraph (data.x,k,tidx,'reciprocal',true);
[g_simp, members, ~] = filtergraph(g,d,'reciprocal',true); % Perform simplification
[a1, a2] = plotgraphtcm(g_simp,data.x_label,data.t,members,'nodesizerange',[1,25]);

%% Comparison to other methods

figure('position',[10,10,1000,400]);
% Plot regular recurrence plot
subplot(1,2,1)
recur_reg = pdist2(data.x,data.x);
imagesc(t,t,recur_reg)
colormap hot
cb=colorbar;
axis square
xlabel('time (s)')
ylabel('time (s)')
xlabel(cb,'distance between S_E')
title('Recurrence of S_E')

% Compute dynamic (windowed) functional connectivity
winsize = 30; %default 30
lag = 1;
[~,dfc_fz] = dynFC(data.x,winsize,lag);
% -- plot recurrence plot of dFC
recur_dfc = pdist2(dfc_fz,dfc_fz);
subplot(1,2,2)
imagesc(t,t,recur_dfc)
colormap hot
cb=colorbar;
axis square
xlabel('time (s)')
ylabel('time (s)')
xlabel(cb,'distance between dFC')
title('Recurrence of dFC')

% Recurrence of control parameter G
figure
imagesc(pdist2(data.G,data.G))
colormap hot
cb=colorbar;
axis square
xlabel('time (s)')
ylabel('time (s)')
xlabel(cb,'distance between G')
title('Recurrence of G')
