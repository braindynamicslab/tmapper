% visualtrial_knn.m
% -----------------
% visualize dynamics with various types of knn
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
% x = a.X(1:Nds:end,1:a.N); % subsampled data
x = B.B(1:Nds:end,:); % BOLD
t = a.t(1:Nds:end); % subsampled time
tidx = (1:length(t))'; % index of t, s.t. could be used to en
x_label = state_regime(1:Nds:end); % regime of attractor
distMat = pdist2(x,x);
%% ===== or load human data
ss = 99;
load(['../data/CME/SBJ' num2str(ss,'%02g') '_Shine_375.mat'])
%% ===== construct nearest neighbor graph
% [~, mapping] = bdl_isomap(x,2,num_k);
k = 17;700;5;
delta = [];2;
tic
g = tknndigraph (x,k,tidx,'reciprocal',true);
% g = tknngraph (x,k,tidx,'reciprocal',true);
% g = knngraph(x,k,'reciprocal',false);
% g = cknngraph(x,k,delta,'average',0);
toc

% -- plotting
[a1, a2] = plotgraphtcm(g,x_label,t,[],'nodesizerange',5);
title(a1,['k=' num2str(k)])
% title(sprintf('$w_{EE}=%g, w_{EI}=%g, k=%g, \\delta=%g$',...
%     a.w_EE,a.w_EI,k,delta),'Interpreter','latex')
%% ===== transition detection 
figure('position',[10 10 500 200]);
plot(t,mean(D_geo)) % average sink
hold on
plot(t,mean(D_geo,2)) % average source
legend('as sink','as source')
xlabel('time (s)')
ylabel('average path length')
title('directed graph')
%% ===== compare to traditional methods
figure('position',[10,10,1000,400]);
% -- plot regular recurrence plot
subplot(1,2,1)
imagesc(t,t,pdist2(x,x))
colormap hot
cb=colorbar;
axis square
xlabel('time (s)')
ylabel('time (s)')
xlabel(cb,'distance between S_E')
title('Recurrence of S_E')
% -- compute dynamic (windowed) functional connectivity
winsize = 30;
lag = 1;
[dfc,dfc_fz] = dynFC(x,winsize,lag);
dfc = cellfun(@(x) belowdiag(x)',dfc,'uniformoutput',false);
dfc = cell2mat(dfc);
% -- plot regular recurrence plot
subplot(1,2,2)
imagesc(t,t,pdist2(dfc,dfc))
colormap hot
cb=colorbar;
axis square
xlabel('time (s)')
ylabel('time (s)')
xlabel(cb,'distance between dFC')
title('Recurrence of dFC')

% -- plot average distance
figure('position',[10 10 800 600]);
% plot control variable
hsp1=subplot(3,1,1);
scatter(t,a.G(t),1,x_label)
ylabel('G')
colormap(gca,'jet')
% plot regime
hsp2=subplot(3,1,2);
plot(t,x_label)
ylabel('attractor index')
% plot distance
hsp3=subplot(3,1,3);
plot(t,normrange(mean(D_geo))) % average sink
hold on;
plot(t,normrange(mean(pdist2(x,x))))
plot(t(1:lag:size(dfc,1)),normrange(mean(pdist2(dfc,dfc))))
legend('Geodesic sink','S_E','dFC','location','northeastoutside')
xlabel('time (s)')
ylabel('average distance')
[hsp1.Position(3), hsp2.Position(3)]= deal(hsp3.Position(3));

% -- recurrence of G
figure
imagesc(pdist2(a.G(t),a.G(t)))
colormap hot
cb=colorbar;
axis square
xlabel('time (s)')
ylabel('time (s)')
xlabel(cb,'distance between G')
title('Recurrence of G')
%% ===== graph laplacian ===== %%
