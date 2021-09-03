% tmapper_vsothers.m
% ---------------
%{
Compare temporal mapper to other methods (regular recurrence, dFC)

%}
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
return
% -- subsample data
Nds = 720; % step of down sampling
% x = a.X(1:Nds:end,1:a.N); % subsampled data (S_E)
x = B.B(1:Nds:end,:); % subsampled BOLD
t = a.t(1:Nds:end); % subsampled time
Nt = length(t);
tidx = (1:Nt)'; % index of t, s.t. could be used to en
x_label = state_regime(1:Nds:end); % regime of attractor
distMat = pdist2(x,x);
%% ===== construct simplify the graph
k = 14;
d = 2;
g = tknndigraph (x,k,tidx,'reciprocal',true);
tic
[g_simp, members, ~] = filtergraph(g,d,'reciprocal',true);
toc
% -- plotting
[a1, a2] = plotgraphtcm(g_simp,x_label,t,members,'nodesizerange',[1,25]);
title(a1,['k=' num2str(k) ', d=' num2str(d)]);
% -- compute TCM
tcm_knn = normtcm(distances(g));
tcm = normtcm(TCMdistance(g_simp,members));
%% ===== compare to traditional methods: recurrence
figure('position',[10,10,1000,400]);
% -- plot regular recurrence plot
subplot(1,2,1)
recur_reg = pdist2(x,x);
imagesc(t,t,recur_reg)
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
[~,dfc_fz] = dynFC(x,winsize,lag);
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

%% ===== compare to traditional methods: average
figure('position',[10 10 1200 600]);
% plot control variable
hsp1=subplot(3,1,1);
scatter(t,a.G(t),1,x_label)
ylabel('G')
colormap(gca,'jet')
% plot regime
hsp2=subplot(3,1,2);
scatter(repmat(t,size(x,2),1),x(:),[],repmat(x_label,size(x,2),1),'.')
ylabel('BOLD')
ylim([0 0.05])
colormap jet
% plot distance
hsp3=subplot(3,1,3);
% stem(t(1:end-1),diff(x_label)~=0,'k')
hold on;
plot(t,normrange(mean(tcm,2))) % average source
plot(t,normrange(mean(recur_reg)))
plot(winsize*Nds*a.dt + t(1:lag:size(dfc_fz,1)),normrange(mean(recur_dfc)))
xlim([min(t) max(t)])
legend('TCM_o','S_E','dFC','location','northeastoutside')
xlabel('time (s)')
ylabel('average distance')
[hsp1.Position(3), hsp2.Position(3)]= deal(hsp3.Position(3));
%% ===== compare to traditional methods: change of recur
d_tcm = sqrt(sum([diff(tcm,[],1),diff(tcm,[],2)'].^2,2));
d_recur_reg = sqrt(sum(diff(recur_reg).^2,2));
d_recur_dfc = sqrt(sum(diff(recur_dfc).^2,2));

figure;
hold on
stem(t(1:end-1),diff(x_label)~=0,'k')
plot(t(1:lag:size(dfc,1)-1),normrange(d_recur_dfc))
plot(t(1:end-1),normrange(d_recur_reg))
plot(t(1:end-1),normrange(d_tcm))
