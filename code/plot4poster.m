% plot4poster.m
% -------------
% plotting for SfN2019

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

%% ===== plot neural activity
% -- subsampling
Nsp = 720;
tidx = 1:Nsp:a.Nt;
% -- rank activity
[~,orderByMeanSE] = sort(mean(a.X(:,1:a.N),1));

figure('position',[10 10 1200 750])
subplot(3,1,1)
plot(a.t(tidx),a.G(a.t(tidx)))
ylabel('global coupling')
ylim([1 5])

subplot(3,1,2)
hse= plot(a.t(tidx),a.X(tidx,orderByMeanSE));
colorlines(hse,[linspace(0,1,a.N)', zeros(a.N,2)])
ylabel('activity of excitatory populations')
ylim([0 1])

subplot(3,1,3)
hsi= plot(a.t(tidx),a.X(tidx,orderByMeanSE + a.N));
colorlines(hsi,[zeros(a.N,2), linspace(0,1,a.N)'])
ylabel('activity of inhibitory populations')
xlabel('time (s)')
ylim([0 1])
