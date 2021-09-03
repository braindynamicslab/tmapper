% tmapper_CME_task.m
% ------------------
%{
differentiate between tasks. 
created by MZ, 2/24/2020

%}

addpath(genpath('../'))
% addpath(genpath('../../GWTLBgpu/'))
addpath('plot_utils','tmapper_tools','other_utils')
%% ===== configure paths, import meta-data ===== %%
clear all
close all
clc

datapath = '../data/CME/'; % path of fMRI data
sumdatpath = [datapath 'behavior_cme.xlsx']; % path of summary data
savepath = '../results/CME/'; % path for saving results

sumdat = readtable(sumdatpath); % behavioral and meta-data
N_ss = height(sumdat); % number of subjects
%% ===== construct tmapper ===== %%
% --- construction parameters --- %
k=6;% default 5
d=3;% deftaul 2

tstart = tic;
for n = 1:N_ss
    % --- load data
    disp(sumdat.id(n));
    load(sprintf('../data/CME/%s_Shine_375.mat',sumdat.id{n})); 
    
    % --- construct tmapper
    tic
    sumdat.g{n} = tknndigraph (x,k,tidx,'reciprocal',true,'timeExcludeSpace',true); % tknn
    [sumdat.g_simp{n}, sumdat.members{n}, sumdat.nodesize{n}, sumdat.D_simp{n}]...
        = filtergraph(sumdat.g{n},d,'reciprocal',true); % simplified graph
    toc
end
telapsed = toc(tstart)

save([ savepath 'sumdat_tmapper_' par2filename({'k',k,'d',d}) '.mat'],'sumdat')
%% ===== modularity (classical) ===== %%
for n = 1:N_ss
    sumdat.Qmod(n) = full(Qasym(sumdat.g_simp{n}.adjacency, findnodelabel(sumdat.members{n},x_label)));
    sumdat.Qmod_knn(n) = full(Qasym(sumdat.g{n}.adjacency, x_label));
end

gscore = sumdat.Qmod;% graph score
bscore = sumdat.AveragePC;% behavioral score
figure;
scatter(gscore, bscore,'filled')
[RHO,P]=corrcoef(gscore, bscore);
disp(['rho=' num2str(RHO(1,2)) ', p=' num2str(P(2,1))])
