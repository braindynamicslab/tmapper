% bifur_path.m
% ------------
% construct the "correct" path of bifurcation-induced phase transitions for
% simulated neural data.
%{
created by MZ, 8-21-2019

%}

addpath('../','tmapper_tools')

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
t = a.t(1:Nds:end); % subsampled time
tidx = (1:length(t))'; % index of t, s.t. could be used to en
x_label = state_regime(1:Nds:end); % regime of attractor
distMat = pdist2(x,x);
%% ===== create a priori transition map
% -- construct transition map and rank attractors by dwelltime
[dg,dwelltime,nodemembers] = symDyn2digraph(x_label);
[~, idx] = sort(dwelltime);
attractor_size_rank(idx,1) = (1:dg.numnodes);

% -- plot transition map
figure('position',[10 10 800 800])
plot(dg,'layout','force','NodeCData',str2double(dg.Nodes.Name),...
    'MarkerSize',10,'EdgeColor','k','NodeLabel','','linewidth',1)
colormap jet
axis square 
axis off
title({'theoretically predicted transition map', sprintf('for w_{EE}=%g, w_{EI}=%g',...
    a.w_EE,a.w_EI)})

% -- plot transition map with node size
figure('position',[810 10 800 800])
plot(dg,'layout','force','NodeCData',str2double(dg.Nodes.Name),...
    'MarkerSize',attractor_size_rank*0.5,'EdgeColor','k','NodeLabel',1:height(dg.Nodes),'linewidth',1)
colormap jet
axis square 
axis off
title({'theoretically predicted shape graph', sprintf('for w_{EE}=%g, w_{EI}=%g',...
    a.w_EE,a.w_EI)})
%% ===== create a priori recurrence plot
geod = TCMdistance(dg,nodemembers);
geod_undir = TCMdistance(digraph2graph(dg),nodemembers);
figure('position',[10 10 1000 400]);
% -- plot undirected geodesic
subplot(1,2,1)
imagesc(t,t,geod_undir)
cb = colorbar;
cb.Label.String = 'shortest path length';
colormap hot
axis square
xlabel('time (s)')
ylabel('time (s)')
title('geodesic on undirected graph')
% -- plot directed geodesic
subplot(1,2,2)
imagesc(t,t,geod)
cb = colorbar;
cb.Label.String = 'shortest path length';
colormap hot
axis square
xlabel('time (s)')
ylabel('time (s)')
title('geodesic on directed graph')
%% ===== transition detection 
figure('position',[10 10 1000 200]);
subplot(1,2,1)
plot(t,mean(geod_undir))
xlabel('time (s)')
ylabel('average path length')
title('undirected graph')
subplot(1,2,2)
plot(t,mean(geod)) % average sink
hold on
plot(t,mean(geod,2)) % average source
legend('as sink','as source')
xlabel('time (s)')
ylabel('average path length')
title('directed graph')
%% ===== cycle decomposition ===== %%
[cc,cl,cp,allcycles]=CycleCount2p(dg.adjacency);
T=CycleCluster(allcycles,0.5); 

%% ===== critical points: boundary between cycles ===== %%
[cluster_conn, cluster_conn_dir, ...
    clusters_nodes, clusters_boundary, clusters_interior,...
    clusters_crtpts, clusters_intcrtpts] ...
    = CycleClusterConn(dg, allcycles, T);

% -- plot transition map with labelled cycles
figure('position',[810 10 800 800])
ho=plot(dg,'layout','force','NodeCData',zeros(dg.numnodes,1),...
    'MarkerSize',attractor_size_rank*0.5,'EdgeColor','k','NodeLabel',1:height(dg.Nodes),'linewidth',1);
colormap jet
hbc=colorbar;
hbc.Label.String = 'cycle cluster index';
axis square 
axis off
title({'theoretically predicted shape graph', sprintf('for w_{EE}=%g, w_{EI}=%g',...
    a.w_EE,a.w_EI)})

for n=max(T)+1-unique(T')
    ho.NodeCData(clusters_nodes{n}) = n;
end
ho.NodeCData(unique(cell2mat(clusters_boundary))) = 0;


%% ===== connectivity between critical points ===== %%
allbd = unique(cell2mat(clusters_boundary));
allupath=Cycles2Paths(allcycles,allbd);
Tp=CycleCluster(allupath,0.5);
[pcluster_conn, pcluster_conn_dir, ...
    pclusters_nodes, pclusters_boundary, pclusters_interior,...
    pclusters_crtpts, pclusters_intcrtpts] ...
    = CycleClusterConn(dg, allupath, Tp);

disp(['# clusters = ' num2str(max(Tp))])

% -- plot transition map with labelled paths
figure('position',[810 10 800 800])
ho=plot(dg,'layout','force','NodeCData',zeros(dg.numnodes,1),...
    'MarkerSize',attractor_size_rank*0.5,'EdgeColor','k','NodeLabel',1:height(dg.Nodes),'linewidth',1);
colormap jet
hbp=colorbar;
hbp.Label.String = 'path cluster index';
axis square 
axis off
title({'theoretically predicted shape graph', sprintf('for w_{EE}=%g, w_{EI}=%g',...
    a.w_EE,a.w_EI)})

for n=max(Tp)+1-unique(Tp')
    ho.NodeCData(pclusters_nodes{n}) = n;
end
ho.NodeCData(unique(cell2mat(clusters_boundary))) = 0;

%% ===== all-in-one path decomposition ===== %%
[allbd,pclusters_nodes,pclusters_interior,pclusters_boundary,...
    pcluster_conn_dir,pclusters_intcrtpts, allupath, Tp] ...
    = CyclePathDecomp(dg);

% -- plot transition map with labelled paths
figure('position',[810 10 800 800])
ho=plot(dg,'layout','force','NodeCData',zeros(dg.numnodes,1),...
    'MarkerSize',attractor_size_rank*0.5,'EdgeColor','k','NodeLabel',1:height(dg.Nodes),'linewidth',1);
colormap jet
hbp=colorbar;
hbp.Label.String = 'path cluster index';
axis square 
axis off
title({'theoretically predicted shape graph', sprintf('for w_{EE}=%g, w_{EI}=%g',...
    a.w_EE,a.w_EI)})

for n=max(Tp)+1-unique(Tp')
    ho.NodeCData(pclusters_nodes{n}) = n;
end
ho.NodeCData(unique(cell2mat(clusters_boundary))) = 0;