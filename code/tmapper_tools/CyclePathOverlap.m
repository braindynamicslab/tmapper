function [CO,grpnames] = CyclePathOverlap(c,varargin)
%CYCLEPATHOVERLAP calculate the overlap between cycles or paths.
%   [CO] = CyclePathOverlap(c,'type','edge','cycle',true)
% input:
%   c: a N-by-1 cell array of cycles or paths. The n-th cell contains a
%   vector of M_n elements, which are indices of nodes along a path.
% output:
%   CO: a N-by-N matrix of percent overlap between any two cycles.
%   C0(ii,jj) is the the percent overlap between the ii-th and jj-th cycle.
%   or a P-by-P matrix for P groups of cycles. 
% parameters:
%   cycle: whether the paths in `c` are cycles (default true). 
%   type: calculate overlap between edges (`edge`) or nodes (`node`).
%   Default `edge`.
%   grpvar: a grouping variable of the paths, a N-by-1 vector of integers,
%   assigning a group index for each cycle. By default, each cycle is
%   its own group.
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 10-13-2020 ~

%}

p=inputParser;
p.addParameter('cycle',true) % whether input paths are cycles
p.addParameter('type','edge') % edge overlap or node overlap
p.addParameter('grpvar',[]) % a grouping variable for paths/cycles
p.parse(varargin{:})
par = p.Results;

% -- adjust path representation
switch par.type
    case 'edge'
        c = cellfun(@(x) [x(:), circshift(x(:),-1,1)], c,'uniformoutput',0);
        if ~par.cycle % remove link back to the first node if not a cycle
            c = cellfun(@(x) x(1:end-1,:), c, 'uniformoutput',0);
        end
    case 'node'
        c = cellfun(@(x) x(:), c,'uniformoutput',0);
end

% -- grouping cycles
Nc = length(c);
if isempty(par.grpvar)
    par.grpvar = (1:Nc);
end
grpnames = unique(par.grpvar);
Ngrp = length(grpnames);

% -- compute overlap
CO = nan(Ngrp,Ngrp);

for ii=1:Ngrp
    for jj=1:ii-1
        tmpii = unique(vertcat(c{par.grpvar==grpnames(ii)}),'rows');
        tmpjj = unique(vertcat(c{par.grpvar==grpnames(jj)}),'rows');
        CO(ii,jj) = size(intersect(tmpii,tmpjj,'rows'),1)...
            /size(union(tmpii,tmpjj,'rows'),1);
        CO(jj,ii)=CO(ii,jj);
    end
end

CO(logical(eye(size(CO))))=1;
end

