function [dg,dwelltime,nodemembers] = symDyn2digraph(symDyn)
%SYMDYN2DIGRAPH construct a digraph (matlab object for directed graph)
%based on from a single time series of symbolic dynamics.
%   [dg, dwelltime, nodemembers] = symDyn2digraph(symDyn)
% input:
%   symDyn: a vector of N integers. Each integer is the name of a unique
%   state of the system. This vector is purely nominal, i.e. the numeric
%   difference between any two integers are irrelevant. 
% output:
%   dg: matlab object digraph. Each node is a state in symDyn, named by a
%   string converted from the integer. Each edge indicates that a
%   transition exists between two states.
%   dwelltime: how many time points belong to each state.
%   nodemembers: a N-by-1 cell array of a digraph of N nodes, the n-th cell
%   includes a vector of indices of the time points associated with each
%   state. 
%{
created by MZ, 8-21-2019
modifications:
(8-23-2019) adding option to return temporal members of each node, i.e.
"nodemembers".

%}

tidx = (1:length(symDyn))';
% -- find source and target of each transition 
tidx_trans = tidx(diff(int32(symDyn))~=0);
state_s = symDyn(tidx_trans);
state_t = symDyn(tidx_trans+1);

% -- extract unique transitions (define edges)
trans = unique([state_s state_t],'row');

% -- define states (nodes)
[state_names,~,state_idx] = unique(symDyn);

% -- calculate # time points in state (dwell time)
dwelltime = accumarray(state_idx,1);

% -- create digraph
dg = digraph(cellstr(num2str(trans(:,1))),cellstr(num2str(trans(:,2))),[],cellstr(num2str(state_names)));

% -- if asked, find members of each node (time points where )
if nargout>2
    nodemembers = cell(dg.numnodes,1);
    for n = 1:dg.numnodes
        nodemembers{n} = tidx(state_idx == n);
    end
end

end

