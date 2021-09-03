function [cycCount,cycLen,cycPath,allcycles] = CycleCount2p(A,varargin)
%CYCLECOUNT2P estimate the number of cycles of different length by finding
%the smallest cycle that passes through every two vertices in the graph or
%digraph. The function only counts unique cycles, and by default it only
%counts unique simple cycles. 
%   [cycCount,cycLen,cycPath,CYC] = CycleCount2p(A,...)
% input: 
%   A: N-by-N matrix, adjacency matrix for a simple directed graph. 
% output:
%   cycCount: M-by-1 vector, number of cycles per length. m-th element is
%   the number of cycles of m-th shortest length.
%   cycLen: M-by-1 vector, unique cycle length (smallest to largest). m-th
%   element is the m-th shortest cycle length.
%   cycPath: M-by-1 cell array. Each cell contains a Nm-by-Lm matrix, where
%   Nm is the number of cycles of length Lm. 
%   allcycles: C-by-1 cell array, each cell contains the path of one cycle.
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 3-23-2020 ~
modifications:
(9-1-2020) add output for all cycles in one cell array
%}


p = inputParser;
p.addParameter('simple',true); % only count simple cycles
p.parse(varargin{:});
par = p.Results;

G = digraph(A);
CYC = cell(nchoosek(G.numnodes,2),1);% preallocate storage for cycles

n = 0;
tic
for s = 1:G.numnodes
    for t = s+1:G.numnodes
        spf = shortestpath(G,s,t);% forward path
        spb = shortestpath(G,t,s);% backward path
        if ~isempty(spf) && ~isempty(spb)% path closed
            cycle = [spf, spb(2:end-1)];% cycle
            % -- check if cycle is simple
            if ~par.simple || (length(unique(cycle))==length(cycle))
                n = n + 1;
                % -- record cycle
                [~,rootIdx] = min(cycle);
                CYC{n} = circshift(cycle,1-rootIdx);% start cycle with smallest idx
            end
        end
    end
end
toc
CYC(n+1:end)=[];% remove empty cells

CYC_LEN = cellfun(@length,CYC);% length of each cycle
cycLen = unique(CYC_LEN);% unique cycle length

N_len = length(cycLen);% number of unique cycle length
cycPath = cell(N_len,1);% each cell contains all closed path for cycles of a particular length
cycCount = zeros(N_len,1);% number of cycles of a particular length

% -- find unique cycles for each cycle length
for lenIdx = 1:N_len
    cycPath{lenIdx}=unique(cell2mat(CYC(CYC_LEN == cycLen(lenIdx))),'rows');
    cycCount(lenIdx)=size(cycPath{lenIdx},1);
end


% -- put all cycles in a cell array
Nc=sum(cycCount);
allcycles=cell(Nc,1);
idx=1;
for ncl=1:N_len% cycle length
    for nc=1:cycCount(ncl)% cycle path
        allcycles{idx}=cycPath{ncl}(nc,:);
        idx=idx+1;
    end
end

end
