function nodepath = CycleCutter(cyc,nodeName)
%CYCLECUTTER cut a cycle into multiple paths at given nodes.
%   nodepath = CycleCutter(cyc,nodeName)
% input:
%   cyc: 1-by-N array of nodeNames (usually integers). 
%   nodeNames: vector of any size, containing names of nodes that are the
%   cutting points of the cycle. Names not contained in the cycle is
%   ignored.
% output:
%   nodepath: M-by-1 cell array. Each cell contains a path from one cutting
%   point to the next. 
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 9-2-2020 ~

%}

nodeIdx = find(ismember(cyc,nodeName));
if isempty(nodeIdx)
    nodepath=cyc;
    return
end

N_path = length(nodeIdx);
nodepath = cell(N_path,1);

if N_path>1% when there is more than one cutting point
    for n=1:N_path-1
        nodepath{n}=cyc(nodeIdx(n):nodeIdx(n+1));
    end
else% when there's only one cutting point
    n=0;
end

nodepath{n+1}=[cyc(nodeIdx(n+1):end), cyc(1:nodeIdx(1))];
end

