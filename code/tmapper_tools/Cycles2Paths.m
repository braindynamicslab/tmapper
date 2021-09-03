function allupath = Cycles2Paths(allcycles,cutpts)
%CYCLES2PATHS cut a set of cycles to obtain a set of path connecting a set
%of cutting points (nodes). 
%   allupath = Cycles2Paths(allcycles,cutpts)
% input:
%   allcycles: N-by-1 cell array, of N cycle paths (row vectors).
%   cutpts: a vector of M elements, integer indices of nodes at which each
%   cycle should be cut -- cutting points. 
% output:
%   allupath: M-by-1 cell array, of M paths linking the cutting points. 
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 9-3-2020 ~
%}

% -- cut all the cycles at the cutting points
allpath=cellfun(@(x) CycleCutter(x,cutpts),allcycles,'UniformOutput',0);
allpath=vertcat(allpath{:});

% -- sort path by length and find unique path for each length
pathlen=cellfun(@length,allpath);
upathlen=unique(pathlen);
N_len = length(upathlen);
sortedpath=cell(N_len,1);
for n=1:N_len
    sortedpath{n}=unique(cell2mat(allpath(pathlen==upathlen(n))),'row');
end

% -- make a single list of unique paths
N_path = sum(cellfun(@(x) size(x,1),sortedpath));
allupath = cell(N_path,1);
idx = 1;
for n=1:N_len    
    Nc=size(sortedpath{n},1);
    for nc=1:Nc
        allupath{idx}=sortedpath{n}(nc,:);
        idx = idx + 1;
    end
end

end

