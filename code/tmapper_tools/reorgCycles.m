function allcycles = reorgCycles(cycPath)
% reorgCycles unpack `cycPath` output of `CycleCount2p` to a cell array,
% where each cell only contains one path. 
%   allcycles = reorgCycles(cycPath)

N_len = length(cycPath);
cycCount = cellfun(@(x) size(x,1), cycPath);

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