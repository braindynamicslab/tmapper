function [Div,Conv] = localDivConv(A,X,tidx)
%LOCALDIVCONV local divergences and convergence of the dynamics, based on
%nearest neighbors on a graph. 
%   Detailed explanation goes here

N = length(A);
A_t = mat2cell(A,ones(N,1),N);
% -- divergence at each node
Div = cell2mat(cellfun(@(tar) avgtargetdist(tar,X,tidx), A_t, 'UniformOutput',0));

if nargout>1
    Conv = localDivConv(A',X,tidx);
end

end
function d = avgtargetdist(tar,X,tidx)
    idx_wafter = tidx+1 == circshift(tidx,-1);
    d = mean(pdist(X(find(idx_wafter(logical(tar)))+1)));
end

