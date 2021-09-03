function cmap = CMEcmap(idx)
%CMECMAP Colormap for the CME data
%   
if nargin<1 || isempty(idx)
    idx=1:5;
end

cmap=lines(5);
cmap=[0 0 0; cmap([5,3,1,2],:); ];
cmap=cmap(idx(:),:);
end

