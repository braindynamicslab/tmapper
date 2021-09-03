function ax = xytickBlock(ax,pt,ptnames,varargin)
%XYLABELBLOCK change the tick labels of the x- and y-axis to name of the
%block that each block. 
%   ax = xytickBlock(ax,pt,ptnames,...)
% input:
%   ax: axes handle, where a N-by-N matrix is plotted as an image. 
%   pt: a N-by-1 vector of integer indices that partition rows/columns of
%   the matrix into M blocks. 
%   ptnames: a cell array of M elements (strings or numbers), which are
%   ticklabels to be displayed. 
% output:
%   ax: returen the axes handle.
%{
Created by: Mengsen Zhang <mengsenzhang@gmail.com> 3-3-2020
modifications:
%}

if isempty(ax)
    ax = gca;
end

if nargin<3 || isempty(ptnames)
    ptnames = unique(pt);
end

p = inputParser;
p.addParameter('axis','both'); % which axis to label: "x", "y", or "both"
p.addParameter('xlim',ax.XLim);
p.addParameter('ylim',ax.YLim);
p.parse(varargin{:});
par = p.Results;

% --- find ticks at the center of blocks
edges = find([true; diff(pt) ~= 0; true]);  % TRUE if values change
ticks = (edges(2:end)+edges(1:end-1))/2;
ticks = remapRange(ticks,1,length(pt),par.xlim(1), par.xlim(2));

% --- check names of blocks
if length(ticks) ~= length(ptnames)
    error('The number of partitions does not match the number of partition names.');
end

% --- label axes
if par.axis ~= 'y'
    set(ax,'xtick',ticks,'xticklabel',ptnames)
end
if par.axis ~= 'x'
    set(ax,'ytick',ticks,'yticklabel',ptnames)
end

end

