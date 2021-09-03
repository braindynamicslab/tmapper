function [hm,hse] = plotSpectraSEM(x,y,se,varargin)
%PLOTSPECTRASEM plot spectra with standard errors
%   [hm,hse] = plotSpectraSEM(x,y,se,'cmap',...,'alpha',...)
% input: 
%   x: N-by-Nd matrix, N is the number of lines to be plotted, Nd is the
%   number of points in each line. The (n,nd)-th element is the
%   x-coordinate of the nd-th point in the n-th line.
%   y: y-coordinates, N-by-Nd matrix, similar to "x".
%   se: standard error or some kind of interval around y. A N-by-Nd matrix,
%   similar to "x" and "y".
% output:
%   hm: N-by-1 graphical object array, handles of the lines defined by x,
%   y.
%   hse: N-by-1 graphical object array, handles of the error band defined
%   by se. 
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 03/27/2020 ~
modifications:

%}

% -- check inputs
[Nx,Ndx]=size(x);
[Ny,Ndy]=size(y);
[Nse,Ndse]=size(se);
if Nx==1 && Ny>1
    x = repmat(x,Ny,1); % match the number of rows in x to that of y
    Nx = Ny;
end

if (Nx~=Ny) || (Ndx~=Ndy) || (Nx~=Nse) || (Ndx~=Ndse)
    error('Input matrix dimensions mismatch.')
end

% -- parse parameters
p = inputParser;
p.addParameter('cmap',lines(Nx)) % colormap, default "lines"
p.addParameter('alpha',0.3) % transparency of the error band
p.parse(varargin{:})
par = p.Results;

% -- graphical handles
[hm,hse] = deal(gobjects(Nx,1));
% -- plot means and the error bands
hold on
for n = 1:Nx
    hse(n) = plotColorBand(x(n,:),y(n,:)-se(n,:),y(n,:)+se(n,:),...
            par.cmap(n,:),par.alpha);% standard error
    hm(n) = plot(x(n,:), y(n,:),'color',par.cmap(n,:)); % mean
end


end

function h = plotColorBand(x, lowerbd, upperbd, cl, al)
    h = fill([x fliplr(x)]',...
        [lowerbd,fliplr(upperbd)]',...
        cl,'FaceAlpha',al,'EdgeAlpha',al);
end