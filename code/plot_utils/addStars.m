function [varargout] = addStars(ha,x,y,pval,varargin)
%ADDSTARS add asterisks to indicate significance level. For each greater
%level of significance, an marker is stacked on top. 
%   h = addStars(x,y,pval,...)
% input:
%   ha: axes handle.
%   x: x coordinates of potential markers. Vector of N elements.
%   y: y coordinates of potential markers (markers are aligned to the
%   bottom). Vector of N elements.
%   pval: p-values of each potential marker. Vector of N elements.
% output: 
%   h: m handles of markers, where m < M, given M levels of significance.
% parameters:
%   siglevels: significance levels, a vector of length M. (default: [0.05, 0.01,0.001])
%   sigsymbol: marker to label significance (default: '*')
%   symbolgap: y-distance between markers (default: range(ylim)/50)
%   symbolcolor: color of markers (default: 'k')
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 9-17-2020 ~

%}

p = inputParser;
p.addParameter('siglevels',[0.05, 0.01,0.001])
p.addParameter('sigsymbol','*')
p.addParameter('symbolgap',range(ylim)/50)
p.addParameter('symbolcolor','k')
p.parse(varargin{:})
par = p.Results;

if isempty(ha)
    ha = gca;
end
% -- check input dimensions
x=x(:);
y=y(:);
pval=pval(:);

if ~(length(x)==length(y) && length(y)==length(pval))
    error('first three arguments must be vectors of the same length')
end

% -- plotting
hold on
N_levels = length(par.siglevels);

for n=1:N_levels
    idx = pval<par.siglevels(n);
    if sum(idx)>0
        h(n)=plot(ha,x(idx),y(idx)+(n-1)*par.symbolgap,par.sigsymbol,...
            'color',par.symbolcolor,'linewidth',1);
    end
end

if nargout>0
    varargout{1} = h;
end

end

