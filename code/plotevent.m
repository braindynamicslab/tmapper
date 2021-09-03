function varargout=plotevent(x,cl,lw,varargin)
% function h=plotevent(x,color, linewidth)
%   use: 
%       input: a vector of time of the events, color of the plot
%       return: line handler
%   (created by MZ 2013-12-19)
%{
Modifications: 
2014-07-02: allow multiple color
2015-07-19: return pointers
2016-08-30: length(cl) => size(cl,1)
%}
%%%%%
%--- check 
    if nargin<3,lw=2;end
    if nargin<2,cl='r';end
%     if isempty(cl),cl='r';end
%--- plot
    hold on;
    n=length(x); % number of events
%     h=[];
    for i=1:length(x)
        h(i)=plot(x(i)*[1 1],ylim,'--','color',cl(min(i,size(cl,1)),:),'linewidth',lw,varargin{:});
    end
% --- output
if nargout>0
    varargout{1}=h;
end

end