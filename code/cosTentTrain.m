function y = cosTentTrain(t,varargin)
%COSTENTTRAIN generate a train of tent functions punctuated by epochs of
%plateaus. Each period of the train consists of a ramp-up and a ramp-down
%epoch:
%  _/¯\
%   where / and \ are of the form [1-cos(2pi*wt)]/2 and [1+cos(2pi*wt)]/2 respectively, while
%   _ and ¯ are the min and max of the output for the duration of 1/(2w). 
% if t exceeds the period w, it will repeat itself:
%  _/¯\_/¯\_/¯\...
% usage:
%   y = cosTentTrain(t,...)
% input:
%   t: time, can be either a scalar or a vector.
% parameters:
%   yrange: [y_min, y_max], if only a scalar is given, it is equivalent to
%   [0, yrange]. Default is [0,1];
%   period: period of the tent train. Default is 1.
% output:
%   y: the tent train. It is of the same size as t. 
%{
created by MZ, 8-8-2019
%}

p = inputParser;
p.addParameter('yrange',[0 1]);
p.addParameter('period', 1);
p.parse(varargin{:})

par = p.Results;

% -- construct function
t = rem(t,par.period);
y = nan(size(t));

t_e1_Idx = t < par.period/6;% low plateau
t_e2_Idx = par.period/6 <= t & t < par.period/2;% ramp-up
t_e3_Idx = par.period/2 <= t & t < 2/3*par.period;% high-plateau
t_e4_Idx = 2/3*par.period <= t;% ramp-down

y(t_e1_Idx) = 0;
y(t_e2_Idx) = (1-cos(3*pi*t(t_e2_Idx)/par.period-pi/2))/2; 
y(t_e3_Idx) = 1; 
y(t_e4_Idx) = (1+cos(3*pi*t(t_e4_Idx)/par.period-2*pi))/2; 

% -- scaling the function to a range
switch numel(par.yrange)
    case 1
        y = y*par.yrange;
    case 2
        y = y*range(par.yrange);
        y = y + par.yrange(1);
end

end

