function CI = confint(P,a)
%CONFINT confidence interval (2-tails) given a distribution, and
%significance level
%   CI = confint(P,a)
% input:
%   P: N-by-M matrix, N is the number of observations per data vector, M is
%   the number of data vectors.
%   a: significance level, default 0.05
% output:
%   CI: confidence intervals, 2-by-M matrix. 
%{
created by MZ, 9/25/2019
%}

if nargin<2 || isempty(a)
    a = 0.05;
end

CI = prctile(P, [a*100/2 (1-a/2)*100]);
end

