function H = entropyDiscrete(X)
%ENTROPYDISCRETE entropy of discrete variables. 
%   H = entropyDiscrete(X)
% input:
%   X: a N-by-M matrix of M discrete variables, each column contains N
%   observations of a discrete variable.
% output:
%   Y: a 1-by-M vector, containing entropy (bits) for each column of X.
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 10-15-2020 ~

%}

[~,ncol] = size(X);

if ncol>1
    H = cellfun(@entropyDiscrete,num2cell(X,1));
else
    edges = unique(X);
    edges = [edges(:); edges(end)+1];

    % -- calculate distribution
    P = histcounts(X,edges,'Normalization','probability');

    % -- calculate entropy
    h = P .* log2(P);
    h(isnan(h)) = 0;
    H = - sum(h);
end
end

