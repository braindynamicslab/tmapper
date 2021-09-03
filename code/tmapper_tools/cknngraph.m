function g = cknngraph(XorD,k,delta,varargin)
%CKNNGRAPH construct graph based on continuous k-nearest neighbors, taking
%ideas from Berry & Sauer (2019).
% g = cknngraph(XorD,k,delta)   
% -----------------------------------------------------------------------
% Create graph by connecting points x,y if
%       D(x,y) < delta * (sqrt(D(x,x_k))*D(y,y_k))
% where D(~,~) is the distance, and ~_k is the k-nearest neighbor of ~.
% -----------------------------------------------------------------------
% input:
%   XorD: a N-by-d matrix (X) of the coordinates of N points in d-sim
%       space, or a N-by-N distance matrix (D). If the input is X, D is
%       computed from X based on Euclidean distances. 
%   k: k-nearest neighbors, used to normalize distance based on local point
%   density, i.e. D_norm = D(x,y) ./ (sqrt(D(x,x_k))*D(y,y_k))
%   delta: threshold for linking two nodes, i.e. D_norm < delta.
% parameters:
%   average: whether x_k is the average distance over the first k-nearest
%   neighbors (=true) or just the k-th nearest neighbor (=false). Default
%   is false.
% output:
%   g: matlab graph object (unweighted, undirected).
%{
created by MZ, 8-20-2019

%}


% -- check parameters
p=inputParser;
p.addParameter('average',false) % whether average the first knns, or just use the k-th nn (default)
p.parse(varargin{:});
par=p.Results;


% -- check input and obtain distance matrix D
[nr,nc]=size(XorD);
if nr~=nc || any(any(XorD~=XorD'))
    D = pdist2(XorD,XorD);
else
    D = XorD;
end
Nn = length(D); % number of nodes

% -- find nearest neighbors and compute distance
D_sorted=sort(D,2);
if par.average
    Dk = mean(D_sorted(:,2:1+k),2); % note: first col is dist to self, (1+k)th col is dist to knn
else
    Dk = D_sorted(:,1+k);
end

% -- normalize D by local density
D_norm = D./sqrt(Dk*Dk');

% -- construct graph from adjacency matrix
A = D_norm < delta;
A(logical(eye(Nn,Nn)))=0; % set diagonal to zero
g = graph(A);

end