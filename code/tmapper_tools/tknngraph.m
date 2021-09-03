function g = tknngraph(XorD,k,tidx,varargin)
%TKNNGRAPH construct graph based on k-nearest neighbors which much include
%its temporal neighbors. Here we do so by simply set the distance between
%consecutive time points to zero before running knn.
%   g = tknngraph(XorD,k,tidx)
% input:
%   XorD: a N-by-d matrix (X) of the coordinates of N points in d-sim
%       space, or a N-by-N distance matrix (D). If the input is X, D is
%       computed from X based on Euclidean distances. 
%   k: # nearest neighbors
%   tidx: a vector of N integers, two points x, y are considered temporal
%   neighbors iff tidx[x]+1 = tidx[y] or tidx[x]-1 = tidx[y].
% output:
%   g: matlab graph object (unweighted, undirected).
%{
created by MZ, 8-16-2019
modifcations:
(8-20-2019) add option to not enforce reciprocity.


%}
p = inputParser;
p.addParameter('reciprocal',true)
p.parse(varargin{:})
par = p.Results;

% -- check input and obtain distance matrix D
[nr,nc]=size(XorD);
if nr~=nc || any(any(XorD~=XorD'))
    D = pdist2(XorD,XorD);
else
    D = XorD;
end
Nn = length(D); % number of nodes

D(logical(eye(Nn))) = Inf; % exclude self-loops

% -- make sure temporal neighbors are neareat neighbors (d=0)
t_wbefore = circshift(tidx,1,1) + 1 == tidx; % for which time points there exist a time point before
t_wafter = circshift(tidx,-1,1) - 1 == tidx; % for which time points there exist a time point after
D(circshift(diag(t_wbefore),-1,2)) = 0;
D(circshift(diag(t_wafter),1,2)) = 0;

% -- compute adjacency matrix
A = zeros(Nn,Nn);
[~,Ic]=sort(D,2);
I = sub2ind([Nn Nn], repmat((1:Nn)',1,k), Ic(:,1:k));
A(I(:))=1;
if par.reciprocal
    A = A & A';
else
    A = A | A';
end

% -- convert to graph
g = graph(A);
end

