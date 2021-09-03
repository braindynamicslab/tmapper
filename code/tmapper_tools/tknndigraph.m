function g = tknndigraph(XorD,k,tidx,varargin)
%TKNNDIGRAPH construct a directed graph based on k-nearest neighbors which
%much include its temporal neighbors. Here we do so by simply set the
%distance between consecutive time points to zero before running knn.
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
(2-5-2020) add option to determine whether temporal link can be a spatial
link. parameter: timeExcludeSpace

%}
p = inputParser;
p.addParameter('reciprocal',true)% spatially reciprocal
p.addParameter('timeExcludeSpace',true)% whether temporal links are allow to be spatial links
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

% -- find indices for temporal links D_{i(t),i(t+1)}
t_wafter = circshift(tidx,-1,1) - 1 == tidx; % for which time points there exist a time point after
t_after_idx = circshift(diag(t_wafter),1,2);
if par.timeExcludeSpace
    D(t_after_idx) = 0;
end

% -- compute adjacency matrix
A = zeros(Nn,Nn);
[~,Ic]=sort(D,2);
I = sub2ind([Nn Nn], repmat((1:Nn)',1,k), Ic(:,1:k));
A(I(:))=1;

% -- exclude or retain temporal links as spatial links 
if par.timeExcludeSpace
    A_space = A.* (~t_after_idx); % remove temporal links
else
    A_space = A;
end

% -- enforce symmetry of spatial links
if par.reciprocal
    A_space = A_space & A_space';
else
    A_space = A_space | A_space';
end

% -- (re-)incoporate temporal links
A = t_after_idx | A_space; 

% -- convert to graph
g = digraph(A);
end

