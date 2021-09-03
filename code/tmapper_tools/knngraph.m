function g = knngraph(XorD,k,varargin)
%KNNGRAPH construct graph based on k-nearest neighbors. Here we take a
%naive approach by define each point as a node and adding a link between
%two nodes if one is the other's k-nearest neighbors.
%   g = knngraph(XorD,k)
%   input: 
%       XorD: a N-by-d matrix (X) of the coordinates of N points in d-sim
%       space, or a N-by-N distance matrix (D). If the input is X, D is
%       computed from X based on Euclidean distances. 
%       k: # nearest neighbors
%   output:
%       g: matlab graph object (unweighted, undirected).
%{
by MZ, 8/16/2019
modifications:
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

% -- compute adjacency matrix
A = zeros(Nn,Nn);
[~,Ic]=sort(D,2);
I = sub2ind([Nn Nn], repmat((1:Nn)',1,k), Ic(:,2:1+k));
A(I(:))=1;
if par.reciprocal
    A = A & A';
else
    A = A | A';
end

% -- convert to graph
g = graph(A);
end

