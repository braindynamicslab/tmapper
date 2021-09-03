function [geod_n, nm] = normgeo(geod, nsize, varargin)
%NORMGEO normalize geodesic distance matrix, given that the measure of
%nodes may be different. 
%   [geod_n, nm] = normgeo(geod, nsize)
% input:
%   geod: N-by-N matrix, element(i,j) is the path length from node i or to
%   node j, in a network of N nodes.
%   nsize: N-by-1 vector, size of each node (e.g. how many sample points
%   included in each node). If this is not given by the user, it's a vector
%   of ones. 
% parameter:
%   excludeDiag: whether to excude the diagonal elements while computing
%   the normalizing factor or diameter of the network. Default = false. 
% output:
%   geod_n: N-by-N matrix, normalized geodesic distances by a factor
%               |sum_i,j [L(i,j)^2*m(i)*m(j)] |^(-1/2)
%           where L(i,j) is the path length from node i to node j, m(.) is
%           the probability measure of a node.
%   nm: a N-by-1 vector, probability measure of each node with sum(nm)=1.
%{
created by Mengsen Zhang, 9/24/2019
%}


p = inputParser;
p.addParameter('excludeDiag',0)% whether to exclude the diagonal of "geod" when calculating normalizing factor
p.parse(varargin{:})
par = p.Results;

N_nodes = length(geod);% number of nodes

if nargin < 2 || isempty(nsize)% if nsize not provided, assume = 1
    nsize = ones(N_nodes); 
end

% -- handle inf
geod(geod==Inf) = max(geod(geod<Inf));

% -- weight nodes and geodesics
nm = nsize/sum(nsize);% node measure
geod_n = (nm*nm').*geod.^2;% weight geodesic distance with node measures. 

% -- normalizing factor : the sum of weighted geodesics
if par.excludeDiag
    normfactor = sqrt(sum(geod_n(~logical(eye(N_nodes)))));
else
    normfactor = sqrt(sum(geod_n(:)));
end

% -- normalize
if normfactor~=0
    geod_n = geod/normfactor;
end
end

