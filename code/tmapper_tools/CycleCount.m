function [Primes,elapsedTime] = CycleCount(A,L0)
% CycleCount: counts all simple cycles of length up to L0 included on a
%             graph whose adjacency matrix is A.
%
% USAGE: [Primes,elapsedTime]=CycleCount(A,L0);
%
% INPUTS: A: full or sparse (faster) adjacency matrix of the graph G 
%         L0: maximum length of the simple cycles to be counted
%
% OUTPUTS: Primes: a list, Primes(i) is the number of primes of length i>=1
%                  up to L0.
%          elapsedTime: time taken in sec. by the algorithm
%
% Authors: P.-L. Giscard, N. Kriege, R. Wilson, Septembre 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if ~issparse(A)            % Test whether the matrix is sparse
   A = sparse(A);
end

Primes = zeros(1,L0);      % Initialisation of the list of primes
Primes(1) = trace(A);      % Counts the self-loops
A = A - diag(diag(A));     % Gets rid of the self-loops   
A = CleanThis(A);          % "Cleans" the matrix: removes isolated vertices, sinks and sources

if  issymmetric(A)         % Test whether the graph is directed
    Aundir = +(A~=0);      % Unweighted skeleton of the graph
    directed = false;      % Graph is not directed
else  
    Aundir = +((A~=0)|(A~=0)'); % Undirected, unweighted skeleton of the graph
    directed = true;       % Graph is directed
end

Size = size(A,1);          % Number of vertices    
if L0 > Size               % Checks that maximum length does not exceed 
   L0 = Size;              % the total number of vertices    
end

%%%%%%%%%%%%%%%%%%%%% SUBGRAPHS RECURSIVE CONSTRUCTION %%%%%%%%%%%%%%%%%%%%  
% Indicator vector of vertices that one may consider adding to a subgraph, initially all
AllowedVert = true(1,Size); 
    
    tic;
    for i = 1:(Size-1)                 % Loops over the vertices of the graph, last vertex can be disregard as at this point all others are forbidden
        AllowedVert(i) = false;        % Vertex visited is now forbidden
        Neighbourhood = zeros(1,Size); % Initialize indicator vector for neighbourhood
        Neighbourhood(i) = 1;          % current subgraph
        Neighbourhood = Neighbourhood + Aundir(i,:); % vertices reachable via one edge

        Primes = RecursiveSubgraphs(A,Aundir,directed,L0,i,AllowedVert,Primes,Neighbourhood); % Forms larger subgraphs containing the current one
    end
    elapsedTime = toc;

end



function Primes = RecursiveSubgraphs(A,Aundir,directed,L0,Subgraph,AllowedVert,Primes,Neighbourhood)
% RecursiveSubgraphs: Finds all the connected induced subgraphs of size up 
%                     to "L0" of a graph known through its adjacency matrix
%                     "A" and containing the subgraph "Subgraph"
%
% INPUTS: A: adjacency matrix of the graph, preferably sparse
%         Aundir: adjacency matrix of the undirected version of the graph
%         directed: true if graph is directed, false otherwise
%         L0: maximum subgraph size, an integer
%         Subgraph: current subgraph, a list of vertices, further
%                   vertices are added to this list 
%         AllowedVert: indicator vector of pruned vertices that may be 
%                      considered for addition to the current subgraph to 
%                      form a larger one
%         Neighbourhood: indicator vector of the vertices that are contained
%                        in the current subgraph or reachable via one edge
%
% OUTPUTS: Primes: list regrouping the contribution of all the subgraphs
%                  found so far
%
% Authors: P.-L. Giscard, N. Kriege, R. Wilson, Septembre 2016


%%%%%%%%%%%%%%%%% FINDING AND COUNTING THE NEIGHBOURS %%%%%%%%%%%%%%%%%%%%%
% Finding allowed neighbours and counting all neighbours (including
% non-allowed ones) of the current subgraph   

% Counts the neighbours of the subgraph that are not contained in the
% subgraph itself
L = length(Subgraph);
NeighboursNumber = nnz(Neighbourhood) - L;

% Gets the subgraph contribution to the prime list
Primes = PrimeCount(A,directed,L0,Subgraph,NeighboursNumber,Primes);

% If the subgraph has reached the target size, stop looking for larger subgraphs
if L==L0 
    return
end    

% Indices of vertices that are both, in the neighbourhood and still allowed
% the vertices in the current subgraph are not allowed
Neighbours = find(Neighbourhood & AllowedVert);

%%%%%%%%%%%%%%%%%%%%%%%% SUBGRAPH CONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for j=1:length(Neighbours) % Adds each neighbour found above to Subgraph to form a new subgraph
    
    v = Neighbours(j);      % The vertex to be added to the subgraph
    Subgraph(L+1) = v;      % New subgraph                      
    AllowedVert(v) = false; % Vertex just added to new subgraph becomes forbidden
    newNeighbourhood = Neighbourhood + Aundir(v,:); % Calculate new neighbourhood
  
    % Continue forming larger subgraphs
    Primes = RecursiveSubgraphs(A,Aundir,directed,L0,Subgraph,AllowedVert,Primes,newNeighbourhood);
end
   
end


function Primes = PrimeCount(A,directed,L0,Subgraph,NeighboursNumber,Primes)
% PrimeCount: Calculates the contribution to the combinatorial sieve of a
%             given subgraph. This function is an implementation of the 
%             Eq. (2), extracting prime numbers from connected 
%             induced subgraphs
%
% INPUTS: A: adjacency matrix of the graph, preferably sparse
%         Aundir: undirected version of the graph adjacency matrix
%         directed: true if graph is directed, false otherwise
%         L0: maximum subgraph size, an integer
%         Subgraph: current subgraph, a list of vertices, further
%                   vertices are added to this list 
%         Primes: list regrouping the contributions of all subgraphs
%                 considered earlier
%
% OUTPUTS: Primes: list with the contributions of all subgraphs so far and 
%                  now including the contribution of the subgraph
%                  passed to this function
%
% Authors: P.-L. Giscard, N. Kriege, R. Wilson, Septembre 2016

SubgraphSize = size(Subgraph,2); % Evaluates Subgraph size                         
x = A(Subgraph,Subgraph);        % Extracts the adjacency matrix x of the subgraph

if directed
    xeig = eig(full(x)); % On directed graphs, fastest way to the eigenvalues is by making x full as x is always going to be small in size
else
    xeig = eig(x);       % On undirected graphs, x is kept sparse, this is faster but not possible if x is directed
end

xS = xeig.^SubgraphSize;                        % List of eigenvalues, each to power |H| to later compute Trace(A_H^k) recursively                            
mk = min(L0 , NeighboursNumber + SubgraphSize); % Maximum value of k yielding a relevant non-zero Binomial(N(H),k-|H|)
BinomialCoeff = 1;                              % The initial value of Binomial(N(H),k-|H|), computed iteratively to avoid the overhead of calling nchoosek

for k = SubgraphSize : mk-1
    Primes(k) = Primes(k) + (-1)^k/k * BinomialCoeff * (-1)^SubgraphSize * sum(xS); % Combinatorial sieve
    
    % Preparation of next loop: list of eigenvalues to power k is put to
    % power k+1 and value of k in binomial coefficient is increased by 1
    xS=xS.*xeig; 
    BinomialCoeff = BinomialCoeff * (SubgraphSize - k + NeighboursNumber) / (1 - SubgraphSize + k);
end

    % Contribution at maximum k value is done separately to avoid
    % incrementing xS and the binomial coefficient one time too many
    Primes(mk) = Primes(mk) + (-1)^mk/mk * BinomialCoeff * (-1)^SubgraphSize * sum(xS);
    
end


function CleanedMatrix = CleanThis(A)
% CleanThis: removes all the vertices of graph which do not sustain any cycle
%            by recursively removing isolated vertices, sinks and sources 
%            until the matrix is invariant under such removal.
% 
%
% INPUT: A: adjacency matrix (possibly weighted) of a (possibly directed) 
%           graph G
%
% OUTPUT: CleanedMatrix: submatrix of A representing a subgraph of G where 
%                        all vertices sustain at least one cycle 
%                        (backtracks and self-loops included)
%
% Authors: P.-L. Giscard, N. Kriege, R. Wilson, Septembre 2016

CleanedMatrix = A;

while size(CleanedMatrix,1)~=size(MatrixCleaning(CleanedMatrix),1)
      
    CleanedMatrix = MatrixCleaning(CleanedMatrix);

end

end



function A = MatrixCleaning(Adj)
% MatrixCleaning: removes all the vertices of graph which do not sustain any cycle
%            by recursively removing isolated vertices, sinks and sources 
%            until the matrix is invariant under such removal.
% 
%
% INPUT: Adj: adjacency matrix (possibly weighted) of a (possibly directed) 
%           graph G
%
% OUTPUT: A: submatrix of A representing a subgraph of G where 
%            all vertices sustain at least one cycle 
%            (backtracks and self-loops included)
%
% Authors: P.-L. Giscard, N. Kriege, R. Wilson, Septembre 2016

    A = Adj;

    x = ~any(A,1); % Detects the sources
    % Removes them
    A(x,:)=[];
    A(:,x)=[];
    
    x = ~any(A,2); % Detects the sinks
    % Removes them
    A(x,:)=[];
    A(:,x)=[];
    
    clear('x');

end

% Copyright (c) 2017, Pierre-Louis Giscard
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:

%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of York nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
