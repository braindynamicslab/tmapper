% WGGG
%
%
%
% aim: to calculate modularity score, given a matrix and its community
% assignment
%
% author: saggar@stanford.edu

% Extra comments added by Samir
% Inputs: 
% W         -- a graph adjacency matrix
% m0        -- categorical node labels//community assignment

function mod = calMod(W, m0)



   s = sum(sum(W)); % calculate sum of degrees of nodes (assuming symmetric)
   
   if s == 0
       mod = 0;
       % samir: if there are no edges, consider the network to be neither
       % modular nor non-modular
   else
   
       gamma = 1; % scaling parameter

       B   = (W-gamma*(sum(W,2)*sum(W,1))/s)/s; % normalized form of "modularity matrix"

       B = (B+B.')/2;    % symmetrize
       mod = sum(B(bsxfun(@eq, m0, m0'))); %calculate Q score
   end



end