function [res,res_out,res_in,gamma_out,gamma_in] = emd2RTLB(A,B,mA,mB)
  % Script to compute 2-TLB on the real line directly from the distortion expression
  % Vectorized, no loops
  % Weight matrices A and B
  % mA, mB are the measures
  % No entropic regularization, just EMD
  
  
   
  
  % RTLB consists of an LOP over a cost matrix, each entry of which is given
  % by a formula.
  % For inner LOPs, need to get inner cost matrices.
  % Do this by looping over rows of A and rows of B (ecc-out case)


% Copyright (c) 2017, Samir Chowdhury and Facundo Memoli. The Ohio State University. All rights reserved. 
%Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
%following conditions are met: 
%• Redistributions of source code must retain the above copyright notice, this list of
% conditions and the following disclaimer. 
%• Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the followingdisclaimer in the documentation and/or other materials provided with the distribution. 
%• Neither the name of the {organization} nor the names of its contributors may be used to endorse or promote products derived 
% from this software without specific prior writtenpermission. 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
  n     = size(A,1);
  m     = size(B,1);
  
  [r,c] = meshgrid(1:n,1:m);
  
  getOuterCostAnon    = @(ii,jj) getOuterCost(ii,jj,A,B,mA,mB);
  
  [eccout_cst,eccin_cst] = arrayfun(getOuterCostAnon,r,c);
  eccout_cst             = eccout_cst';
  eccin_cst              = eccin_cst';
  
  % using 2-TLB, get squared cost matrix
  eccout_cst             = eccout_cst.^2;
  eccin_cst              = eccin_cst.^2;
  
  
  % EMD
  [dist_out,gamma_out] = mexEMD(mA,mB,eccout_cst);
  [dist_in,gamma_in] = mexEMD(mA,mB,eccin_cst);
  
  %gamma_out
  %gamma_in
  
  % square root
  dist_out = sqrt(dist_out);
  dist_in  = sqrt(dist_in);
  
  % calculate best lower bound
  res_out = dist_out/2;
  res_in  = dist_in/2;
  res = max(res_out,res_in);
    
  
  end






