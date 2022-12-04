function [res,res_out,res_in,gamma_out,gamma_in] = emd2RTLB_hetero(A,B,mA,mB)
  % Script to compute 2-TLB on the real line directly from the distortion expression
  % Vectorized, no loops
  % Weight matrices A and B
  % mA, mB are the measures
  % No entropic regularization, just EMD
  
  % RTLB consists of an LOP over a cost matrix, each entry of which is given
  % by a formula.
  % For inner LOPs, need to get inner cost matrices.
  % Do this by looping over rows of A and rows of B (ecc-out case)
% -------------------------- what is this ------------------------------- %
% This is a modified version by Mengsen Zhang (Sep, 2019) of the original
% script "emd2RTLB.m" by Chowdhury and Memoli (2018;
% https://github.com/samirchowdhury/GWnets.git). Here the original function
% "compareRealDistributions.m" for calculating the cost matrices has been
% substituted by its GPU counterpart ("compDist_merge"). 
% This accelerated version, however, requires the size of the matrices
% (n+m) <= 3071. This is limited by common cache size of consumer-grade
% NVIDIA graphic cards. 
% In addition, the GPU version uses single precision operations, so that
% the agreement to the CPU counterpart is around the 8-th significant
% digit. 
% The usage of this function is otherwise identical with the original
% version. 
% ----------------------------------------------------------------------- %
  
  n     = size(A,1);
  m     = size(B,1);
  
  % -- outgoing 
  [sorted_A_out, CMA_out] = sorteachrow(A,mA);
  [sorted_B_out, CMB_out] = sorteachrow(B,mB);
  eccout_cst = compDist_merge(sorted_A_out',sorted_B_out',CMA_out',CMB_out',n,m);
  % -- incoming
  [sorted_A_in, CMA_in] = sorteachrow(A',mA);
  [sorted_B_in, CMB_in] = sorteachrow(B',mB);
  eccin_cst = compDist_merge(sorted_A_in',sorted_B_in',CMA_in',CMB_in',n,m);
  
  % EMD
  [dist_out,gamma_out] = mexEMD(mA,mB,eccout_cst'.^2);
  [dist_in,gamma_in] = mexEMD(mA,mB,eccin_cst'.^2);
  
  % calculate best lower bound
  res_out = sqrt(dist_out)/2;
  res_in = sqrt(dist_in)/2;
  res = max(res_out,res_in);
  
end

function [sorted_A, CMA] = sorteachrow(A,mA)
    [sorted_A, sorted_idA] = sort(A,2);
    sorted_mA = mA(sorted_idA);
    CMA = cumsum(sorted_mA,2);
end





