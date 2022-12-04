function [res,res_out,res_in,gamma_out,gamma_in] = emd2RTLB_simple(A,B,mA,mB)
  % Script to compute 2-TLB on the real line directly from the distortion expression
  % Vectorized, no loops
  % Weight matrices A and B
  % mA, mB are the measures
  % No entropic regularization, just EMD
  
  % RTLB consists of an LOP over a cost matrix, each entry of which is given
  % by a formula.
  % For inner LOPs, need to get inner cost matrices.
  % Do this by looping over rows of A and rows of B (ecc-out case)
  
  %{
  This is a simplified version by Mengsen Zhang (circa 9/6/2019). Two
  modifications are made comparing to the original version by Chowdhury and
  Memoli (2018; https://github.com/samirchowdhury/GWnets.git). 
  1. sorting matrices A, B and calculating the accumulative probability
  distributions outside of the for-loop in the original version. This is
  done using a new function "sorteachrow".
  2. vectorize the index-finding for-loop in the original function
  "compareRealDistribution.m", which becomes
  "compareRealDistribution_simple.m" used here.
  %}
  
  n     = size(A,1);
  m     = size(B,1);
  
  [r,c] = meshgrid(1:n,1:m);
  
  % -- outgoing 
  [sorted_A_out, CMA_out] = sorteachrow(A,mA);
  [sorted_B_out, CMB_out] = sorteachrow(B,mB);
  eccout_cst = arrayfun(@(ii,jj) compareRealDistributions_simple(...
      sorted_A_out(ii,:)',sorted_B_out(jj,:)',CMA_out(ii,:)',CMB_out(jj,:)'),r,c);
  % -- incoming
  [sorted_A_in, CMA_in] = sorteachrow(A',mA);
  [sorted_B_in, CMB_in] = sorteachrow(B',mB);
  eccin_cst = arrayfun(@(ii,jj) compareRealDistributions_simple(...
      sorted_A_in(ii,:)',sorted_B_in(jj,:)',CMA_in(ii,:)',CMB_in(jj,:)'),r,c);
  
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





