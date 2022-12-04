function [res,res_out,res_in,gamma_out,gamma_in] = emd2RTLB_unih(A,B)
  % Script to compute 2-TLB on the real line directly from the distortion expression
  % Vectorized, no loops
  % Weight matrices A and B
  % Here we assume node distributions are uniform.
  % No entropic regularization, just EMD
  
  % RTLB consists of an LOP over a cost matrix, each entry of which is given
  % by a formula.
  % For inner LOPs, need to get inner cost matrices.
  % Do this by looping over rows of A and rows of B (ecc-out case)
  
  %{
  This is a simplified version by Mengsen Zhang (circa 9/6/2019), adapted
  from the original code of Chowdhury and Memoli (2018;
  https://github.com/samirchowdhury/GWnets.git). This a special,
  accelerated version where we assume the measure of each node of each
  network is identical. This allows us to only compute the merged
  distribution between two network once, thus, greatly shortens the
  computation time. 
  In constrast to the function "emd2RTLB_uni.m", this version compare
  distributions in GPU via the mex function "compDist_uni.m".
  Due to GPU cache limit, this function is not designed to handle two
  matrices greater than 4000x4000.
  %}
  
  NA     = size(A,1);
  NB     = size(B,1);
  
  [idxA, idxB, dCM, mA, mB] = findIdx(NA,NB); 
  NdCM = length(dCM);
  % -- outgoing 
  eccout_cst = compDist_uni(sort(A,2)',sort(B,2)', idxA, idxB, dCM, NA, NB, NdCM);
  % -- incoming
  eccin_cst = compDist_uni(sort(A,1),sort(B,1), idxA, idxB, dCM, NA, NB, NdCM);
  
  % EMD
  [dist_out,gamma_out] = mexEMD(mA,mB,eccout_cst'.^2);
  [dist_in,gamma_in] = mexEMD(mA,mB,eccin_cst'.^2);

  % calculate best lower bound
  res_out = sqrt(dist_out)/2;
  res_in = sqrt(dist_in)/2;
  res = max(res_out,res_in);
    
  
end

function [idxA, idxB, dCM, mA, mB] = findIdx(NA,NB)
    mA = ones(NA,1)/NA;
    mB = ones(NB,1)/NB;
    cmA = cumsum(mA);
    cmB = cumsum(mB);
    sorted_all = union(cmA,cmB);
  
    %next step prevents bugs from rounding errors, needs Matlab to run
    sorted_all = uniquetol(sorted_all, 10^(-10));
    sorted_all = [0;sorted_all];

    % -- Mengsen's version
    idxA=sum(cmA(:)'<=sorted_all(1:end-1),2)+1;
    idxB=sum(cmB(:)'<=sorted_all(1:end-1),2)+1;  
    dCM = diff(sorted_all);
end




