function [res,res_out,res_in,gamma_out,gamma_in] = emd2RTLB_uni(A,B)
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
  %}
  
  NA     = size(A,1);
  NB     = size(B,1);
  
  [r,c] = meshgrid(1:NA,1:NB);
  [idxA, idxB, dCM, mA, mB] = findIdx(NA,NB); 
  % -- outgoing 
  sorted_A_out = sort(A,2);
  sorted_B_out = sort(B,2);
  eccout_cst = arrayfun(@(ii,jj) compdist(...
      sorted_A_out(ii,:)',sorted_B_out(jj,:)',idxA, idxB, dCM),r,c);
  % -- incoming
  sorted_A_in = sort(A',2);
  sorted_B_in = sort(B',2);
  eccin_cst = arrayfun(@(ii,jj) compdist(...
      sorted_A_in(ii,:)',sorted_B_in(jj,:)',idxA, idxB, dCM),r,c);
  
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

function dist = compdist(sorted_vA,sorted_vB,idxA, idxB, dCM)
    dist = sqrt(sum(abs(sorted_vA(idxA) - sorted_vB(idxB)).^2.*dCM));
end




