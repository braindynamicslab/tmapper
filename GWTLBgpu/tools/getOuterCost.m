function [eccout_cst_ii_jj,eccin_cst_ii_jj] = getOuterCost(ii,jj,...
                                            A,B,mA,mB)
                                    
  vA_out = A(ii,:);
  vB_out = B(jj,:);
  
  vA_in  = A(:,ii);
  vB_in  = B(:,jj);

  % fill in outer cost matrix
  eccout_cst_ii_jj = compareRealDistributions(vA_out,vB_out,mA,mB);
  eccin_cst_ii_jj  = compareRealDistributions(vA_in,vB_in,mA,mB);
end