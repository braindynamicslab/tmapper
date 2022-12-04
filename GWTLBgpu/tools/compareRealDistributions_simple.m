function dist = compareRealDistributions_simple(sorted_vA,sorted_vB,cmA,cmB)
  % computes the 2-TLB linear program over the real line using generalized inverse     
  % Integrating over the interval [0,1], 
  % functions of interest are the CDFs.
  % Only care about the places where cmA, cmB change
  sorted_all = union(cmA,cmB);
  
  %next step prevents bugs from rounding errors, needs Matlab to run
  sorted_all = uniquetol(sorted_all, 10^(-10));
  sorted_all = [0;sorted_all];
 
  % -- Mengsen's version
  idxA=sum(cmA(:)'<=sorted_all(1:end-1),2)+1;
  idxB=sum(cmB(:)'<=sorted_all(1:end-1),2)+1;  
  dist = sqrt(sum(abs(sorted_vA(idxA) - sorted_vB(idxB)).^2.*diff(sorted_all)));
end