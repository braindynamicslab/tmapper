function dist = compareRealDistributions(vA,vB,mA,mB)
  % computes the 2-TLB linear program over the real line using generalized inverse     
  [sorted_vA,sorted_idA] = sort(vA);
  [sorted_vB,sorted_idB] = sort(vB);
  
  sorted_mA = mA(sorted_idA);
  sorted_mB = mB(sorted_idB);
  
    
  cmA = cumsum(sorted_mA);
  cmB = cumsum(sorted_mB);
  
  % Integrating over the interval [0,1], 
  % functions of interest are the CDFs.
  % Only care about the places where cmA, cmB change
  sorted_all = union(cmA,cmB);
  
  %next step prevents bugs from rounding errors, needs Matlab to run
  sorted_all = uniquetol(sorted_all, 10^(-10));
  sorted_all = [0;sorted_all];
  N_all = length(sorted_all);
  summand = 0;

 
  for ii = 1:N_all - 1
    % find first index where cmA exceeds sorted_all(ii)
    % use the index to get the generalized inverse of the distribution function
    idxA = find(cmA > sorted_all(ii) , 1);
    ginvA = sorted_vA(idxA);
    %ginvA = cmA(idxA);
    % likewise for B
    idxB = find(cmB > sorted_all(ii) , 1);
    ginvB = sorted_vB(idxB);
    %ginvB = cmB(idxB);
    
    riem_int = (abs(ginvA - ginvB))^2; %2-TLB here
    summand  = summand + riem_int*(sorted_all(ii+1) - sorted_all(ii));
  end
  dist = sqrt(summand);
  
end