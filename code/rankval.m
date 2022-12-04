function r = rankval(val)
%RANKVAL return ranks of values, from 1 to N, the smallest to largest
%number. 
%  
[~,idx] = sort(val);
r = (1:length(val))';
r(idx) = r;
end

