function FCz = fisherzFC(FC)
%FISHERZFC Fisher-Z transform for functional connectivity matrix. 
%   which is simple atanh of subdiagonal elements

FCz = atanh(FC(tril(true(size(FC)),-1)));
end

