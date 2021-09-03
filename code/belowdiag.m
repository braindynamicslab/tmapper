function bld = belowdiag(X)
%BELOWDIAG extract elements of matrix below diagonal.

bld=X(tril(true(size(X)),-1));
end

