function M = zerodiag(M)
%ZERODIAG set diagonal to zero
if size(M,1)~=size(M,2), error('This function only takes square matrices.');end
M(logical(eye(length(M))))=0;
end

