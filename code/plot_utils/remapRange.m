function x_ = remapRange(x,l1,u1,l2,u2)
% map function to a different range: from [l1,u1] to [l2,u2]
    x_ = l2+(x-l1)./(u1-l1).*(u2-l2);
end