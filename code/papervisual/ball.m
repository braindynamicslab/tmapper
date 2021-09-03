function [bXX,bYY,bZZ] = ball(x,y,z,r,res)
%BALL coordinates of a ball centered at (x,y,z) with radius r at resolution
% res. 

if nargin<5
    res = 100;
end
[bXX,bYY,bZZ]=sphere(res);
% :: cooridnates of the main ball
if length(r)==1
    r=repmat(r,1,3);
end

[bXX,bYY,bZZ]=deal(r(1)*bXX+x,r(2)*bYY+y,r(3)*bZZ+z);
end

