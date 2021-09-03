function [C,pt_C] = blockMatrix(A,pt)
%BLOCKMATRIX convert matrix to a cell array of matrices, each of which is a
%block of the orignal matrix, given a partition of the rows/columns of the
%matrix.
%   [C,pt_C] = blockMatrix(A,pt)
% input:
%   A: N-by-N square matrix.
%   pt: N-by-1 vector. Integer labels of each column. 
% output:
%   C: M-by-M cell array. C_ij contains the (i,j) block of A.
%   pt_C: N-by-1 vector. Integer labels for each block. 
%{
Author: Mengsen Zhang <mengsenzhang@gmail.com> 2-27-2020
modifications:

%}


d = [true; diff(pt) ~= 0; true];  % TRUE if values change
ptsizes = diff(find(d)); 
C=mat2cell(A,ptsizes,ptsizes);
pt_C = pt(cumsum([1; ptsizes(1:end-1)]));
end

