function [C, up] = subMatrix(A,p)
%SUBMATRIX convert matrix to a cell array of matrices, each of which is a
%submatrix of the orignal matrix, given a partition of the rows/columns of
%the matrix.
%   [C, up] = subMatrix(A,p)
% input:
%   A: N-by-N square matrix.
%   p: N-by-1 vector. Integer labels of each column. 
% output:
%   C: M-by-M cell array. M is the number of unique labels in "p". C_ij
%   contains the submatrix A(p==p_i, p==p_j), where p_k is the k-th
%   smallest unqiue label. 
%   up: N-by-1 vector of unique labels. 
%{
Author: Mengsen Zhang <mengsenzhang@gmail.com> 02-26-2020
%}

up = unique(p);
N_up = length(up);

C = cell(N_up);

for i=1:N_up
    for j=1:N_up
        C{i,j} = A(p==up(i),p==up(j));
    end
end

end

