function [block_start, block_end, block_size] = findtaskn(ind)
% FINDTASKN find the beginning and end of blocks based on an indicator
% function "ind" of when a specific task occurs in the entire duration of
% recording.
%   [block_start, block_end, block_size, block_type] = findtaskn(ind)
% input:
%   ind: N-by-1 vector, indicator function of 0s and 1s, where 1 indicate
%   the presence of a specific task N is number of samples.
% output:
%   block_start: M-by-1 vector, first sample of each block. M is the number
%   of blocks of the task (ind=1).
%   block_end: M-by-1 vector, last sample of each block. 
%   block_size: M-by-1 vector, number of samples in each block.
%{
Author: Mengsen Zhang <mengsenz@stanford.edu> 2-26-2020
%}

block_start = find(diff([0; ind(:)])==1);
block_end = find(diff([ind(:);0])==-1);

if length(block_start) ~= length(block_end)
    error('Block structure incomplete.')
end

block_size = block_end - block_start + 1;
end