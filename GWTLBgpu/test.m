% test.m
addpath(genpath('.'))

format long
clear all
close all
clc

N = 300; 
M = 237; 
A = rand(N); 
B=rand(M); 

% -- uniformly distribution over nodes 
mA = ones(N,1)/N; 
mB = ones(M,1)/M;

% -- random distribution over nodes
% mA = rand(N,1); mA = mA/sum(mA);
% mB = rand(M,1); mB = mB/sum(mB);

% -- CPU ALG's of different complexity and speed
% disp("--- CPU original version ---")
% tic
% disp(emd2RTLB(A,B,mA,mB))
% toc

disp("--- CPU general ---")
tic
disp(emd2RTLB_simple(A,B,mA,mB))
toc

disp("--- CPU uniform ---")
tic
disp(emd2RTLB_uni(A,B))
toc

% -- GPU ALG's of different complexity

disp("--- GPU general ---")
tic
disp(emd2RTLB_hetero(A,B,mA,mB));
toc

disp("--- GPU uniform ---")
tic
disp(emd2RTLB_unih(A,B));
toc
