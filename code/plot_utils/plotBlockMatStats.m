function [blk_mean,blk_std,blk_meanAsym,blk_stdAsym] = plotBlockMatStats(A,pt,varargin)
%PLOTBLOCKMATSTATS plot the statistics of blocks of a square matrix, which
%includes block mean, block standard deviation, asymmetry of the mean,
%asymmetry of the standard deviation. Note that the statistics are plotted
%proportional to the size of the original blocks. 
%   [blk_mean,blk_std,blk_meanAsym,blk_stdAsym] = plotBlockMatStats(A,pt,...)
% input:
%   A: a N-by-N matrix.
%   pt: a N-by-1 vector of integer indices, creating a partition of the
%   rows/cols of A. 
% output:
%   blk_mean: a N-by-N matrix, each element in a block is equal to the
%   average within block
%   blk_std: a N-by-N matrix, each element in a block is equal to the
%   standard deviation (std) within block
%   blk_meanAsym: a N-by-N matrix, each element in a block is equal to the
%   asymmetry of block-mean (asymOpt=afterStats), or the mean of the
%   asymmetry (asymOpt=beforeStats). See NOTE for the measure of asymmetry.
%   blk_stdAsym: a N-by-N matrix, each element in a block is equal to the
%   asymmetry of block-std (asymOpt=afterStats), or the std of the
%   asymmetry (asymOpt=beforeStats). See NOTE for the measure of asymmetry.
% -------------------------------------------------------------------------
% NOTE:
% - "asymmetry" of a matrix X is measure by the positive part of X-X'. The
% negative part is ignored because X-X' is antisymmetric, thus the positive
% part contains all the information. 
%{
created by Mengsen Zhang <mengsenzhang@gmail.com> 3-3-2020
%}
p = inputParser;
p.addParameter('ptnames',unique(pt)); % name of blocks
p.addParameter('x',[]);% x-coordinates
p.addParameter('y',[]);% y-coordinates
p.addParameter('xlabel',''); 
p.addParameter('ylabel','');
p.addParameter('tickOpt','coordinate');% how to label ticks: coordinate, or block
p.addParameter('blockcmap',[]);% colormap for blocks, if non-empty, diagonal blocks will be highlighted with colored boxes
p.addParameter('asymOpt','beforeStats');% when to calculate the asymmetry: 'beforeStats', 'afterStats' (before or after calculating the mean & std)
p.parse(varargin{:});

par = p.Results;

% --- compute block-aggregated statistics
[C, pt_C]=blockMatrix(A, pt);
blk_mean = cell2mat(cellfun(@(x) repmat(mean(x(:)),size(x)),C,'UniformOutput',0));
blk_std = cell2mat(cellfun(@(x) repmat(std(x(:)),size(x)),C,'UniformOutput',0));
mean_title = 'mean path length (a.u.)';
std_title = 'std of path length (a.u.)';
switch par.asymOpt
    case 'beforeStats'% calculating asymmetry before block-wise statistics
        AAsym = A-A';
        AAsym = AAsym .* (AAsym>0);
        CAsym = blockMatrix(AAsym,pt);
        blk_meanAsym = cell2mat(cellfun(@(x) repmat(mean(x(:)),size(x)),CAsym,'UniformOutput',0));
        blk_stdAsym = cell2mat(cellfun(@(x) repmat(std(x(:)),size(x)),CAsym,'UniformOutput',0));
        meanAsym_title = 'mean of asymmetry (a.u.)';
        stdAsym_title = 'std of asymmetry (a.u.)';
    case 'afterStats'% calculating asymmetry after block-wise statistics
        blk_meanAsym = blk_mean - blk_mean';
        blk_meanAsym = blk_meanAsym .* (blk_meanAsym>0);
        blk_stdAsym = blk_std - blk_std';
        blk_stdAsym = blk_stdAsym .* (blk_stdAsym>0);
        meanAsym_title = 'asymmetry of mean (a.u.)';
        stdAsym_title = 'asymmetry of std (a.u.)';
end


% --- plotting
figure('position',[10,10,1000,1000]); 
% - average path length
subplot(2,2,1)
imagesc(par.x,par.y,blk_mean)
colormap gray
colorbar
axis square
xlabel(par.xlabel)
ylabel(par.ylabel)
title(mean_title)
addDiagBlock(gca,pt,par.blockcmap);
switch par.tickOpt
    case 'block'
        xytickBlock(gca,pt,par.ptnames);
end

% - std of path length
subplot(2,2,2)
imagesc(par.x,par.y,blk_std)
colormap gray
colorbar
axis square
xlabel(par.xlabel)
ylabel(par.ylabel)
title(std_title)
addDiagBlock(gca,pt,par.blockcmap);
switch par.tickOpt
    case 'block'
        xytickBlock(gca,pt,par.ptnames);
end

% - average path length
subplot(2,2,3)
imagesc(par.x,par.y,blk_meanAsym)
colormap gray
colorbar
axis square
xlabel(par.xlabel)
ylabel(par.ylabel)
title(meanAsym_title)
addDiagBlock(gca,pt,par.blockcmap);
switch par.tickOpt
    case 'block'
        xytickBlock(gca,pt,par.ptnames);
end

% - std of path length
subplot(2,2,4)
imagesc(par.x,par.y,blk_stdAsym)
colormap gray
colorbar
axis square
xlabel(par.xlabel)
ylabel(par.ylabel)
title(stdAsym_title)
addDiagBlock(gca,pt,par.blockcmap);
switch par.tickOpt
    case 'block'
        xytickBlock(gca,pt,par.ptnames);
end


end

