function ax = addBackBlock(ax,labels,cmap,varargin)
%ADDBACKBLOCK add blocks in the background of a plot to label task
%structures.
%   ax = addBackBlock(ax,labels,cmap)
% input:
%   ax: axes handle, default current axes 
%   labels: N-by-1 vector of positive integers, for a N-by-N matrix to be
%   labelled. 
%   cmap: M-by-3 matrix, provide a colormap. M is the number of unique
%   labels.
% output:
%   ax: axes handle
%{
created by MZ, 2-13-2020

%}

if isempty(ax)
    ax = gca;
end

p = inputParser;
p.addParameter('alpha',0.3) % transparency of blocks
p.addParameter('ignore',[]) % a set of labels to be ignored
p.addParameter('xlim',[])
p.addParameter('ylim',[])
p.addParameter('ulabelnames','') % cell array of names of each numeric label
p.addParameter('fontsize',16) % fontsize for task names
p.parse(varargin{:});

par = p.Results;

ulabels = unique(labels);
N_ulabels = length(ulabels);
N_samples = length(labels);

cmap = [cmap, repmat(par.alpha,N_ulabels,1)];

for n=1:N_ulabels
    if all(ulabels(n)~=par.ignore)% if user did not ask to ignore label
        [block_start, ~, block_size] = findtaskn(labels==ulabels(n)); % boundaries of a single task

        for n_ts = 1:length(block_start)
            block_pos = remapRange([block_start(n_ts) block_size(n_ts)],1,N_samples,par.xlim(1), par.xlim(2));
            % -- draw blocks
            rectangle('position', [block_pos(1), par.ylim(1), block_pos(2), range(par.ylim)],...
                'edgecolor', 'none','facecolor',cmap(ulabels==ulabels(n),:),'linewidth',1);
            % --- add text labels (if any)
            if ~isempty(par.ulabelnames)
                text(block_pos(1)+block_pos(2)/2, par.ylim(2)+0.075*range(par.ylim),...
                    par.ulabelnames{n},...
                    'HorizontalAlignment','center','FontWeight','bold',...
                    'FontSize',par.fontsize,'Color',cmap(ulabels==ulabels(n),:))
            end
        end
    end
end


end



