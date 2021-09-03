function tasknames = CMEtasknames(tasknum,varargin)
%CMETASKNAMES Name of the tasks for the CME data.
%   
p=inputParser;
p.addParameter('short',false);% use shortened names
p.parse(varargin{:})
par=p.Results;

if par.short
    tasknames = {'inst','rest','mem','video','math'};
else
    tasknames = {'instruction', 'rest', 'memory', 'video', 'math'};
end

if nargin>0 && ~isempty(tasknum)
    tasknames = tasknames(tasknum);
end


end

