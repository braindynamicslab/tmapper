function s = par2filename(par)
%PAR2FILENAME from parameter assignment (cell array of Name-Value pairs) to
%proper filename as a string
%   

Narg=length(par)/2;
s='';
for n=1:Narg
    s=[s, sprintf('%s%s_',par{2*n-1},num2fname(par{2*n}))];
end
s(end)=[];

end

function n=num2fname(n)
    n=strrep(num2str(n,'%g'),'.','-');
end

