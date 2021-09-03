function cb = CMEcbar(cb)
%CMECBAR Customize colorbar for CME data

cb.Label.String = '';
cb.FontWeight='bold';
set(cb,'xtick',1.4+0.8*(0:4),...
    'xticklabel',{'instruction','rest','memory','video','math'})
end

