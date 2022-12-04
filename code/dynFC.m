function [dfc,dfc_fisherz] = dynFC(X,winsize,lag)
%DYNFC Dynamic Functional Connectivity (DFC)
%   input:
%       X: a N_samples-by-N_channels matrix. 
%       winsize: window size for computing FC (in samples)
%       lag: lag between windows.
%   output:
%       dfc: cell array of DFC
%       dfc_fisherz: N_windows-by-M matrix, DFC after fisher-z
%       transformation. M is the number of elements below diagonal in each
%       FC matrix, M = N_channels*N_channels/2 - N_channels.

[Nsp,Nchn] = size(X);

Nwin = round((Nsp-winsize+lag)/lag);

% -- store DFC
dfc = cell(Nwin,1);
for n = 1:Nwin
    dfc{n}= corrcoef(X(lag*(n-1)+1:lag*(n-1)+winsize,:));
end

dfc_fisherz= cell2mat(cellfun(@(x) fisherzFC(x)', dfc,'uniformoutput',false));
end

