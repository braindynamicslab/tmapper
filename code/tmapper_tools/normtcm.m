function tcm_n = normtcm(tcm,varargin)
%NORMTCM normalize temporal connectivity matrix by its maximal value that
%is not Inf.
%   tcm_n = normtcm(tcm)
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 2019 ~
modifications:
(10-13-2020) add more options for normalization
%}
p = inputParser;
p.addParameter('normtype','max'); % normalize by max (default) or norm
p.addParameter('infreplace','max'); 
p.parse(varargin{:});
par = p.Results;

% -- handle inf
switch par.infreplace
    case 'max'
        tcm(tcm==Inf) = max(tcm(tcm<Inf));
    case 'nan'
        tcm(tcm==Inf) = nan;
end

% -- normalize
switch par.normtype
    case 'max'
        normfactor = max(tcm(:));
    case 'norm'
        normfactor = norm(tcm(:));
end

if normfactor~=0
    tcm_n = tcm/normfactor;
else
    tcm_n = tcm;
end
end

