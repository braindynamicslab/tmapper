function prX = phaseRand(X,seed)
%PHASERAND Phase randomization of signals, following Liegeois et al (2017,
%NeuroImage).
%   prX = phaseRand(X,seed)
% input:
%   X: orginal signals, a N-by-M matrix, where N is the number of sample
%   points and M is the number of channels. 
%   seed: seed for random number generator. If left empty, the global
%   randStream is used. 
% output:
%   prX: phase randomized signal. 
%{
created by MZ, 02/07/2020
modifications:

%}

if ~(nargin<2 || isempty(seed))
    rng(seed);
end

Nsp = size(X,1); %number of samples
fX = fft(X,Nsp,1); % FFT of X

% -- generate random phase
Npr = floor((Nsp-1)/2);
phi = 2*pi*rand(Npr,1);% random phase
fX(2:Npr+1,:) = abs(fX(2:Npr+1,:)).*exp(1i*(angle(fX(2:Npr+1,:))+phi));% add random phase to half of the spectra
fX(Nsp:-1:Nsp-Npr+1,:) = abs(fX(Nsp:-1:Nsp-Npr+1,:)).*exp(1i*(angle(fX(Nsp:-1:Nsp-Npr+1,:))-phi));% add complex conjugate

% -- get back time series
prX=ifft(fX,Nsp,1); 

end

