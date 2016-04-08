function [f,H]=clutterfilterrespons(Fm);

%CLUTTERFILTERRESPONS Clutter Filter Respons
%
%   [f,H]=clutterfreqresp(Fm)
%
%   Fm - Filter matrix (REGRESJONSFILTER)
%   Nf - Number of frequency samples, default=256;
%
%   f  - Frequency axis of frequency response
%   H  - Power of frequency response
%

%Written by : Hans Torp
%Last update: 26.03.2001, JK

Nf=256;
f=[-0.5 :1/Nf:0.5];
Fmf=fft(Fm,Nf)';
H=mean(conj(Fmf).*Fmf);
H=abs(fftshift(H'))';
H=[H,H(1)];
