function [P,f]=puls_fft(RF,Ts,Nfft,freq_range)
%Enkel FFT for estimering av FFT spekter for hydrofonpuls
%
% function [P,f]=puls_fft(RF,Ts,Nfft,frequency_range)
%
% Input:
%  RF:         Pulsform
%  Ts:         Samplingsintervall i sekund
%  Nfft:       Lengda på FFT, typisk lik lengda på signal
%  freq_range: Frekvensområde, typisk [0 5e6]
%
% Output:
%  P: Amplitudespekter, ikkje dB skala.
%  f: Frekvensakse
%
% JK, 30jan01
%

if nargin<4,
	freq_range=[0 1/Ts/2];
end;

%Fast fourier transform of the RF signal
P=abs(fft(RF,Nfft));

%The corresponding frequency axis
f=([0:Nfft-1]/Nfft)/Ts;

%Limit the spectrum to the desired frequency range
minind=min(find(f>min(freq_range)));
maxind=max(find(f<max(freq_range)));

P=P(minind:maxind);
f=f(minind:maxind);