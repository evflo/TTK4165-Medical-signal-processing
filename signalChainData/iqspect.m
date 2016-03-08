function [Pxlin,PxdB,f]=iqspect(x,fs,fdemod,Nfft);

%IQSPECT    Estimate frequency spectrum of IQ-data
%
%  [Pxlin,PxdB,f]=iqspect(x,fs,fdemod,Nfft);
% or:
%  [Pxlin,PxdB,f]=iqspect(x,info,Nfft);
%
% Estimates the power spectrum of IQ-data matrix x.
% Power spectra for each beam (column) is estimated,
% and spectra for all beams are averaged.
% The estimation uses Hamming weighting on the data
% and FFT.
%            
%  Inputs: 
%    x       - IQ-data matrix (see READIQ)
%    fs      - Radial sampling frequency of x 
%    fdemod  - Demodulation frequency of IQ-data
%              (default: fdemod=0)
%    Nfft    - Length of FFT. 
%              (default: Zeropad x to next power of 2)
%    info    - info vector (see READINFO)
%
%  Outputs: 
%    Pxlin   - Power spectrum of x
%    PxdB    - Power spectrum of x in dB
%    f       - frequency axis
%
% If no output arguments are specified, the dB spectrum is
% plotted in a new figure window.
%

%(c) GE Vingmed Ultrasound 1997-99. 
% Language   : Matlab 5.2
% Written by : MM 
% Version no.: 1.4
% Date       : 22.02.99, JK
% Calling    : HAMMING (In Matlab Signal Processing Toolbox)
%              GETPARAM


if nargin<2,
    error('Too few input arguments');
end;

% if getparam('isinfo',fs,0),
%     info=fs;
%     if nargin>3,
%         error('Too many input arguments for this synopsis');
%     end;
%     if nargin<3,
%         Nfft=0;
%     else
%         Nfft=fdemod;
%     end;
%     fdemod=getparam('fdemod',info,0);
%     fs=getparam('fradsampiq',info,0);
% else
%     if nargin<3,
%         fdemod=0;
%     end
%     if nargin<4,
%         Nfft=0;
%     end
% end;

%Append zeroes to the next power of 2
N=max(Nfft , 2^(ceil(log10(size(x,1))/log10(2))) );

%Hamming window on data
w= hamming(size(x,1))*ones(1,size(x,2));
xw = w.*x;

%Estimate FFT spectrum for all beams:
Fx=fftshift(fft(xw,N) );    

% frequency axis
f=fdemod+fs*[-0.5:1/(N-1):0.5];

%Square to get power spectrum:
Px = abs(Fx).^2;

% Normalizing scale factor ==> asymptotically unbiased
KMU = norm(hamming(size(x,1)))^2;
Px = Px/KMU;            

%Average spectrum over all beams:
%Pxlin = mean(Px.',2);
Pxlin = mean(Px.');

%dB spectrum:
PxdB=10*log10(Pxlin);

if ~nargout,
    figure;
    plot(f,PxdB-max(PxdB),'k-');
    set(gca,'xlim',1e6*[floor(min(f)/1e6) ceil(max(f)/1e6)]);
    grid on;
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
    title(['Estimated power spectrum [dB]']);
end;

    


function w = hamming(n)
%Symmetric Hamming window of length N
w = .54 - .46*cos(2*pi*(0:n-1)'/(n-1));
    



