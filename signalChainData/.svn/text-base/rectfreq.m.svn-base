function iqf=rectfreq(iq,info,f0,B,N,dum);

%RECTFREQ   Bandpass FIR filter for IQ data.
%
%     iqf=rectfreq(iq,info,f0,B,N);
%
%   Rectangular frequency bandpass filter in range direction
%   (first index) of iq components
%   N'th order FIR filter with impulse resp. 
%   = hamming weighted sinc of length N
%
%   iq   iq signal matrix (quadrature demodulated RF signal)
%        (see READIQ)
%   info info vector from readinfo
%   f0   filter center frequency (default=transmit frequency)
%   B    filter bandwidth, (dualsided -6dB, default 1MHz)
%   N    filterorder,should be an even number (default 10)
%
% Alternative synopsis (No defaults allowed!):
%
%   iqf=rectfreq(iq,frs,fm,f0,B,N); 
%
%   frs  radial sampling frequency (frsiq)
%   fm   mixing frequency (fdemod)
%

%(c) GE Vingmed Ultrasound 1997-99. 
% Calling     :HAMMING (Signal Processing Toolbox)
% Language    :Matlab 5.2
% Written by  :14.10.96 HT
% Updated     :19.02.99 JK

if nargin==6, 
    frs=info;
    fm=f0;
    f0=B;
    B=N;
    N=dum;
else
    if nargin<5, N=10;end;
    if nargin<4, B=1e6;end;
    if nargin<3, f0=getparam('txfreqiq',info,0);end;
    fm=getparam('fdemodiq',info);
    frs=getparam('frsiq',info);
end;

shape=[];
if size(iq,3)>1,
    shape=size(iq);
    iq=reshape(iq,[shape(1) prod(shape(2:3))]);
end;

[rpb beams]=size(iq);
dw=2*pi*(f0-fm)/frs;
t=1:rpb;
mix=exp(-i*dw*t');
T=floor((N+1)/2);
N=2*T+1;
b=sinc([-T:1:T]*B/frs).*hamming(N)';
b=b/sum(b);
mixM=mix*ones(1,beams);%down mixing
iqf=filter2(b',iq.*mixM);%lp filtering
iqf=iqf.*conj(mixM);%up mixing

if any(shape), 
    iqf=reshape(iqf,shape);
end;


function w = hamming(n)
%Symmetric Hamming window of length N
w = .54 - .46*cos(2*pi*(0:n-1)'/(n-1));

function y=sinc(x)
y=ones(size(x));
i=find(x);
y(i)=sin(pi*x(i))./(pi*x(i));

