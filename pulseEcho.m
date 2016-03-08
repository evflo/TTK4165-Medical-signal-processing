%% Pulse echo 1D imaging - simulation
% Prerequisite for exercise 3 in MEDT4165 
% Simulate received signal from an 1 D object defined by variation in acoustic impedance
% Calculates the total pulse echo response h = pel x pxd x hd x hxd ; where
% "x" means convolution
% Original version: 19.01.02  Hans Torp
% Updated:          27.01.03  HT
%                   10.12.09  Sten Roar Snare

clear all;clc;

%% physical constants
c=1540;  %speed of sound

zmax=0.04;%max simulation depth
dz=10e-6; %depth increment
z=0:dz:zmax;%depth axis
z=z';%make z into a vector

t=2/c*z;%time axis
dt=t(2)-t(1);
fs=1/dt;%sampling freq. 
tmax=t(end);


%% define and plot transmitted pulse
f0=4e6;
lambda=c/f0;
Tp=2e-6;%pulse length %1/f0*Nperiods (1/f0*8)
tp=0:1/fs:Tp;tp=tp';
pel=sin(2*pi*f0*tp);% Electrical transmit pulse

figure(1);clf;
plot(tp,pel);
title('Transmit pulse');xlabel('Time [s]');ylabel('Signal amplitude');
%% define and plot transducer impulse response
fc=2.5e6;B=1.2e6;%center frequency and bandwidth of transducer

[bxd,axd]=butter(2,2*[fc-B/2,fc+B/2]/fs);%butterworth bandpass filter

txd=0:dt:3e-6;txd=txd'; %time axis for impulse
impulse=zeros(size(txd));impulse(1)=1; %create impulse

hxd=filter(bxd,axd,impulse); %calculate probe impulse response

figure(2);clf;
plot(txd,hxd);
title('Probe impulse response');xlabel('Time [s]');ylabel('Signal amplitude');

%% calculate the txrx impulse response, including scattering

pAc=conv(pel,hxd);% transmitted acoustical pulse
hd = 0.5*[-1;1];  % differentiation operator to account for scattering
h=conv(pAc,hd);   % differentiation of Ac. impedance
h=conv(h,hxd);    % transducer receiver response = transmit response hxd 

figure(3);clf;
subplot(3,1,1);
t_h = dt*(0:(length(h)-1));
plot(t_h*1e6,h);
title('Txrx impulse response');xlabel('Time [\mus]');ylabel('Amplitude');

%% simulate 1D objects using the txrx impulse response

% Object defined by acoustic impedance Z as a function of depth. Example:
% Object 1: water/fat interface at depth z1
%           fat/water interface at depth z2
z1=0.02;
z2=0.03;

Zw=1.48; % [kg/m^2/s] acoustic impedance of water
Zf=1.66; % [kg/m^2/s] acoustic impedance of human fat tissue
Zair=0.0004;
Zbone=5.5;

Z1=Zw;
Z2=Zf;

Z=ones(size(z))*Z1;
Z(1+round(z1/dz):1+round(z2/dz))=Z2;


R_c = (Z2-Z1)/(Z1+Z2);
R_p = 1;

obj=Z/mean(Z);% 1D object function for the pulse echo system

s=conv(h,obj);%received signal

%Plot results
figure(3);%Make sure figure 3 is active
subplot(3,1,2);
plot(1000*z,obj);
title('Object function');xlabel('Distance [mm]');ylabel('Object "strength"');
subplot(3,1,3);
t_s = dt*(0:length(s)-1);
plot(t_s*1e6,s);
title('Received scattering signal');xlabel('Time [\mus]');ylabel('Amplitude');

%% Add noise and present final results
s=s(1:end-length(h)+1); %Remove time samples corresponding to system txrx response 
an=0.25e-3; %noise amplitude
sn=an*randn(size(s));%Gaussian white noise
s=s+sn;%add thermal noise

figure(4);clf;
subplot(3,1,1);
plot(1000*z,Z);
axis([0 40 1 2])
xlabel('depth z [mm]');ylabel('Object "strength"');title('Object function');

smax=max(abs(s));
subplot(3,1,2);
plot(1e6*t,s);
xlabel('time [\mus]');ylabel('Amplitude');title('Received signal');
axis([1e6*t(1),1e6*t(end),-smax,smax]);

amp=abs(hilbert(s));%amplitude of echo signal s
logamp=20*log10(amp);
subplot(3,1,3);
plot(1e6*t,logamp);
axis([1e6*t(1),1e6*t(end),-70,0]);
xlabel('time [\mus]');ylabel('Power [dB]');title('Received signal power');

