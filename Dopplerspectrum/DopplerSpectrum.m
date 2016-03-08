% Calculate and display  spectrum of Doppler signal
% 05.09.02  Hans Torp

%% Settings (play with these)
Nw=32;              % Window size
Nfft=128;           % FFT length
useHpFilter=1;      % Clutter filter, [0=off, 1=on]
overlap = 0.75;      % Fraction of overlap [0-1]
win = 'hamming';    % Window shape ['rect', 'hamming', 'hanning']

%%
load dopplersignal1.mat;% signal: iq, time: t, xmit frequency f0
figno=1;

% Show IQ signal
figure(figno);set(figno,'DoubleBuffer','on');%prevent flickering
clf;
subplot(3,1,1);
plot(t,real(iq),t,imag(iq));axis('tight');
prf=1/(t(2)-t(1));
pause;

% Clutter filtering (wall filtering)
if useHpFilter,
    [b,a]=butter(4,0.2,'high');
    iq=filter(b,a,iq);
    subplot(3,1,1);
    plot(t,real(iq),t,imag(iq));axis ('tight');
   pause;
end;
s=real(iq);s=s/max(s);sound(s,prf);%listen to Doppler signal

if(strcmpi(win,'rect')) 
    w=boxcar(Nw)';
elseif(strcmpi(win,'hamming'))
    w=hamming(Nw)';%window spectrum analyzer
elseif(strcmpi(win,'hanning'))
    w=hanning(Nw)';%window spectrum analyzer
end

dt=round(Nw*(1-overlap));
Nt=length(iq);
maxiq=max(abs(iq));
Nsp=round((Nt-Nw)/dt);
spectrum=zeros(2*Nfft,Nsp);
dyn=40;gain=-40;
wdisp=0.5*maxiq*[0,w,0];%window for display
t1=1;
for k=1:Nsp,
    subplot(3,1,1);
    hold on;plot(t(t1:t1+Nw+1),wdisp);axis('tight');hold off;
    z=iq(t1:t1+Nw-1).*w;
    G=20*log10(abs(fft(z,Nfft))');%power spectrum in dB
    G=[G;G];%stack 2 spectra
    spectrum(:,k)=G;
    subplot(3,1,2);
    image(64*(spectrum+gain)/dyn);colormap(gray);
    subplot(3,1,3);plot(G-110);grid;axis([1,2*Nfft,-80,0]);
    pause(0.001);
    %pause;
    t1=t1+dt;
end;

%temporal averaging to reduce variance
Nt=1;
if Nt>1,
    ht=ones(1,Nt)/Nt;
    spectrumF=filter2(ht,spectrum);
    figure(2);
    image(64*(spectrumF+gain)/dyn);colormap(gray);
end;
