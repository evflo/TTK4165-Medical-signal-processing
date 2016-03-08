load lasse_1;


f0=s.iq.TxFreqIQ_Hz;
fh = 2*f0;
frs=s.iq.frsIQ_Hz;
fdemod=s.iq.fDemodIQ_Hz;
B=0.8e6; %0.8 MHz filter bandwidth
N=10; %filterorder
iqf=rectfreq(iq,frs,fdemod,f0,B,N);

iqf2=rectfreq(iqf,frs,fdemod,fh,B,N);

%theta axis
th1=s.tissue.StartAngle_rad;
dth=s.tissue.AngleIncrement_rad;
Nbeams=s.tissue.Beams;
th = th1:dth:((Nbeams-1)*dth+th1);
%range axis
r1=s.tissue.StartDepth_m;
dr=s.tissue.DepthIncrement_m;
Nr=s.tissue.Samples;
r=r1:dr:((Nr-1)*dr+r1);

[tsc,xax,zax]=scanconvert(iqf2,r,th);
image(xax,zax,tsc);
colormap(gray(256));
axis image;
xlabel('cm');ylabel('cm');
title('Scanconverted image');