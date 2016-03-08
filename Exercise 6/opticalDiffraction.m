% simulation of optical diffraction in a photo camera
% 21.01.03 Hans Torp


%% object
o=double(imread('robaat.jpg'));
if ndims(o)>2,o=mean(o,3);end;%convert to grayscale if color image
pix=5e-6;%pixel size
[M,N]=size(o);
x=(0:N-1)*pix;
y=(0:M-1)*pix;
figure(1);clf;colormap(gray(255));
subplot(2,2,1);image(x,y,o);axis image;

%% Object-Fourier
O=fft2(o);
O=fftshift(O);
fmax=1/pix;
fx=fmax/N*(0:N-1)-fmax/2;%spatial frequency x direction
fy=fmax/M*(0:M-1)-fmax/2;%spatial frequency y direction
gain=-60;dyn=40;
ObjdB=20*log10(abs(O));
subplot(2,2,2);image(fx,fy,255*(ObjdB+gain)/dyn);axis image;

%% Point spread function Fourier
%create circular mask in Fourier domain
kx=ones(M,1)*2*pi*fx; %matrix with wavenumber in x-direction
ky=2*pi*fy'*ones(1,N); %matrix with wavenumber in x-direction

R=10e-3;%focal depth
a=2.5e-3;%aperture
fnumber=R/a
lambda=0.8e-6;%wavelength for light
kaperture=2*pi/fnumber/lambda/2;% radius of aperture function in k-space

k=sqrt(kx.^2+ky.^2);%matrix with magnitude of wavenumber
mask=zeros(M,N);
mask(find(k<kaperture))=1;
%defocusing
r=1.0e-3;%displacement of of focal plane. set to zero for max focussing

deltaMax=(a/2)^2*r/2/R^2;%N. Wright equation 2.24
fiMax=2*pi*deltaMax/lambda;%Maximum phase error

mask=mask.*exp(i*fiMax*(k/kaperture).^2);
subplot(2,2,4);imagesc(fx,fy,real(mask));axis image;
pause;
%% Image in Fourier domain
Odiff=mask.*O;
ObjdB=20*log10(abs(Odiff));
subplot(2,2,4);image(fx,fy,255*(ObjdB+gain)/dyn);axis image;

 %% Image in spatial domain
odiff=ifft2(fftshift(Odiff));
subplot(2,2,3);image(x,y,abs(odiff));axis image;

