%linear image system demo
%20.01.03 Hans Torp
obj =[... %object 2D matrix
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     4     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     8     16    0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0;
     0     0     0     0     0     0     0     0     0     0     0     0     0];

psf = [ 0   0.5  0; %point spread function
        0.5  1  0.5;
        0  0.5  0];

im=conv2(obj,psf,'same');
disp(im);
figure(1);colormap(gray(16));
subplot(2,1,1);image(obj);colorbar;
axis image;
subplot(2,1,2);image(im);colorbar;
axis image;

O=fft2(obj);
[N,M]=size(obj);
Psf=fft2(psf,N,M);
figure(2);
imagesc(fftshift(abs(Psf)));colormap(gray);