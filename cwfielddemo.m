%Continuous wave from point sources
%15.01.01  Hans Torp
% last update 27.12.03

f0=2e6;%frequency
c=1540;%speed of sound
k=2*pi*f0/c;%wave number

xmax=10e-3;% +/- size  x direction
zmax=30e-3;%size  z direction
Nx=300;%number of points in x direction
Nz=400;%antall punkter i zretning

xaxis=linspace(-xmax,xmax,Nx);
zaxis=linspace(0,zmax,Nz)';
x=ones(Nz,1)*xaxis;
z=zaxis*ones(1,Nx);

p=0*x;
image(xaxis,zaxis,real(p));
colormap(gray);axis('image')
disp('Click to set point source position');
disp('press return to finnish');
[x0,z0]=ginput(1);
while ~isempty(x0),
   r=sqrt((x-x0).^2+(z-z0).^2);
   p=p+exp(i*k*r)./r;
   image(xaxis,zaxis,0.05*real(p));
   axis('image');
   [x0,z0]=ginput(1);
end;

%animation
set(1,'Doublebuffer','on');%to prevent flickering
dt=0.05e-6;
w0=2*pi*f0;
for n=1:200,
    p=p.*exp(-i*w0*dt);
    image(xaxis,zaxis,0.05*real(p));
   axis('image');
   pause(0.01);
end;

return;
%plane source
z0=5e-3;
dx0=0.3e-3;
p=0*x;
for x0=-10e-3:dx0:10e-3,
   r=sqrt((x-x0).^2+(z-z0).^2);
   p=p+exp(i*k*r)./r;
   image(xaxis,zaxis,0.05*real(p));
   axis('image');
   pause(0.01);
end;    
    