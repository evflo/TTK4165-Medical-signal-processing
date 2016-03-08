%2D fouriertransform demo
% 08.02.02  Hans Torp
% Impulse in the 2D Fourier plane transformes to plane waves in spaceread(
Nf=64;
S=zeros(Nf,Nf);
clf;
figure(1); whitebg(0);
colormap(gray);
subplot(1,2,1);
S=zeros(Nf,Nf);
ax=[-Nf/2:Nf/2-1]
image(ax,ax,S);axis('image');grid;

y=1;
while y>0,
    subplot(1,2,1);
    [x,y]=ginput(1);
    x=round(x+Nf/2+1);y=round(y+Nf/2+1);
    S=zeros(Nf,Nf);
    if y>0, S(y,x)=1;end;
    imagesc(ax,ax,S);axis('image');grid;
    s=ifft2(S);
    s=ifft2(fftshift(S));
    subplot(1,2,2);
    imagesc(real(s));axis('image');
end;


