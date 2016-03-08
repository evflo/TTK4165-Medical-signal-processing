function [im,xax,zax]=scanconvert(imIn,r,th,nx,nz)

% Scanconversion
% [im,xax,zax]=scanconvert(imIn,r,th,nx,nz)
%
% feb 2010, Torbj?rn Hergum
%
if nargin<5,
	nz=size(imIn,1);
end
if nargin<4,
	nx=300;
end
if nargin<3,
	error('You must specify the axis of the data')
end

[r2,th2] = meshgrid(r,th);
[ySc,xSc] = pol2cart(th2,r2);

%Specify position in cartesian coordinates
xax=linspace(min(xSc(:)),max(xSc(:)),nx);
zax=linspace(r(1),r(end),nz);
[y2,x2]=meshgrid(zax,xax);

warning off; % If origo is included griddata will return a useless warning
im = griddata(xSc',ySc',imIn,x2',y2');
warning on;
