function rf=iq2rf(iq,fdemod,frsiq,frsrf);

%IQ2RF      Convert IQ-data to RF-data
%
% function rf=iq2rf(iq,info,frs);
%
% iq   - matrix of IQ-data (see READIQ)
% info - info vector (see READINFO)
% frs  - desired radial sampling frequency
%        (default = 20MHz)
%

%(c) GE Vingmed Ultrasound 1999
%Written    : 13.07.1998, JK
%Last update: 11.11.1998, JK
%Calling    : INTERP (Matlab Signal Processing Toolbox)
%
%Removing readecho-dependencies. TH 02.2010

if nargin<4, frsrf=20e6;end;
if nargin<3, error('Too few input arguments');end;

I=round(frsrf/frsiq);

order=3;
cutoff=1; %=frsiq/2

[N,M]=size(iq);
rf = zeros(I*N,M);
[N,M]=size(rf);
t=[0:N-1]'/frsrf;
mix=exp(i*2*pi*fdemod*t);
for n=1:M,
   rf(:,n)=sqrt(2)*real(interp(iq(:,n),I,order,cutoff).*mix);
end;





function [odata,b] = interp(idata,r,l,alpha)

if nargin < 3
   l = 4;
end
if nargin < 4
   alpha = .5;
end
if l < 1 | r < 1 | alpha <= 0 | alpha > 1
   error('Input parameters are out of range.')
end
if abs(r-fix(r)) > eps
   error('Resampling rate R must be an integer.')
end
if 2*l+1 > length(idata)
    s = int2str(2*l+1);
    error(['Length of data sequence must be at least ',s, 10,...
'You either need more data or a shorter filter (L).']);
end

% ALL occurrences of sin()/() are using the sinc function for the
% autocorrelation for the input data.  They should all be changed consistently
% if they are changed at all.

% calculate AP and AM matrices for inversion
s1 = toeplitz(0:l-1) + eps;
s2 = hankel(2*l-1:-1:l);
s2p = hankel([1:l-1 0]);
s2 = s2 + eps + s2p(l:-1:1,l:-1:1);
s1 = sin(alpha*pi*s1)./(alpha*pi*s1);
s2 = sin(alpha*pi*s2)./(alpha*pi*s2);
ap = s1 + s2;
am = s1 - s2;
ap = inv(ap);
am = inv(am);

% now calculate D based on INV(AM) and INV(AP)
d = zeros(2*l,l);
d(1:2:2*l-1,:) = ap + am;
d(2:2:2*l,:) = ap - am;

% set up arrays to calculate interpolating filter B
x = (0:r-1)/r;
y = zeros(2*l,1);
y(1:2:2*l-1) = (l:-1:1);
y(2:2:2*l) = (l-1:-1:0);
X = ones(2*l,1);
X(1:2:2*l-1) = -ones(l,1);
XX = eps + y*ones(1,r) + X*x;
y = X + y + eps;
h = .5*d'*(sin(pi*alpha*XX)./(alpha*pi*XX));
b = zeros(2*l*r+1,1);
b(1:l*r) = h';
b(l*r+1) = .5*d(:,l)'*(sin(pi*alpha*y)./(pi*alpha*y));
b(l*r+2:2*l*r+1) = b(l*r:-1:1);

% use the filter B to perform the interpolation
[m,n] = size(idata);
nn = max([m n]);
if nn == m
   odata = zeros(r*nn,1);
else
   odata = zeros(1,r*nn);
end
odata(1:r:nn*r) = idata;
% Filter a fabricated section of data first (match initial values and first derivatives by
% rotating the first data points by 180 degrees) to get guess of good initial conditions
% Filter length is 2*l*r+1 so need that many points; can't duplicate first point or
% guarantee a zero slope at beginning of sequence
od = zeros(2*l*r,1);
od(1:r:(2*l*r)) = 2*idata(1)-idata((2*l+1):-1:2);
[od,zi] = filter(b,1,od);
[odata,zf] = filter(b,1,odata,zi);
odata(1:(nn-l)*r) = odata(l*r+1:nn*r);

% make sure right hand points of data have been correctly interpolated and get rid of
% transients by again matching end values and derivatives of the original data
if nn == m
    od = zeros(2*l*r,1);
else
    od = zeros(1,2*l*r);
end
od(1:r:(2*l)*r) = [2*idata(nn)-(idata((nn-1):-1:(nn-2*l)))];
od = filter(b,1,od,zf);
odata(nn*r-l*r+1:nn*r) = od(1:l*r);
