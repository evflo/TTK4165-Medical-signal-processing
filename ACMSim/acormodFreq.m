function [R,G]=acormodFreq(n,m,w,v,fi,L,D,f0,frs,prf,SNR);
%function [R,G,s1,ro,s2]=acormodFreq(n,m,w,v,fi,L,D,f0,frs,prf);
%	Autocorrelation function model for uniform velocity field
% [R,s1,ro,s2]=acormod(rlags,tlags,v,fi,L,D,f0,frs,prf);
% 6/2-97 Hans Torp
% 2/12-97 returns parameters s1,ro,s2  H.T
% disp('v.1.1');

% Comments by Lasse Løvstakken 04.10.2004:
% Signal model taken from: [1] "Autocorrelation techniques in color flow imaging:
% signal model and statistical properties of the autocorrelation estimates"
% H. Torp et Al, IEEE UFFC vol.41, no.5, 1994.
% 
% Assumtions: 
% 1. Flow region with constant recilinear velocity in space and time.
% 2. The scattering cross-section per unit volume blood (upsilon) is constant.
% 3. The (two-way) beam profile (B) is constant along the pulse length.
% 4. The transversal variation of the (two-way) beam profile (B) can be separated
%    into lateral and elevation components.
% 5. Constant angular beam velocity (not relevant for arrays).

% Just a small number to help prevent division by zero later on
%eps=1e-6;

% Number of lags along and between beams
n1=ones(1,length(n));
m1=ones(1,length(m));

% Radial and lateral velocity
vr=v*cos(fi);
vl=v*sin(fi);

% Sound velocity
c=1540;

% Distance travelled between radial samples (pulse-echo) (Big delta in
% [1]).
d=c/frs/2;

% Correlation length in radial direction (sigma1 in article [1]). Eq.16 in
% [1].
s1=L/sqrt(3)/d;

% Pulse repitition time
T=1/prf;

% Transmit time in radial and lateral direction as a fraction of the PRT.
% Eq.14 in [1].
Tr=L/(vr+eps)/T;
Tl=D/(vl+eps)/T;

% Crosscorrelation coefficient, elevation movement neglected. Interpretation from [1]: Shows to what
% extent the same fluid elements remain inside the volume insonified by the
% ultrasonic beam. If for example the lateral velocity follows the beam
% movement, the same blood cells will be observed in several range cells.
% Eq.16 in [1]
ro=1/sqrt(1+(Tr/Tl)^2);

% Correlation length in lateral direction (sigma2 in article [1]). Eq.16 in
% [1].
s2=ro*Tr/sqrt(3);

% Exponent in Gaussian expression (Covariance matrix times lags) for all
% lags. Eq.15 in [1]
Q=(n'/s1).^2*m1+2*ro*(n'/s1)*(m/s2)+n1'*(m/s2).^2;

% Final correlation function. R(0,0) = 1. Eq.15 in [1].
wd = 4*pi*f0*vr/c*T;
R=10^(SNR/10)*exp(-Q/2).*exp(i*wd*(n1'*m));
% Calculate power spectrum
G = 10^(SNR/10)*sqrt(2*pi)*s2*exp(-0.5*(w-wd).^2*s2^2); % Assume continuous variable => fourier transform pair => analytical expression

