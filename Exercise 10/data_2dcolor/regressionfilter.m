function [y,c,Fm]=regressionfilter(x,N);

%Regression filter for ColorFlow
%
% function [y,c,Fm]=regressionfilter(x,N);
%
% Inputs: x  - 3D matrix with IQ-data 
%         N  - Polynomial order
%
% Output: y  - IQ-signal after filtrering y=x-c
%         c  - clutter komponent, c=x-y
%         Fm - Filter matrix
%

%Johan Kirkhorn, 22mar01


    if nargin<2,
        error('Too few input arguments');
    end;
    if ndims(x)~=3,
        error('IQ-data matrix must be 3-dimensional');
    end;

    [packetsize,samples,beams]=size(x);

    if N~=round(N),
        error('N must be an integer');
    end;

    if N>packetsize-2,
        error('N too large (max. packetsize-2)');
    end;

    if N<0,
        c=zeros(size(x));
        y=x;
        return;
    end;

    Fm=polyreg(packetsize,N,0);

    y=zeros(size(x));

    for beamno=1:beams,
        y(:,:,beamno)=(squeeze(x(:,:,beamno)).'*Fm).';
    end;

    c=x-y;



function Fm=polyreg(N,P,dw,skipendpoints,AttenAtZero)

%POLYREG    Generate NxN filter matrix Fm for polynomial regression filter
%
% Fm=polyreg(N,P,dw,skipendpoints,AttenAtZero);
% N = packet size (number of temporal samples)  N >= 2
% P = max. polynomial order, -1<= P < N-1.  
%            P=-1 gives Fm=identity matrix
% dw [rad] = filter center frequency offset. 
%            If abs(dw)>0, an aditional zero is added at freq.=0
%            In this case P should be less than N-2
% skipendpoints=1 removes first and last output sample in the filter output signal
%            Default: skipendpoints=0
% AttenAtZero is the attenuation at zero frequency in [dB]
% If not spesified AttenAtZero= + inf.


% Written by  : Hans Torp 29/2-96
% Version no. : 1.1  second version
%               1.2 opprydning, virker naa ogsaa for P>6
%               1.3 skipendpoints inkludert 16/1 97 Hans Torp
%               1.4 AttenAtZero inkludert 10/3 97 Hans Torp

if nargin<3,dw=0;end;
if nargin<4,skipendpoints=0;end;
if nargin<5,AttenAtZero=  120;end;

if abs(dw)>0.0001, P=P+1; ,else dw=0;end;
if P> (N-2), disp('Filter order too large');return; end;

Fm= zeros(N,N);

t=1:N;
if P>-1,
%form the polynomial basis vectors by Gramm-Schmith
base=zeros(P+1,N);
for p=0:P,
    b=t.^p;%polynom of order p
    if p>0, if abs(dw)>0, b=t.^(p-1).*exp(i*dw*t);end;end;% polynom of order p-1, incl. mixing
    for k=0:p-1, %ortogonalization
             bk=base(k+1,:);
             b=b-(b*bk')*bk;
    end; 
    b=b/sqrt(b*b');%Normalization
    base(p+1,:)=b;
    Fm=Fm+b'*b;
end;%for p
end;%P>-1

c=1;
if  AttenAtZero < 100,
    c = 1- 10^(-AttenAtZero/20);
end;
c2 = c*(2-c);

Fm = eye(N,N) - c*Fm;

if skipendpoints, Fm=Fm(2:N-1,:);end;



    
