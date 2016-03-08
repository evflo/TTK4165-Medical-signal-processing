function echo2dsim(nodefaults)

%ECHO2DSIM  2D ultrasound image simulation
%
% Object (source) as 8 bit grayscale image from image-file
%
% Covolution in 2D frequency domain, assuming constant F-number
% x is along aperture (azimut)
% z is in ultrasound beam direction (range)
%
% Read an object image from file, set parameters, and look at the result
%


if ~nargin, nodefaults=0;end;

if nodefaults,
    if exist('echo2dsim_default.mat')==2,
        delete('echo2dsim_default.mat');
    end;
end;


% Ultrasound probe and pulse parameters
    f0=8e6;     % Pulse center freq.
    B=3.2e6;      % Pulse bandwidth (50%)

    dBnoise=-55;% White noise level

    %Transmit aperture and apodization:
    FnTx=2;     % F-number, FnTx = distance/aperture diameter
    apodTx=0;   % Apodisation function Tx. 
                % Kaiser window with beta=apodTx.  
                % apodTx=0 gives rectangular window

    %Receive aperture and apodization:
    FnRx=2;     % F-number receive aperture
    apodRx=0;   % Apodisation function Rx 
                % Kaiser window with beta=apodRx.  
                % apodTx=0 gives rectangular window

    % Object parameters
    Ofname='artery.tif'; % Name of object file. 8 bit image with 256 gray levels     
    dxo=0.015e-3;  % Pixel size in x direction (azimut)
    dzo=0.015e-3;  % Pixel size in z direction (range)

    % General parameters
    c=1540;              % Speed of sound
    fdepAttenuation=4.0; % Frequency dependent attenuation [dB/MHz]

    fmax=20e6;     % Ultrasound frequency upper limit. Determines image sampling rate
    lmin=c/fmax;   % Minimum wavelength 
                       % Constraint: dxo<lmin/2 and dzo<lmin/4

    % Display parameters
    gain=54;  % Gain  
    dyn=40;   % Dynamic range

    if exist('echo2dsim_default.mat'),
        disp('Loading previous parameters from file: echo2dsim_default.mat');
        disp('Use >>echo2dsim(1) to return to initial settings');
        disp(' ');
        load('echo2dsim_default.mat','f0','B','dBnoise','FnTx','FnRx',...
             'apodRx','Ofname','dxo','fdepAttenuation','gain');
    end;

    questions={'Frequency [MHz] (1 10)',...
               'Dual-sided, -6dB bandwidth (10-100%)',...
               'F-number Transmit (1-4)',...
               'F-number Receive (1-4)',...
               'Apodization Receive (Kaiser window, 0=rectangular)',...
               'Object file (blank = browse)',...
               'Pixel size [mm]',...
               'Noise level [dB]',...
               'Frequency dependent attenuation [dB/MHz]',...
               'Gain [dB]'};

    defaults={[num2str(f0/1e6,3)],...
              [num2str(B/f0*100,2)],...
              [num2str(FnTx,2)],...
              [num2str(FnRx,2)],...
              [num2str(apodRx)],...
              [Ofname],...
              [num2str(dxo*1e3)],...
              [num2str(dBnoise)],...
              [num2str(fdepAttenuation)],...
              [num2str(gain)]};
    answer=inputdlg(questions,'Parameters',1,defaults);

    if prod(size(answer))==0,
        return;
    end;

    f0=str2num(char(answer(1)))*1e6;
    f0=max(1e6,min(10e6,f0));
    B=str2num(char(answer(2)));
    B=min(100,max(10,B));
    B=B/100*f0;
    FnTx=str2num(char(answer(3)));
    FnTx=min(4,max(1,FnTx));
    FnRx=str2num(char(answer(4)));
    FnRx=min(4,max(1,FnRx));
    apodRx=round(str2num(char(answer(5))));
    Ofname=char(answer(6));
    dxo=str2num(char(answer(7)))/1e3;
    dzo=dxo;
    dBnoise=str2num(char(answer(8)));
    fdepAttenuation=str2num(char(answer(9)));
    gain=str2num(char(answer(10)));

    % Load object image
    eval('Orgim=imread(Ofname);err=0;','err=1;');
    while err,
        disp('Use browser window to select object file');
        [fname,pname]=uigetfile('*.*','Open object file');
        Ofname=[pname,fname];
        eval('Orgim=imread(Ofname);err=0;','err=1;');
    end;        
    if ndims(Orgim)==3,
        disp('Converting to grayscale');
        Orgim=mean(double(Orgim),3);
    end;

    Orgim=double(Orgim);    
    disp('Stretching image to 0-255');
    Orgim=(Orgim-min(Orgim(:)));
    Orgim=255*Orgim/max(Orgim(:));
    Orgim=uint8(round(Orgim));
    
    save('echo2dsim_default.mat','f0','B','dBnoise','FnTx','FnRx',...
         'apodRx','Ofname','dxo','fdepAttenuation','gain');

    if (length(Ofname)>20),
        Ofname=['..',Ofname(end-19:end)];
    end;
    
    %Tex compatibility
    Ofname=strrep(Ofname,'\','\\');
    Ofname=strrep(Ofname,'_','\_');
    
    parameters={['  Object file:     '],...
                ['    ',Ofname],...
                ['  Pixel size:      '],...
                ['    ',num2str(dxo*1e3),' [mm]'],...
                ['  Frequency:       '],...
                ['    ',num2str(f0/1e6),' [MHz]'],...
                ['  Bandwidth:       '],...
                ['    ',num2str(B/1e6),' [MHz] (',num2str(100*B/f0,2),'%)'],...
                ['  F-number Tx:     '],...
                ['    ',num2str(FnTx)],...
                ['  F-number Rx:     '],...
                ['    ',num2str(FnRx)],...
                ['  Apodization Rx:  '],...
                ['    ',num2str(apodRx)],...
                ['  Noise level:     '],...
                ['    ',num2str(dBnoise),' [dB]'],...
                ['  Freq. dep. att.: '],...
                ['    ',num2str(fdepAttenuation),' [dB/MHz]'],...
                ['  Gain:            '],...
                ['    ',num2str(gain),' [dB]']};


%Delete all exisiting waitbars and initialise new one:
delete(findobj(allchild(0),'tag','TMWWaitbar'));
waitbarhandle=waitbar(0,'Please wait');

% Create object
    createobject=1;
    %Do not re-create object if object is identical to previous object:
    if exist('echo2dsim_object_temp.mat')==2,
        load('echo2dsim_object_temp.mat','thisobject','thisdxo');
        if isequal(thisobject,Orgim)&isequal(thisdxo,dxo),
            disp('Re-using object');
            load('echo2dsim_object_temp.mat');
            createobject=0;
        end;
    end;
    
    if createobject,
        waitbar(1/3,waitbarhandle);

        % Dimensions of input image
        zo=(0:size(Orgim,1)-1)*dzo;
        xo=(0:size(Orgim,2)-1)*dxo;

        % Generate k-space image of object
        O=kspaceobject(Orgim,dzo,dxo,lmin);

        %return if O is empty
        if size(O,1)==0,
            return;
        end;

        [Nzi,Nxi]=size(O);
        dzi=zo(end)/Nzi;
        dxi=xo(end)/Nxi;
        zi=[0:(Nzi-1)]*dzi;
        xi=[0:(Nxi-1)]*dxi;

        thisobject=Orgim;
        thisdxo=dxo;   
        save('echo2dsim_object_temp.mat','thisobject','thisdxo',...
             'xo','zo','O','Nzi','Nxi','dzi','dxi','zi','xi');
    end;

    clear('thisobject');
    waitbar(2/3,waitbarhandle);

    % Create aperture frequency response A
    NaTx= round(Nxi/dzi/FnTx*dxi/2); %Na/Nxi*1/dxi*fn=1/dzi
    NaRx= round(Nxi/dzi/FnRx*dxi/2);
    ap=conv(kaiser(NaTx,apodTx),kaiser(NaRx,apodRx));
    A=kspaceaperture(Nzi,Nxi,ap);

    waitbar(7/9,waitbarhandle);

    % Create pulse frequency response H
    f=(0:Nzi-1)'/Nzi*c/2/dzi; % Frequency axis
    Hsc=(f/f(end)).^2;

    % Frequency dependent attenuation response
    Hatt=10.^(-fdepAttenuation*(f*1e-6)/20);

    % Pulse frequency response
    Brms=B/2/sqrt(2*log(2));%relation between -6dB bandwidth B and rms bandwidth
    Hp=exp(-0.5*((f-f0)/Brms).^2);

    %Hp=exp(-0.5*((f-f0)/(B/4)).^2);%Tja..............

    %Combined frequency response
    H=Hp.*Hatt.*Hsc;

    waitbar(8/9,waitbarhandle);

    % Matrix with pulse frequency response as columns
    H=H*ones(1,Nxi); 

    % Create point spread function in the frequency domain by
    % combining Aperture frequency response and pulse frequency response
    Psf=A.*H;
    % Normalize:
    Psf=Psf/sqrt(sum(sum(abs(Psf).^2)));

    %Reshape H as column vector
    H=H(:,1);

    %Create Ultrasound image: 
    % Multiply the pointspread function PSF and the object O in the frequency domain
    U=Psf.*O;
    % Use Inverse Fourier Transform to transform back into spatial domain 
    echoIm=ifft2(U);

    %Add white noise:
        %White noise with unit intensity
        wnoise=randn(size(O))+i*randn(size(O));
        %Scale white noise level according to dBnoise:
        wnoise=wnoise*10^(dBnoise/20);
        %Add white noise
        echoIm=echoIm+wnoise;

    waitbar(9/9,waitbarhandle);


%k-space frequency axes:
    kxaxis=[0:size(O,2)-1]/size(O,2);
    kzaxis=[0:size(O,1)-1]/size(O,1)/2;
    faxis=f/1e6;%Frekvensakse i MHz;

%FFTshift along kx axis:
    Nkx=ceil(length(kxaxis)/2);
    flipkxindex=[Nkx+1:length(kxaxis),1:Nkx];
    flipkxaxis=[kxaxis([Nkx+1:length(kxaxis)])-1, kxaxis(1:Nkx)];

%Plot object and resulting image
figure(1);clf;

    subplot(2,3,1);
    set(gca,'fontsize',8);
    image(xo*1e3,zo*1e3,Orgim);
    colormap(gray(256));
    axis('image');
    title('Object image');
    ylabel('Range z [mm]');
    xlabel('Azimut x [mm]');

    subplot(2,3,2);
    set(gca,'fontsize',8);
    Olog=256*logabs(abs(O),-60,50);
    colormap(gray(256)); 
    image(flipkxaxis,faxis,Olog(:,flipkxindex));
    title('Object image, frequency domain');
    xlabel('k_x [fs]');
    ylabel('Frequency [MHz]');
    
    subplot(2,3,4);
    set(gca,'fontsize',8);
    colormap(gray(256));
    image(xi*1e3,zi*1e3,255*logabs(echoIm,gain,dyn));
    axis('image');
    title('Ultrasound image');
    ylabel('Range z [mm]');
    xlabel('Azimut x [mm]');

    %Gain controls:
    uicontrol('style','text','string','Gain','units','normal',...
              'pos',[0.0 0.15 0.05 0.05]);
    uicontrol('style','push','string','+','units','normal',...
             'pos',[0.0 0.10 0.05 0.05],...
             'callback',...
             ['gain=gain+3;figure(1); subplot(2,3,4);',...
              'image(xi*1e3,zi*1e3,255*logabs(echoIm,gain,dyn));',...
              'colormap(gray(256));axis(''image'');',...
              'title(''Ultrasound image'');',...
              'ylabel(''Range z [mm]'');xlabel(''Azimut x [mm]'');',...
              'save(''echo2dsim_default.mat'',''gain'',''-append'');']);
    uicontrol('style','push','string','-','units','normal',...
             'pos',[0.0 0.05 0.05 0.05],...
             'callback',...
             ['gain=gain-3;figure(1); subplot(2,3,4);',...
              'image(xi*1e3,zi*1e3,255*logabs(echoIm,gain,dyn));',...
              'colormap(gray(256));axis(''image'');',...
              'title(''Ultrasound image'');',...
              'ylabel(''Range z [mm]'');xlabel(''Azimut x [mm]'');',...
              'save(''echo2dsim_default.mat'',''gain'',''-append'');']);

    subplot(2,3,5);
    set(gca,'fontsize',8);
    Ulog=255*logabs(abs(U),10,50);
    colormap(gray(256)); 
    imagesc(flipkxaxis,faxis,Ulog(:,flipkxindex));
    title('Ultrasound image, frequency domain');
    xlabel('k_x [fs]');
    ylabel('Frequency [MHz]');


    subplot(1,3,3);
    set(gca,'fontsize',8);
    xlim([0 1]);
    ylim([0 1]);
    set(gca,'xtick',[],'ytick',[]);
    text(0,1,parameters,'fontsize',8,...
         'verticalalignment','top',...
         'horizontalalignment','left');
    %Frame
    hold on;
    plot(xlim,min(ylim)*[1 1],'k-');
    plot(xlim,max(ylim)*[1 1],'k-');
    plot(min(xlim)*[1 1],ylim,'k-');
    plot(max(xlim)*[1 1],ylim,'k-');

    uicontrol('style','push','string','New simulation','callback','echo2dsim',...
              'units','normal','pos',[0.80 0.02 0.18 0.09]);

figure(2);clf;
    subplot(2,3,1);
    set(gca,'fontsize',8);
    Alog=255*logabs(A,15,50);
    colormap(gray(256)); 
    image(flipkxaxis,faxis,Alog(:,flipkxindex));
    title('Aperture frequency response');
    xlabel('k_x [fs]');
    ylabel('Frequency [MHz]');

    subplot(2,3,2);
    set(gca,'fontsize',8);
    Hlog=255*logabs(H,80,50);
    colormap(gray(256)); 
    image(flipkxaxis,faxis,Hlog*ones(1,length(flipkxaxis)));
    title('Pulse frequency response');
    xlabel('k_x [fs]');
    ylabel('Frequency [MHz]');

    subplot(2,3,4);
    set(gca,'fontsize',8);
    Psflog=255*logabs(Psf,80,50);
    colormap(gray(256)); 
    image(flipkxaxis,faxis,Psflog(:,flipkxindex));
    title('Point spread function, frequency domain');
    xlabel('k_x [fs]');
    ylabel('Frequency [MHz]');


    %Plot point spread function in spatial domain
    psf=fftshift(ifft2([Psf]));

    subplot(2,3,5);
    set(gca,'fontsize',8);
    image(xi*1e3,zi*1e3,255*logabs(psf,140,80));
    colormap(gray(256));
    title('PSF, spatial response');
    xlabel('Azimut x [mm]');
    ylabel('Range z [mm]');
    axis image;


    %Frequency and time response of pulse
    subplot(2,3,3);
    set(gca,'fontsize',8);
    p=20*log10(abs(H)+1e-30);%Add small nonzero number to avoid log og zero
    p=p-max(p);
    plot(p,f/1e6,'k-');
    title('Pulse frequency response');
    xlim([-60 0]);
    ylim([0 20]);
    xlabel('[dB]');
    ylabel('Frequency [MHz]');
    set(gca,'ydir','reverse');
    set(gca,'ylim',[min(f/1e6) max(f/1e6)]);
    
    subplot(2,3,6);
    set(gca,'fontsize',8);
    p=real(fftshift(ifft([H;zeros(size(H))])));
    z=interp(zi,2)*1e3;

    p=p./max(abs(p));
    plot(p,z,'k-');
    title('Pulse response');
    ylabel('z [mm]');
    xlabel('Normalized amplitude');
    xlim([-1.1 1.1]);
    set(gca,'ydir','reverse');
    set(gca,'ylim',[min(z) max(z)]);


    uicontrol('style','push','string','Zoom +','units','normal',...
             'pos',[0.92 0.10 0.08 0.05],...
             'callback',...
             ['subplot(2,3,6);set(gca,''ylim'',mean(ylim)+diff(ylim)*[-0.25 0.25]);']);
    uicontrol('style','push','string','Zoom -','units','normal',...
             'pos',[0.92 0.05 0.08 0.05],...
             'callback',...
             ['subplot(2,3,6);if(min(ylim)<=0), return;end;',...
              'set(gca,''ylim'',mean(ylim)+diff(ylim)*[-1 1]);']);


    %Transfer some variables to main workspace
    save('echo2dsimtemp.mat','gain','xi','zi','echoIm','dyn');
    evalin('base','load(''echo2dsimtemp.mat'');');
    delete('echo2dsimtemp.mat');

delete(waitbarhandle);


%===================%
%   Sub functions   %
%===================%


function A=kspaceaperture(Nz,Nx,a);
% Generate frequency response of the aperture 
    A=zeros(Nz,Nx);
    a=a(round((length(a)+1)/2):end);
    Na=length(a);
    x=0:Na-1;
    for z=2:Nz-1,
       xi=x*(Nz-1)/z;
       xi=xi(find(xi<(Na-1)));
       w=interp1(x,a,xi);
       A(z+1,1:length(xi))= w;
       A(z+1,end-length(xi)+2:end)=w(end:-1:2);
    end;



function O=kspaceobject(im,dzo,dxo,lmin);
% Generate image of object in k-space (frequency domain)
%
% FFT2, lowpass filter and decimation 
%
% im   - spatial image of object 
% dzo  - pixel size in z-direction of object image [m]
% dxo  - pixel size in x-direction of object image [m]
% lmin - minimum wavelength 
% 

    c=1540;
    [Nzo,Nxo]=size(im);
    fso=c/2/dzo;        % Sampling frequency in range direction
    dfso=(c/2/dzo)/Nzo; % Frequency increment
    
    %decimate to dxi~=lmin/2 and dzi~=lmin/2
    Nzi=round(Nzo*dzo/(lmin/2));
    Nxi=round(Nxo*dxo/(lmin/2)); 
    Nxi2=round((Nxi-1)/2);
    
    if Nzi>Nzo/2,
        disp('Pixel size too large');
        return;
    end;

    if Nxi>Nxo/2,
        disp('Pixel size too large');
        return;
    end;

    %Convert image to double precision
    im=double(im);
    
    %2D FFT:
    O=fft2(im);

    %Decimation
    O=O(1:Nzi,[1:Nxi2,Nxo-Nxi2+2:Nxo]);


function w = kaiser(n,beta)
%KAISER Kaiser window.
%   W = KAISER(N,beta) returns the BETA-valued N-point Kaiser window.

    nw = round(n);
    bes = abs(besseli(0,beta));
    odd = rem(nw,2);
    xind = (nw-1)^2;
    n = fix((nw+1)/2);
    xi = (0:n-1) + .5*(1-odd);
    xi = 4*xi.^2;
    w = besseli(0,beta*sqrt(1-xi/xind))/bes;
    w = abs([w(n:-1:odd+1) w])';



function y=logabs(x,gain,dyn);
%LOGABS   Logaritmisk kompresjon av absoluttverdi

    y=abs(x); %Ta absoluttverdi 
    y=y+1e-30; %Legg til eit lite positivt tal for å unngå logaritme av 0
    y=20*log10(y); %Gjer amplitude om til dB
    y=y+gain; %Legg til gain (i dB)
    y=y/dyn;  %Skaler dynamisk område
    y=max(0,y); %Sett negative verdiar lik 0
    y=min(1,y); %sett verdiar >1  lik 1
    
    
function [odata,b] = interp(idata,r)
%INTERP Resample data at a higher rate using lowpass interpolation.
%   Y = INTERP(X,R) resamples the sequence in vector X at R times
%   the original sample rate.  The resulting resampled vector Y is
%   R times longer, LENGTH(Y) = R*LENGTH(X).

l = 4;
alpha = .5;

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
    