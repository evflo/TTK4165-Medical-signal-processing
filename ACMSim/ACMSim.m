
%% Autocorrelation method simulation
%  Authors: Lasse Løvstakken, Hans Torp
%  Last update: 
%  Calls: sim1dSignal.m (simulate signals), acormodFreq.m (generate model)

%% Settings
    % Add / remove signal components
    addNoise = 1;
    addClutter = 1;    
    addBlood = 1;
    % Turn clutter filtering on / off
    clutterFilter = 1;    
    % Number of signal realizations in simulation
    NVects = 1000;
    % Amount of averaged signal realizations, mimicks spatial averaging assuming independent realizations
    NAvg = 9;    
    % Velocity and angle of blood and clutter [m/s] / [rad]
    vBlood = 20/100;
    angBlood = 15*pi/180;
    vClutter = 1/100;
    angClutter = 90*pi/180;
    % Signal-to-noise [dB]
    SNRBlood = 25;
    SNRClutter = 40;
    % Ultrasound imaging setup
    centerFreq = 5e6;       % Center frequency
    pulsePeriods = 4.5;     % Number of periods in pulse
    packetSize = 12;        % Number of pulse emissions per estimate
    FNTransmit = 2.5;       % F-number on transmit
    FNReceive = 1.4;        % F-number on receive
    % Plot distribution of velocity estimates
    plotDist = 1;
    % Plot ensemble average
    plotEns = 1;

%% Setup up simulation    
    % Imaging system properties
    p1 = struct();
    p1.f0 = centerFreq;                                % Center frequency of pulse   
    p1.NPulse = pulsePeriods;                            % Number of periods in pulse
    p1.PRF = 4e3;                               % Pulse repitition frequency (temporal sampling frequency)    
    p1.c = 1540;                                % Speed of sound in tissue
    p1.frs = 1e6;                               % Radial sampling frequency (not used in these simulations, ignore) 
    p1.packetsize = packetSize;                         % Packetsize
    p1.FNTransmit = FNTransmit;                        % F-number on transmit
    p1.FNReceive = FNReceive;                         % F-number on receive
    p1.FNTwoway = 1 / (1/p1.FNTransmit + 1/p1.FNReceive); % Two-way F-number
    p1.lambda = p1.c / p1.f0;                   % Wavelength of emitted pulse
    p1.pulseLength = p1.NPulse*p1.lambda;       % Length of emitted pulse
    p1.beamWidth = p1.FNTwoway*p1.lambda;       % Beam width (rayleigh)
    p1.vNyquist = p1.c*p1.PRF/(4*p1.f0);     

    % Simulation settings
    p1.rlags = 0;                               % Number of lags in radial direction
    p1.tlags = 512;                             % Number of lags in temporal direction
    p1.NVects = NVects;
    p1.Nfft = p1.packetsize*p1.NVects;

    % Blood signal
    p1.BNR1 = SNRBlood;                               % Blood signal to noise ratio
    p1.vBlood1 = vBlood;                        % Blood velocity [m/s]
    p1.angBlood1 = 0*pi/180;                   % Angle of blood flow compared to beam         

    % Clutter signal
    p1.CNR1 = SNRClutter;                               % Clutter signal to noise ratio
    p1.vClutter1 = vClutter / 100;                     % Clutter velocity [m/s]
    p1.angClutter1 = 90*pi/180;                 % Angle of clutter movement compared to beam               
        
%% Generate signals and models
    m = 0:p1.tlags;
    n = -p1.rlags:p1.rlags;             % Not really used in this model, ignore
    Nfft = p1.packetsize*p1.NVects;
    w = linspace(-pi,pi,Nfft);

    % Blood flow signal correlation
    [Rbfunc1, Gb1]=acormodFreq(n,m,w,p1.vBlood1,p1.angBlood1,p1.pulseLength,p1.beamWidth,p1.f0,p1.frs,p1.PRF,p1.BNR1);
    % Clutter signal correlation
    [Rcfunc1, Gc1]=acormodFreq(n,m,w,p1.vClutter1,p1.angClutter1,p1.pulseLength,p1.beamWidth,p1.f0,p1.frs,p1.PRF,p1.CNR1);
    % Thermal noise correlation function, white noise => delta correlated
    [Nr,Nt]=size(Rcfunc1);
    n0=(Nr+1)/2;m0=(Nt+1)/2;
    Rnfunc1=zeros(length(n),length(m));
    Rnfunc1(n0,m0)=1;

    % Theoretical spectrum
    maxVal1 = max(max(Gb1),max(Gc1));    
    Gbnorm1 = Gb1 / maxVal1;
    GbLog1 = 10*log10(Gbnorm1 + realmin);
    Gcnorm1 = Gc1 / maxVal1;
    GcLog1 = 10*log10(Gcnorm1 + realmin);
    Gn1 = abs(fftshift(fft(Rnfunc1,Nfft)));
    Gnnorm1 = Gn1 / maxVal1;
    GnLog1 = 10*log10(Gnnorm1 + realmin);      
    
    % Clutter filter
    hHp =[-0.4909    0.6436   0.0773   -0.1226    -0.1077];% min fase orden 4 FIR high pass
    H = fftshift(abs(fft(hHp,p1.Nfft)));
    HLog = 20*log10(abs(H)/max(abs(H)));
    
%% Do simulations

    % Set up total correlation function
    Rxfunc = zeros(size(Rcfunc1));
    if(addClutter), Rxfunc = Rxfunc + Rcfunc1; end
    if(addNoise), Rxfunc = Rxfunc + Rnfunc1; end
    if(addBlood), Rxfunc = Rxfunc + Rbfunc1; end

    % Simulate NVect signal realizations => sx = [packet size x NVects]
    sx = sim1Dsignal(Rxfunc,p1.packetsize,p1.NVects);    

    % Clutter filter
    if(clutterFilter)
        sx = filter(hHp,1,sx);
        % Throw away FIR initialization samples
        sx = sx((length(hHp)-1):end,:,:);
    end
    
    % Estimate correlation functions
    R0 = mean(conj(sx).*sx);
    R1 = mean(conj(sx(1:end-1,:,:)).*sx(2:end,:,:));

    % "Spatial" average (independent realizations)
    b = ones(1,NAvg)/NAvg;
    R1 = filter2(b,R1);
    R1 = R1(round(NAvg/2):NAvg:end); % Extract independent and averaged signal realizations
    R0 = filter2(b,R0);
    R0 = R0(round(NAvg/2):NAvg:end); % Extract independent and averaged signal realizations
    
    % Estimate mean frequency
    P = ((R0/max(R0))/p1.Nfft);
    f1 = angle(R1)/(2*pi);
        
    % Estimate ensemble average
    P_ens = mean(P,2);
    f1_ens = mean(f1,2);
    
    % Estimate bias and standard deviation
    f_ref = 2*p1.f0*p1.vBlood1*cos(p1.angBlood1)/(p1.c*p1.PRF);
    bias_f1 = f1_ens-f_ref
    std_f1 = std(f1)

%% Plot results
    lineWidth = 2;
    faxis = w/(2*pi);
    hfig = figure(200);
    clf;
    ylim([-80 5]);
    set(gca,'FontSize',12);
    hold on;
    if(plotDist)
        hPDist = plot(repmat(get(gca,'XLim')',1,length(P)),10*log10([P;P]),':','Color',[0.8 0.8 0.8]);
        hFDist = plot([f1;f1],repmat(get(gca,'YLim')',1,length(f1)),':','Color',[0.6 0.6 0.6]);
    end
    if(plotEns)
        hFEns = plot([f1_ens f1_ens],get(gca,'YLim'),'--','LineWidth',lineWidth);
        hPEns = plot(get(gca,'XLim'),10*log10([P_ens P_ens]),'--g','LineWidth',lineWidth);
    end
    if(addClutter)
        hC = plot(faxis,GcLog1,'-b','LineWidth',lineWidth);
    end
    if(addBlood)
        hB = plot(faxis,GbLog1,'-r','LineWidth',lineWidth);
    end
    if(addNoise)
        hN = plot(faxis,GnLog1,'-c','LineWidth',lineWidth);    
    end
    if(clutterFilter)
        hClut = plot(faxis,HLog,'k-','LineWidth',lineWidth);
    end
    hold off;
    grid on; box on;
    xlim([0 0.5]);
    xlabel('Normalized Doppler frequency','FontSize',12);
    ylabel('Normalized power [dB]','FontSize',12);
    title('Theoretical Doppler spectrum','FontSize',14);
    legStr = {}; legH = [];
    if(addClutter),legStr{end+1} ='Clutter'; legH(end+1) = hC; end
    if(addBlood),legStr{end+1} = 'Blood'; legH(end+1) = hB; end
    if(addNoise),legStr{end+1} =  'Noise floor'; legH(end+1) = hN; end
    if(clutterFilter),legStr{end+1} = 'Clutter filter'; legH(end+1) = hClut; end
    if(plotEns)
        legStr{end+1} = 'ACM f ens. avg.'; legH(end+1) = hFEns;
        legStr{end+1} = 'ACM P ens. avg.'; legH(end+1) = hPEns;
    end
    if(plotDist)
        legStr{end+1} = 'ACM f dist.'; legH(end+1) = hFDist(1);
        legStr{end+1} = 'ACM P dist.'; legH(end+1) = hPDist(1);
    end
    legend(legH, legStr,'Location','NorthEast');
 