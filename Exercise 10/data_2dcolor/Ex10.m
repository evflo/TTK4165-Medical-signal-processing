%% Exercise 10 TTK4165 Medical Signal Processing
% *Even Florenes Spring 2016*

%% Documentation
% Ex10.m
%
% Purpose: Script answering tasks given in exercise 10 in TTK4165
%
% Made by:
%   Even Florenes NTNU 2016
%
% Related files:
%   imagelog.m: Image a matrix of ultrasound power in log scale
%   clutterfilterrespons.m
%   regressionfilter.m: Removes clutter from iq signal
%
% Last changes:
%   2016-04-06 EF: First attempt
%
% Status:
%   In production

%% Exercise 1: Visualization

load carotid.mat

iq2d = squeeze(iq(6,:,:,1));

gain = -33;
dyn = 40;
imagelog(abs(iq2d).^2,gain,dyn);

%% Exercise 2: Blood tissue ratio

% Set data-points for vessel wall & artery.
r_tissue=48; b_tissue=30;
r_artery=111;b_artery=40;
% Averaging R samples radially and B samples laterally 
R=6; B=3;
% Extract data
iq_tissue =iq(:,r_tissue+[0:R],b_tissue+[0:B],1);
iq_artery = iq(:,r_artery+[0:R],b_artery+[0:B],1);
% Calculate intensity.
I_tissue = abs(iq_tissue).^2;
I_artery = abs(iq_artery).^2;
% Average in all three dimensions.
I_tissue = mean(I_tissue(:));
I_artery = mean(I_artery(:));

intensityRatio = 10*log(I_artery/I_tissue);

fprintf('Intensity ratio between artery and tissue: %0.2f dB\n',intensityRatio);
% How much larger is the intensity in the vessel wall compared to the
% vessel?
% 
% The intensity, power transferred per unit area, is much larger in the 
% vessel wall.
%
%

%% Exercise 3: Regression filter for removing clutter

% Filter data with current order N. 
%y=regressionfilter(squeeze(iq(:,:,:,frameno),N));
frameNo=1;
packetNo=6;
gain = -20;
for N=-1:4
% Create subplot.
figure(2),subplot(2,3,N+2);
% Filter data with current order N. 
y=regressionfilter(iq(:,:,:,frameNo),N);
% Squeeze to remove all other packet-data. 
y=squeeze(y(packetNo,:,:));
% Show image.
imagelog(abs(y).^2,gain,dyn);
% Set title.
title(['N=',num2str(N)]);
end
% Comment on the remaining blood and tissue signal for increasing filter
% order:
%
% As we increase N the image becomes darker and darker. For higher N the
% only remaining part in the image is the fluctuation in the artery and
% vein. Every part of the surrounding vessel wall is removed.
%
%
packetSize=13;
% Set data-points for vessel wall & artery.
r_tissue=48; b_tissue=30;
r_artery=111;b_artery=40;
% Averaging R samples radially and B samples laterally 
R=6; B=3;
% Extract data
N = -1:packetSize-2;
iq_tissue =iq(:,r_tissue+[0:R],b_tissue+[0:B],1);
iq_artery = iq(:,r_artery+[0:R],b_artery+[0:B],1);
I_ratio = zeros(1,length(N));
I_tissueMean = zeros(1,length(N));
I_arteryMean = zeros(1,length(N));
for i = 1:length(N)
    y_tissue = regressionfilter(iq_tissue,N(i));
    y_artery = regressionfilter(iq_artery,N(i));
    I_tissue = abs(y_tissue).^2;
    I_artery = abs(y_artery).^2;
    % Average in all three dimensions.
    I_tissueMean(i) = mean(I_tissue(:));
    I_arteryMean(i) = mean(I_artery(:));
    I_ratio(i) = 10*log(I_arteryMean(i)/I_tissueMean(i));
end
figure(3),plot(N,10*log(I_tissueMean)),title('Intensity of vessel wall');
xlabel('N');ylabel('dB');
intensityClutterComponent = abs(10*log(I_arteryMean(2))-10*log(I_arteryMean(1)));
fprintf('The intensity of the clutter component compared to the \n')
fprintf('echo strength from the red blood cells is: %0.2f dB\n',intensityClutterComponent);

figure(4),plot(N,I_ratio),title('Intensity ratio of artery on tissue');
xlabel('N');ylabel('dB');
%
%
%
%
%
%
figure(5);hold on;
for N=0:4
    [y,c,Fm]=regressionfilter(iq(:,:,:,frameNo),N); 
    [f,H]=clutterfilterrespons(Fm); 
    plot(f,10*log10(H),'k-');
end;
title('Clutter filter response'); xlabel('Normalized Doppler frequency'); ylabel('dB');
ylim([-20 3]); xlim([-0.5 0.5]);
%
%
%
%
%
%% Exercise 4: Complex plots
frameNo = 1;
r = 48;
b = 30;
iq1 = iq(:,:,:,frameNo);
% Plot the complex data. Circle the datapoints. 
figure(6);

subplot(3,2,1);
plot(squeeze(iq1(:,r,b)),'k-o');
% Make the axes scale equally in each direction. 
axis('equal');
% Label the axes.
xlabel('Re');ylabel('Im');
% Save the current axes for later. 
axlimits=axis;
% Turn grid on.
grid on;
% Use the relation between phase shift and movement 
% to explain how the signal looks.
%
%
%
%
%
for N=0:4,
    subplot(3,2,N+2);
    % Filter... 
    [y,c]=regressionfilter(iq1,N); 
    plot(squeeze(c(:,r,b)),'kx'); % Set axis.
    axis(axlimits)
    grid on;
    % Set title. 
    title(['N=',num2str(N)]); 
end;
%How does the clutter component look?
%
%
%
%
%
r = 111;
b= 40;
figure(7)
subplot(2,1,1),plot(squeeze(y(:,r,b)),'k-o');
title('Middle of artery after filtration');
axis('equal'); % Make the axes scale equally in each direction. 
xlabel('Re');
ylabel('Im'); % Label the axes.
%Is the adaption between signal and clutter polynomial as good in blood as in tissue?
%
%
%
%
subplot(2,1,2),plot(squeeze(c(:,r,b)),'k-o');
title('Middle of artery residual');
%What does the residual signal look like?
%
%
%
%

%% Exercise 5: Velocity estimation using the autocorrelation method