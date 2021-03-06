% Ex8.m
% Purpose: 
%   Script to answer tasks in exercise 8 given in course
%   TTK4165 Medical Signal Processing
%
% Made by:
%   Even Florenes NTNU 2016
%
% Last changes:
%   2016-03-08: First attempt
%

%% PART 2
load slowmotion

middleBeamIq = squeeze(iq(:,4,:));

frameRate = s.Framerate_fps; % nFrames/seconds

nFrames = size(middleBeamIq,2); %nFrames

nSamples = size(middleBeamIq,1); 

nSeconds = nFrames/frameRate;   %seconds = (nFrames/(nFrames/seconds))

time = 0:nSeconds/(nFrames-1):nSeconds;

distanceLength = s.iq.DepthIncrementIQ_m;

distance = 0:distanceLength/(nSamples-1):distanceLength;

amplitude = abs(middleBeamIq);
power = amplitude.^2;

gain = -50;
dynamicRange = 50;

logarithimicPower = imagelog(power,gain,dynamicRange);

figure,imagesc(time,10^4*distance,logarithimicPower),colormap(gray),...
    title('M-Mode of IQ-beams'),xlabel('Time[sec]'),ylabel('Distance[cm]');
hold on
%% Subpart 1:Measure parameters
[x,y] = ginput(3);

rotationPeriod = x(2)-x(1); % T [s]

t0 = x(3); 

excenterDistance = (y(3)-mean([y(1),y(2)]))/2; 
%y(1) and y(2) are supposed to
%be the same height, but 
%measurements may misaligne.
%Use mean of the measurements
%to get best possible precision

fprintf('T= %g sec\n',rotationPeriod);
fprintf('t0= %g sec\n',t0);
fprintf('R= %g cm\n',excenterDistance);
%% Subpart 2: Plot point position in M-mode image
pistonLength = 10;
fixedPoint = 4;
pistonAngularFrequency = (2*pi*(time-t0))/rotationPeriod;

%pistonVelocity = pistonVelocityAmplitude*sin(pistonVelocityAngularFrequency);
pistonPosition = pistonLength+excenterDistance*cos(pistonAngularFrequency);
pointPosition = fixedPoint+excenterDistance*cos(pistonAngularFrequency);

plot(time,pointPosition,'y'),hold off;

%% PART 3
pistonVelocityAmplitude = -(2*pi*excenterDistance)/rotationPeriod;
pistonVelocity = pistonVelocityAmplitude*sin(pistonAngularFrequency);

pointVelocity = -pistonVelocity;

figure,plot(time,pointVelocity),title('Point velocity'),xlabel('Time [sec]'),...
    ylabel('Velocity [cm/s]');

%% Subpart 1: Choose frame and calulate parameters
maxVelocity = max(pointVelocity);
maxVariance = (1/30)*maxVelocity;

framesBelowOneThirdVelocity = find(pointVelocity < (1/3)*maxVelocity+maxVariance);

framesWithThirdVelocity = find(pointVelocity(framesBelowOneThirdVelocity) > (1/3)*maxVelocity ...
        - maxVariance);
    
frameWithThirdVelocity = median(framesWithThirdVelocity);

imageFrames = middleBeamIq(:,frameWithThirdVelocity-1:frameWithThirdVelocity);

figure,plot(10^4*distance,imageFrames),title('Two subsequent beams(in time) with 1/3 of max velocity'),...
xlabel('Positon[cm]');

imageFramesVelocity = mean(pointVelocity(frameWithThirdVelocity-1:frameWithThirdVelocity));

radialDisplacement = imageFramesVelocity/frameRate;% [cm]

speedSound = 1540*100; %[cm/s]

timeShiftAnalytically = 2*(radialDisplacement/speedSound);

fprintf('Analytic velocity: %g cm/s\n',imageFramesVelocity);
fprintf('Calculated radial displacement = %g cm \n',radialDisplacement);
fprintf('Calculated time shift = %g nsec\n',10^(9)*timeShiftAnalytically);

samplingFrequencyRF=200e6; %RF sampling freq.

demodulationFrequency=s.iq.fDemodIQ_Hz; %IQ demodulation freq.
samplingFrequencyIQ=s.iq.frsIQ_Hz; %IQ sampling freq

firstRFBeam=iq2rf(imageFrames(:,1),demodulationFrequency,samplingFrequencyIQ,samplingFrequencyRF);
secondRFBeam=iq2rf(imageFrames(:,2),demodulationFrequency,samplingFrequencyIQ,samplingFrequencyRF);

timeRF=[1:length(firstRFBeam)]/samplingFrequencyRF;
timeRF = 10^9*timeRF;
figure,plot(timeRF,firstRFBeam,'b'),title('RF data of two subsequent beams(in time) with 1/3 of max velocity'); hold on;
plot(timeRF,secondRFBeam,'r'); hold off;
ylabel('RF signal');
xlabel('Time[nsec]');

%% PART 4
crossCorrelation = zeros(1,201);
i = 1;
for l = -100:100

    for t = 15000:20000
        crossCorrelation(i) = crossCorrelation(i)+...
            firstRFBeam(t)*secondRFBeam(t+l);
    end
    i = i+1;
end

figure,plot(-100:100,crossCorrelation),title('Cross correlation');

sampleShiftMax = find(crossCorrelation  == max(crossCorrelation));
sampleShiftMax = (sampleShiftMax-101);
timeShiftCalculated = sampleShiftMax/samplingFrequencyRF;
timeShiftAnalytically = 2*abs(pointPosition(frameWithThirdVelocity-1)- ...
    pointPosition(frameWithThirdVelocity))/speedSound;

fprintf('After %d samples (e.g %g nsec) the cross-correlation has its maximum\n',sampleShiftMax,10^9*timeShiftCalculated);
fprintf('Analytical time shift: %g nsec\n',timeShiftAnalytically*10^9);
phaseShiftCalculated = timeShiftCalculated*2*pi*demodulationFrequency;
phaseShiftAnalytically = timeShiftAnalytically*2*pi*demodulationFrequency;

fprintf('Calculated phase shift: %g radians \n',phaseShiftCalculated);
fprintf('Analytical phase shift: %g radians \n',phaseShiftAnalytically);

%% PART 5

phaseShiftAutoCorr = angle(mean(conj(imageFrames(60:80,2)).*imageFrames(60:80,1)));

fprintf('Calculated phase shift using autocorrelation method: %g radians \n',phaseShiftAutoCorr);

firstMiddleBeamIq = middleBeamIq(:,1:nFrames-1);
secondMiddleBeamIq = middleBeamIq(:,2:nFrames);

totalPhaseShift = angle(mean(conj(secondMiddleBeamIq).*firstMiddleBeamIq));

conversionValue = (speedSound*frameRate)/(4*pi*demodulationFrequency);

calculatedVelocity = conversionValue*totalPhaseShift;

figure(9),plot(time(1:630),calculatedVelocity,time,pointVelocity);
title('Estimated velocity using auto correlation method');
xlabel('Time [sec]');
legend('Estimated using the auto correlation method',...
    'Analytically calculated','Location','southeast');

%% PART 6

load fastmotion.mat
middleBeamIq = squeeze(iq(:,4,:));

frameRate = s.Framerate_fps; % nFrames/seconds

nFrames = size(middleBeamIq,2); %nFrames

nSamples = size(middleBeamIq,1); 

nSeconds = nFrames/frameRate;   %seconds = (nFrames/(nFrames/seconds))

time = 0:nSeconds/(nFrames-1):nSeconds;

distanceLength = s.iq.DepthIncrementIQ_m;

distance = 0:distanceLength/(nSamples-1):distanceLength;

amplitude = abs(middleBeamIq);
power = amplitude.^2;

gain = -50;
dynamicRange = 50;

logarithimicPower = imagelog(power,gain,dynamicRange);

figure,imagesc(time,10^4*distance,logarithimicPower),colormap(gray),...
    title('M-Mode of IQ-beams'),xlabel('Time[sec]'),ylabel('Distance[cm]');
hold on
%% Subpart 1:Measure parameters
[x,y] = ginput(3);

rotationPeriod = x(2)-x(1); % T [s]

t0 = x(3); 

excenterDistance = (y(3)-mean([y(1),y(2)]))/2; 
%y(1) and y(2) are supposed to
%be the same height, but 
%measurements may misalign.
%Use mean of the measurements
%to get best possible precision

fprintf('T= %g sec\n',rotationPeriod);
fprintf('t0= %g sec\n',t0);
fprintf('R= %g cm\n',excenterDistance);

%% Subpart 2: Plot position and velocity
pistonLength = 10;
fixedPoint = 4;

pistonAngularFrequency = (2*pi*(time-t0))/rotationPeriod;
pistonPosition = pistonLength+excenterDistance*cos(pistonAngularFrequency);
pointPosition = fixedPoint+excenterDistance*cos(pistonAngularFrequency);

plot(time,pointPosition,'y'),hold off;
pistonVelocityAmplitude = -(2*pi*excenterDistance)/rotationPeriod;
pistonVelocity = pistonVelocityAmplitude*sin(pistonAngularFrequency);

pointVelocity = -pistonVelocity;

figure,plot(time,pointVelocity),title('Point velocity for fastmotion.mat');
xlabel('Time[sec]');

%% Subpart 3: Calculate velocity by auto-correlation method
firstMiddleBeamIq = middleBeamIq(:,1:nFrames-1);
secondMiddleBeamIq = middleBeamIq(:,2:nFrames);
speedSound = 1540;
demodulationFrequency = s.iq.fDemodIQ_Hz;
totalPhaseShift = angle(mean(conj(secondMiddleBeamIq).*firstMiddleBeamIq));

conversionValue = (speedSound*frameRate)/(4*pi*demodulationFrequency);

calculatedVelocity = conversionValue*totalPhaseShift;

figure,plot(time(1:432),calculatedVelocity,time,pointVelocity);
title('Estimated velocity using auto correlation method');
xlabel('Time [sec]');
legend('Estimated using the auto correlation method',...
    'Analytically calculated','Location','southeast');


%% PART 7

load slowmotion.mat

firstIqBeam = squeeze(iq(:,1,:));
nFrames = size(firstIqBeam,2);

firstBeamIq = firstIqBeam(:,1:nFrames-1);
secondBeamIq = firstIqBeam(:,2:nFrames);
speedSound = 1540;
demodulationFrequency = s.iq.fDemodIQ_Hz;
frameRate = s.Framerate_fps;
totalPhaseShift = angle(mean(conj(secondBeamIq).*firstBeamIq));

conversionValue = (speedSound*frameRate)/(4*pi*demodulationFrequency);

calculatedVelocity = conversionValue*totalPhaseShift;

startAngle = s.iq.StartAngleIQ_rad;
angleIncrement = s.iq.AngleIncrementIQ_rad;

velocityInPhantom = zeros(1,nFrames-1);

for i = 1:nFrames-1
    
    velocityInPhantom(i) = abs(calculatedVelocity(i))*cos(startAngle+(i-1)*...
        angleIncrement);
end % for i

nSeconds = nFrames/frameRate;   %seconds = (nFrames/(nFrames/seconds))

time = 0:nSeconds/(nFrames-1):nSeconds;
time = time(1:630);

figure,plot(time,calculatedVelocity,time,velocityInPhantom),...
    title('Piston velocity estimated and in phantom');

legend('Estimated',...
    'Estimated with phase shift','Location','southeast');