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

%%
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
%%
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
%%
pistonLength = 10;
fixedPoint = 4;
pistonAngularFrequency = (2*pi*(time-t0))/rotationPeriod;

%pistonVelocity = pistonVelocityAmplitude*sin(pistonVelocityAngularFrequency);
pistonPosition = pistonLength+excenterDistance*cos(pistonAngularFrequency);
pointPosition = fixedPoint+excenterDistance*cos(pistonAngularFrequency);

plot(time,pointPosition,'y'),hold off;

%%
pistonVelocityAmplitude = -(2*pi*excenterDistance)/rotationPeriod;
pistonVelocity = pistonVelocityAmplitude*sin(pistonAngularFrequency);

pointVelocity = -pistonVelocity;

figure,plot(time,pointVelocity),title('Point velocity'),xlabel('Time [sec]'),...
    ylabel('Position [cm]');


