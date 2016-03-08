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

figure,imagesc(time,1000*distance,logarithimicPower),colormap(gray),...
    title('M-Mode of Iq'),xlabel('Time[sec]'),ylabel('Distance[mm]');


