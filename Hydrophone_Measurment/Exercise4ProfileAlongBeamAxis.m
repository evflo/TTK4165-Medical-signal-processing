function Exercise4ProfileAlongBeamAxis()
%
% Purpose: Plot profile along the beam axis versus the depth
%
% Made by:
%   Even Florenes NTNU 2016
%
% Last changes:
%   2016-01-28: First attempt
%
speedOfSound = 1540;

load('z_axis_12dB.mat');

depth = Position(4,:)*speedOfSound;

intensity = sum(RF.^2,1);
intensitydB = 10*log10(intensity);
intensitydBNormalized = intensitydB-max(intensitydB);


figure,plot(depth,intensitydBNormalized),title('Profile along beam axis'),ylabel('Intensity[dB]'),xlabel('Depth');