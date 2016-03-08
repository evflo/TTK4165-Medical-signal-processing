function Exercise4TwoDimGrayscaleImage()
%
% Purpose: Plot 2D-gray scale image of both the full and reduced profile
%
% Made by:
%   Even Florenes NTNU 2016
%
% Last changes:
%   2016-01-28: First attempt
%
Profiles = {'profile_Az_7cm_full.mat','profile_Az_7cm_reduced.mat'};
ProfileNaming = {'profile Az 7cm full','profile Az 7cm reduced'};

speedofsound=1500;%m/s
TwoDimGrayScale = figure;
for i = 1:2
    load(Profiles{i});
    timescale=Position(4,1)+Ts(1)*[1:size(RF,1)]; %Ts= Sampling frequency 
    depthscale=timescale*speedofsound*1e3;%mm
    azimutpositionscale=Position(1,:)*1e3; %mm 
    figure(TwoDimGrayScale),subplot(1,2,i),imagesc(azimutpositionscale,depthscale,RF),title(ProfileNaming{i}); %Negative pressure> dark, positive > bright 
    colormap(gray);
end