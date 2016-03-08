function Exercise4FullProfiles()
%
% Purpose: Plot azimuth and elevation profiles in both X and Y-direction
%
% Made by:
%   Even Florenes NTNU 2016
%
% Last changes: 
%   2016-01-26: First attempt
%   2016-01-28: Bug fixing
%
profiles = {'profile_AZ_7cm_Full','profile_EL_7cm_Full',};
profileNaming = {'profile AZ 7cm Full','profile EL 7cm Full'};
profile = figure;
subPlotNum = 1;
for i = 1:2
    load(profiles{i});
    % Find intensity in dB and normalize
    Intensity = sum(RF.^2,1);
    Intensity = 10*log10(Intensity);
    Intensity = Intensity-max(Intensity);
    if i == 1
        X = Position(1,:);
        Y = Position(3,:);
    else
        X = Position(2,:);
        Y = Position(3,:);
    end % if i==1
    figure(profile),subplot(2,2,subPlotNum),plot(X,Intensity),title([profileNaming{i} ' in azimuth direction']),ylabel('Intensity[dB]'),xlabel('X[mm]');
    subPlotNum = subPlotNum+1;
    figure(profile),subplot(2,2,subPlotNum),plot(Y,Intensity),title([profileNaming{i} ' in elevation direction']),ylabel('Intensity[dB]'),,xlabel('Y[mm]');
    subPlotNum = subPlotNum+1;
end