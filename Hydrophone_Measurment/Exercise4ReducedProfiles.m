function Exercise4ReducedProfiles()
%
% Purpose: Plot the reduced azimuth and elevation profiles in both X and
% Y-direcetion
%
% Made by:
%   Even Florenes NTNU 2016
%
% Last changes:
%   2016-01-28: First attempt
%
profiles = {'profile_AZ_7cm_reduced','profile_EL_7cm_reduced',};
profileNaming = {'profile AZ 7cm Reduced','profile EL 7cm Reduced'};
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