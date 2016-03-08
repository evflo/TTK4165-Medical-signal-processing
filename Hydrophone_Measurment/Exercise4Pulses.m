function Exercise4Pulses()
%
% Purpose: Plot pulse in both time and frequency (with noise)
%
% Made by: 
%   Even Florenes NTNU 2016
%
% Last changes:
%   2016-01-26: First attempt
%   2016-01-28: Added opportunity for plotting in time with kPa as axis
%
    pulses = {'pulse_7cm_0dB_10x_aver.mat','pulse_7cm_-8dB_10x_aver.mat','pulse_7cm_-14dB_10x_aver.mat','pulse_7cm_-14dB_no_aver.mat'};
    time = figure('Name','Pulses plotted in time');
    frequency = figure('Name','Pulses plotted in frequency');
    for i = 1:4
        load(pulses{i})
        %RF = RF*1e6*(1/25.1);
        figure(time),subplot(2,2,i),plot(RF),xlabel('Time range[\mus]'),ylabel('Amplitude[kPa]');
        % Task 2
        RFAveraged = AverageInput(RF,10);
        [P,f] = puls_fft(RFAveraged(1:200),Ts,2000,[0.5e6 6e6]);
        figure(frequency),subplot(2,2,i),plot(f,20*log10(P)),xlabel('Frequency range[Hz]'),ylabel('Amplitude[dB]');
        hold on
        % Task 2 noise
        [PNoise,f] = puls_fft(RF(end-200:end),Ts,2000,[0.5e6 6e6]);
        figure(frequency),subplot(2,2,i),plot(f,20*log10(PNoise));
        keyboard
    end
    figures = [time frequency];
    for i = 1:2
        figure(figures(i));
        subplot(2,2,1),title(pulses{1});
        subplot(2,2,2),title(pulses{2});
        subplot(2,2,3),title(pulses{3});
        subplot(2,2,4),title(pulses{4});
    end

    hold off