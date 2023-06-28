function [center_frequency,frequency_range,...
          time_window,dt] = find_tw_f0(frequency_range,Nf)
%FIND_TW_F0 This code find the center frequency and the time window
%required for the desired range of frequencies and the number of frequency
%points.

center_frequency = (frequency_range(2)+frequency_range(1))/2; % THz
dt = 1/abs(frequency_range(2)-frequency_range(1)); % ps
time_window = Nf*dt;

% 
df = 1/time_window;
df = center_frequency/ceil(center_frequency/df);
time_window = 1/df;
dt = time_window/Nf;
frequency_range = center_frequency + [floor((Nf-1)/2),-floor(Nf/2)]*df;

end