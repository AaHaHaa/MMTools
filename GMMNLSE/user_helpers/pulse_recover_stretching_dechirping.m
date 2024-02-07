function field = pulse_recover_stretching_dechirping(field0,recover_info,varargin)
%PULSE_RECOVER_STRETCHING_DECHIRPING Summary of this function goes here
%   Detailed explanation goes here

optargs = {false,[]};
optargs(1:length(varargin)) = varargin;
[verbose, time] = optargs{:};

field = field0.*recover_info{1};
field_w = fftshift(ifft(field),1); % spectrum
field_w = field_w.*recover_info{2};
field = fft(ifftshift(field_w,1))./recover_info{1};

% Shift the pulse temporally
field0(abs(field0)<max(abs(field0))/3) = 0; % remove noise for the original field
% 1. Shift the pulse to where the previous peak is first:
% This initial step is important in case the stretched pulse goes
% beyond the time window and fails the intensity-weighted
% computation.
[~,pulse_center0] = max(abs(field0).^2);
[~,pulse_center] = max(abs(field).^2);
index_shift = round(pulse_center - pulse_center0);
field = double(circshift(field,-index_shift,1));
% 2. Then shift the pulse according to intensity-weighted pulse
% center.
% This assumes that the input doesn't have a field that wrap around
% the time window.
field1 = field; field1(abs(field1)<max(abs(field1))/3) = 0; % remove noise for the recovered field
[~,pulse_center0] = calc_RMS((1:length(field0))',abs(field0).^2);
[~,pulse_center] = calc_RMS((1:length(field1))',abs(field1).^2);
index_shift = round(pulse_center - pulse_center0);
field = circshift(field,-index_shift,1);

if verbose && length(varargin) == 2
    show_result(time,field);
end

end

%%
function show_result(time,field)
%SHOW_RESULT

intensity = abs(field).^2;

% Plot only the central part of the time window
intensity_plot = intensity;
threshold_factor = 100;
intensity_plot(intensity<max(intensity)/threshold_factor) = 0;
left = find(intensity_plot~=0,1);
right = find(intensity_plot~=0,1,'last');
center = floor((left+right)/2);
span_factor = 2;
span = floor((right-left)/2)*span_factor;
left = max(1,floor(center-span));
right = min(length(time),ceil(center+span));

figure('Name','recovered pulse');
h = plot(time(left:right)*1e3,intensity(left:right));
xlim([min(time(left:right)) max(time(left:right))]*1e3);
xlabel('Time (fs)');
ylabel('Power (W)');
title('Recovered pulse');

set(h,'linewidth',2);
set(gca,'fontsize',14);

end