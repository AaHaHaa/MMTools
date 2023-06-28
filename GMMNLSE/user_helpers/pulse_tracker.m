function [ output_field,multipulsing,time_delay ] = pulse_tracker( input_field,dt,time_delay,energy )
%PULSE_TRACKER Track the pulse to keep it at the center of the time window
%in case it's too close to the boundaries of the time window and causes the
%difficult convergence of simulations.
%
% "energy" is to start the pulse-centering action only if the pulse
% starts to form and its energy array starts to converge.

switch nargin
    case 1
        time_delay = 0;
    case {3,4}
    otherwise
        error('pulse_tracker:InputArgumentError',...
        'The number of input arguments should be 1, 3, or 4.');
end

if isstruct(input_field)
    field = input_field.fields(:,:,end);
    dt = input_field.dt;
else
    field = input_field;
end

multipulsing = false;

Nt = size(field,1);

%% Track the pulse center
switch nargin
    case 4
        least_rt = 5; % the least number of round trips necessary to take. It should be at least larger than 3.
        tol_close_to_pulse = 0.008;

        diff_num = 3; % the size of the considered difference elements
        rt_num = find(energy==0,1)-1;
        if rt_num <= diff_num+2
            last_energy = energy(1:rt_num);
        else
            last_energy = energy( ( rt_num-diff_num-1):rt_num );
        end

        % Determine the time delay and shift the pulse to its center.
        if rt_num > least_rt
            difference = relative_central_diff(last_energy);

            [field,time_delay] = main(field,Nt,dt,time_delay);

            if ~any(abs(difference) > tol_close_to_pulse)
                multipulsing = is_multipulsing(field);
            end
        end
    case 3
        [field,time_delay] = main(field,Nt,dt,time_delay);
        
        multipulsing = is_multipulsing(field);
    case 1
        field = main(field,Nt);
        
        multipulsing = is_multipulsing(field);
end

if isstruct(input_field)
    output_field = struct('dt',dt,...
                          'fields',field);
else
    output_field = field;
end

end

%% relative_central_diff
function diff_array = relative_central_diff(array)
% DIFF_ARRAY calculate the relative slope from central difference and
% divide it by the center value.

diff_array = (array(3:end)-array(1:end-2))/2./array(2:end-1);

end

function [field,time_delay] = main(field,Nt,dt,time_delay)

% Single pulse
tCenter = floor(sum(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(field).^2),2)/sum(sum(abs(field).^2),2));
if ~isnan(tCenter) && tCenter ~= 0 % all-zero fields
    % Because circshift is slow on GPU, I discard it.
    %field = circshift(field,-tCenter);
    if tCenter > 0
        field = [field(1+tCenter:end,:);field(1:tCenter,:)];
    elseif tCenter < 0
        field = [field(end+1+tCenter:end,:);field(1:end+tCenter,:)];
    end
    if nargin > 2
        time_delay = time_delay + tCenter*dt;
    else
        time_delay = 0;
    end
end

end

%% is_multipulsing
function multipulsing = is_multipulsing(field)
%IS_MULTIPULSING
%
%   Find the most prominent peak and its FWHM, and calculate its pulse energy within a certain range.
%   If this energy takes up only a small portion of the total energy, there's definitely another pulse.

multipulsing = false;

total_field = sum(abs(field).^2,2);

% Smooth the pulse to remove some jagged feature
smooth_range = 10;
total_field = conv(total_field,ones(smooth_range,1),'same');

min_peak_height = max(total_field)/1.001;
[pks,locs,width] = findpeaks(total_field,'MinPeakHeight',min_peak_height,'WidthReference','halfheight','MinPeakProminence',min_peak_height/1.1); % "minPeakProminence" is to get rid of small peaks of the fluctuating pulse near the main peak
[~,highest_peak_idx] = max(pks);
fwhm = width(highest_peak_idx);
peak_loc = locs(highest_peak_idx);
span_left = max(1, floor(peak_loc-fwhm/2*3));
span_right = min(length(total_field), ceil(peak_loc+fwhm/2*3));
pulse_energy_ratio = trapz(total_field(span_left:span_right))/trapz(total_field);
if pulse_energy_ratio < 0.7 % there should be a second pulse that takes up a non-negligible energy
    [~,locs] = findpeaks(total_field,'MinPeakHeight',min_peak_height/2,'MinPeakDistance',2*fwhm);
    if length(locs) > 1
        multipulsing = true;
        warning('off','backtrace');
        if length(locs) == 2
            warning('Double-pulsing!');
        else % More than 2 pulses
            warning('Multi-pulsing!');
        end
        warning('on','backtrace');
        %{
        if length(locs) == 2
            message = 'Double-pulsing!';
        else % More than 2 pulses
            message = 'Multi-pulsing!';
        end
        if exist('cprintf','file')
            cprintf('*[1 0.5 0.31]','%s\n',message);
        else
            fprintf('%s\n',message);
        end
        %}
    end
end
        
end