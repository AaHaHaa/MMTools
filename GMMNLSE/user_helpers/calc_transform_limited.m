function [transform_limited_field,varargout] = calc_transform_limited( input_field,num_insert_points,t )
%DECHIRP_TRANSFORM_LIMITED It gives the dechirped result which is transform-limited.
%
% Input:
%   input_field: (Nt,...); the electric field in time domain (sqrt(W))
%   num_insert_points: (optional input argument)
%                      a scalar;
%                      the number of points inserted between two original data points
%                      (default: zero)
%
%   **Because the transform-limited pulse is really narrow, 
%     the original time spacing can be too large to correctly characterize the pulse.
%     This will be a problem with calculating the transform-limited duration.
%     Therefore, increasing the number of time grid points, i.e., reduing the time spacing, is highly recommended.
%
%   t: (optional input argument)
%      (Nt,1); time (ps)
%
% Output:
%   transform_limited_field
%   t_insert: the interpolated time grid points based on "num_insert_points"
%   transform_limited_FWHM: transform-limited pulse duration (fs)
%   pulse_FWHM: current pulse duration (fs)
% =========================================================================
% "num_insert_points" and "t" are optional input arguments, but "t" is
% required to calculate "t_insert" and "transform_limited_FWHM" as output.
% =========================================================================
% Usage:
%   transform_limited_field = calc_transform_limited(input_field);
%   transform_limited_field = calc_transform_limited(input_field,num_insert_points);
%   [transform_limited_field,t_insert,transform_limited_FWHM] = calc_transform_limited(input_field,num_insert_points,t);

calc_TL_duration = false;
switch nargin
    case 1
        num_insert_points = 0;
    case 3
        calc_TL_duration = true;
end

sE = size(input_field);
Nt = sE(1);
insert_idx = [linspace(1,Nt,(num_insert_points+1)*(Nt-1)+1)'; Nt+(1:num_insert_points)'/(num_insert_points+1)];

% Interpolation
if num_insert_points ~= 0
    input_freq = ifft(input_field);
    input_freq = cat(1,input_freq(1:ceil(Nt/2),:),zeros(Nt*num_insert_points,1),input_freq(ceil(Nt/2)+1:end,:));
    input_field = reshape(fft(input_freq),[Nt*(num_insert_points+1),sE(2:end)]);
end

input_freq_TL = abs(fftshift(ifft(ifftshift(input_field,1)),1));
transform_limited_field = fftshift(fft(ifftshift(input_freq_TL,1)),1);

if calc_TL_duration
    % Inserted time
    t_insert = interp1(t,insert_idx,'linear','extrap');
    
    pulse_FWHM = zeros([1,sE(2:end)]);
    transform_limited_FWHM = zeros([1,sE(2:end)]);
    for i = 1:prod(sE(2:end))
        % Current duration
        threshold = max(abs(input_field(:,i)).^2)/1.01;
        [~,~,tmp_pulse_width,~] = findpeaks(abs(input_field(:,i)).^2,t_insert*1e3,'MinPeakHeight',threshold,'WidthReference','halfheight');
        pulse_FWHM(i) = tmp_pulse_width(1);
        
        % Transform-limited duration
        threshold = max(abs(transform_limited_field(:,i)).^2)/1.0001;
        [~,~,transform_limited_FWHM(i),~] = findpeaks(abs(transform_limited_field(:,i)).^2,t_insert*1e3,'MinPeakHeight',threshold,'WidthReference','halfheight');
    end
    
    varargout = {t_insert,transform_limited_FWHM,pulse_FWHM};
else
    varargout = {};
end

end