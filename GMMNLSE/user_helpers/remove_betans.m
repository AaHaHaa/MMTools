function [dechirped_field,varargout] = remove_betans( input_field,n_remove,n_all_fit,num_interp_points,t,verbose )
%REMOVE_BETANS It gives the dechirped result without specified beta_n.
%
% Input:
%   input_field: (Nt,...); the electric field in time domain (sqrt(W))
%   n_remove: up to the order of betas to remove
%             0 means to remove the constant (DC) phase
%             1 means to remove the first-order phase
%             2 means to remove the second-order phase
%             ......
%   n_all_fit: the order of betas in finding the Taylor-series coefficients
%             My experience is that the best reasonable value is 4.
%   num_interp_points: (optional input argument)
%                      a scalar;
%                      the number of points interpolated between two original data points
%                      (default: zero)
%
%   **Because the transform-limited pulse is really narrow, 
%     the original time spacing can be too large to correctly characterize the pulse.
%     This will be a problem with calculating the transform-limited duration.
%     Therefore, increasing the number of time grid points, i.e., reduing the time spacing, is highly recommended.
%
%   t: (optional input argument)
%      (Nt,1); time (ps)
%   verbose: (optional input argument)
%            plot or not; true or false
%
% Output:
%   dechirped_field
%   t_interp: the interpolated time grid points based on "num_interp_points"
%   dechirped_FWHM: dechirped pulse duration (fs)
%   pulse_FWHM: current pulse duration (fs)
% =========================================================================
% "num_interp_points" and "t" are optional input arguments, but "t" is
% required to calculate "t_interp" and "dechirped_FWHM" as output.
% =========================================================================
% Usage:
%   dechirped_field = calc_dechirped(input_field);
%   dechirped_field = calc_dechirped(input_field,num_interp_points);
%   [dechirped_field,t_interp,dechirped_FWHM] = calc_dechirped(input_field,num_interp_points,t);

calc_dechirped_duration = false;
switch nargin
    case 3
        num_interp_points = 0;
        verbose = false;
    case {5,6}
        calc_dechirped_duration = true;
end

sE = size(input_field);
Nt = sE(1);
interp_idx = linspace(1,Nt,(num_interp_points+1)*(Nt-1)+1)';

% The interpolation can't be done with only
% "interp1(input_field,interp_idx)" because an interpolation of complex
% numbers isn't correct. It should be done with two interpolations of their
% absolute values and phases.
if num_interp_points ~= 0
    abs_input_field = interp1(abs(input_field),interp_idx);
    phase_input_field = interp1(unwrap(angle(input_field)),interp_idx);
    input_field = abs_input_field.*exp(1i*phase_input_field);
end

input_field_f = fftshift(ifft(ifftshift(input_field,1)),1);
tmp_f = (1:size(input_field_f,1))';
[~,~,fitted_param] = characterize_spectral_phase( (1:size(input_field_f,1))',input_field_f,n_all_fit );
p = fitted_param{1}; p = [p(1:n_all_fit-n_remove),zeros(1,n_remove+1)];
s = fitted_param{2};
mu = fitted_param{3};
phase_no_betans = polyval(p,2*pi*tmp_f,s,mu);
dechirped_field = fftshift(fft(ifftshift(abs(input_field_f).*exp(1i*phase_no_betans),1)),1);

if calc_dechirped_duration
    % Interpolated time
    t_interp = interp1(t,interp_idx);
    
    pulse_FWHM = zeros([1,sE(2:end)]);
    dechirped_FWHM = zeros([1,sE(2:end)]);
    for i = 1:prod(sE(2:end))
        % Current duration
        threshold = max(abs(input_field(:,i)).^2)/1.01;
        [~,~,tmp_pulse_width,~] = findpeaks(abs(input_field(:,i)).^2,t_interp*1e3,'MinPeakHeight',threshold,'WidthReference','halfheight');
        pulse_FWHM(i) = tmp_pulse_width(1);
        
        % Dechirped duration
        threshold = max(abs(dechirped_field(:,i)).^2)/1.0001;
        [~,~,dechirped_FWHM(i),~] = findpeaks(abs(dechirped_field(:,i)).^2,t_interp*1e3,'MinPeakHeight',threshold,'WidthReference','halfheight');
    end
    
    varargout = {t_interp,dechirped_FWHM,pulse_FWHM};
    
    if verbose
        intensity = abs(dechirped_field).^2;
        intensity_plot = intensity;
        threshold_factor = 100;
        intensity_plot(intensity<max(intensity)/threshold_factor) = 0;
        left = find(intensity_plot~=0,1);
        right = find(intensity_plot~=0,1,'last');
        center = floor((left+right)/2);
        span_factor = 2;
        span = floor((right-left)/2)*span_factor;
        left = floor(center-span);
        right = ceil(center+span);
        if left < 1, left = 1; end
        if right > Nt*(num_interp_points+1), right = Nt*(num_interp_points+1); end

        [TL,TL_t] = calc_transform_limited(input_field,0,t_interp);
        figure;
        h = plot(t_interp(left:right)*1e3,abs(dechirped_field(left:right)).^2); set(h,'linewidth',2); hold on;
        h2 = plot(TL_t(left:right)*1e3,abs(TL(left:right)).^2); set(h2,'linewidth',2); hold off;
        xlim([min(TL_t(left:right)) max(TL_t(left:right))]*1e3);
        set(gca,'fontsize',20);
        xlabel('Time (fs)'); ylabel('Power (W)');
        legend('Dechirped','TL');
    end
else
    varargout = {};
end

end