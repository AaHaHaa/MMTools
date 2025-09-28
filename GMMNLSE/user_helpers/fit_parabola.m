function [fitted_profile,deviation,fitted_param] = fit_parabola(t,At,varargin)
%FIT_PARABOLA It fits an intensity profile abs(At)^2 to a parabola
% Input arguments:
%   t: time (ps); (Nt,1)
%   At: pulse's field in time (sqrt(W)); (Nt,1)
% Optional argument:
%   verbose: whether to show the figure or not;
%            true or false
%
% Output arguments:
%   fitted_profile: fitted temporal profile of the parabola
%   deviation: difference between the fitted profile and the parabola
%   fitted_param: details of the fitting data with MATLAB's polyfit


%% Default optional input arguments
% Accept only 1 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('fit_parabola:TooManyInputs', ...
          'It takes only at most 1 optional input');
end

% Set defaults for optional inputs
verbose = true;
optargs = {verbose};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[verbose] = optargs{:};

input_t = t; % keep the original data for generating fitted profile later

%% Consider only the central (non-zero) part of the profile
% Remove the weak spectral signal for correct bandwidth computation with the RMS scheme
intensity = abs(At).^2;
threshold_factor = 20;
intensity(intensity<max(intensity)/threshold_factor) = 0;

% RMS computation
[duration,center_t] = calc_RMS(t,intensity); duration = duration*2;
[~,center_idx] = min(abs(t-center_t));
% Find the left and right edges of the spectrum
duration_idx = ceil(duration/((t(end)-t(1))/(length(t)-1)));
span_idx = ceil(duration_idx*1.5);
left_idx = max(1,center_idx - span_idx);
right_idx = min(length(t),center_idx + span_idx);

t = t(left_idx:right_idx);
intensity = intensity(left_idx:right_idx);

% Temporal fitting region should be small, only around the pulse's center
left_duration_idx = max(1,ceil(length(t)/2) - floor(duration_idx*0.7));
right_duration_idx = min(length(t),ceil(length(t)/2) + floor(duration_idx*0.7));
% Fitting is considered only for strong intensity
left_half_intensity = find(intensity > intensity/2,1);
right_half_intensity = find(intensity > intensity/2,1,'last');
left_duration_idx = max(left_duration_idx,left_half_intensity);
right_duration_idx = min(right_duration_idx,right_half_intensity);
% Final fitting region
fitting_region = left_duration_idx:right_duration_idx;

%% Fit to a parabola
fitted_order = 2;
[p,s,mu] = polyfit(t(fitting_region),intensity(fitting_region),fitted_order);
[fitted_profile,delta] = polyval(p,t,s,mu);
fitted_profile(fitted_profile<0) = 0;

fitted_param = {p,s,mu};

%deviation = sqrt(sum((intensity(fitting_region)-fitted_profile(fitting_region)).^2)/sum(intensity(fitting_region).^2));
deviation = sqrt(sum((intensity-fitted_profile).^2)/sum(intensity.^2));

%% Plot
if verbose
    intensity = abs(At).^2; intensity = intensity(left_idx:right_idx);
    % Normalization
    intensity = intensity/max(intensity);
    fitted_profile = fitted_profile/max(fitted_profile);
    
    % Plot spectrum and phases
    figure;
    hI = plot(t,intensity,'b');
    ylim([0,max(intensity)*1.5]);
    ylabel('Power (norm.)');
    set(gca,'YColor','b');
    hold on;
    hpI = plot(t,fitted_profile,'r');
    hold off;
    xlim([t(1),t(end)]);
    xlabel('Time (ps)');
    legend('Temporal profile','Fitted profile');
    set(hI,'linewidth',2);set(hpI,'linewidth',2);
    set(gca,'fontsize',14);
end

%% Output
[fitted_profile,delta] = polyval(p,input_t,s,mu);
fitted_profile(fitted_profile<0) = 0;

end

