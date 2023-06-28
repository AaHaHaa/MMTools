function output = build_noisy_MMgaussian(tfwhm, tfwhm_noise, time_window, pulse_energy, noise_energy, num_modes, N, varargin)
%BUILD_NOISY_MMGAUSSIAN Generate an initial input pulse with random noise based on "build_MMgaussian".
%
% tfwhm - full width at half maximum of pulse, in ps
% tfwhm_noise - full width at half maximum of the noise, in ps
% time_window - width of entire time window, in ps
% pulse_energy - total energy of the small pulses in all modes, in nJ
% noise_energy - the energy of the background noise, in nJ
% num_modes - number of modes
% N - number of time grid points
%
% Optional inputs:
%	
%   noise_resolution - the resolution of the noise. The smaller it is, the noisier the noise is. Recommended 0.01 or lower.
%   frequency_shift - a cell with two elements:
%                     the amount of shifted frequency (THz) (default is 0)
%                     and
%                     Fourier Transform type: 'fft' or 'ifft' (default is 'ifft')
%   coeffs - the normalized complex amplitude coefficients of the different modes (default is equal across all modes)
%   center - temporal position of the pulse in the time window (default is 0)
%   gaussexpo - supergaussian exponent (~exp(-t^(2*gaussexpo))) (default is 1)

%% Default optional input arguments
if ~isempty(varargin)
    noise_resolution = varargin{1};
    if length(varargin)>1
        varargin = varargin(2:end);
    else
        varargin = {};
    end
else
    noise_resolution = 0.005;
end

if noise_energy == 0
        error('build_noisy_MMgaussian:NoiseEnergyError',...
            '"Noise energy" can''t be zero in this function. Please use "build_MMgaussian" for pure gaussian pulses.');
end

% Set defaults for optional inputs
frequency_shift = {'ifft',0};
coeffs = ones(1,num_modes);
coeffs = coeffs./sqrt(sum(abs(coeffs).^2));
optargs = {frequency_shift,coeffs,0,1};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:length(varargin)) = varargin;

% Place optional args in memorable variable names
[frequency_shift, coeffs, center, ~] = optargs{:};

% "coeffs" needs to be a row vector
if size(coeffs,2) == 1
    coeffs = coeffs.';
end

%% Gaussian main pulses
pulse_gaussian = build_MMgaussian(tfwhm, time_window, pulse_energy, num_modes, N, varargin{:});
dt = pulse_gaussian.dt;

%% Noise
if isinf(tfwhm_noise)
    noise_region = 1;
else
    noise_region = build_MMgaussian(tfwhm_noise,time_window,pulse_energy,1,N,frequency_shift,1,center,1); % pulse_energy isn't important here
    noise_region = noise_region.fields/max(noise_region.fields);
end

span_ratio = 2;

% Noise fields
noise_resolution_pulse = build_MMgaussian(time_window*noise_resolution,time_window,pulse_energy,1,N,frequency_shift,1,center,1); % pulse_energy isn't important here
noise_energy = noise_energy*coeffs.^2;
noise = exp(1i*rand(N,num_modes)*2*pi);
decay_region = build_MMgaussian(tfwhm*span_ratio,time_window,pulse_energy,1,N,frequency_shift,1,center,3);
decay_region = 1 - decay_region.fields/max(decay_region.fields);
noise = noise.*decay_region;
noise = fftshift(fft(ifft(noise).*ifft(noise_resolution_pulse.fields)),1).*noise_region;
current_noise_energy = trapz(abs(noise).^2)*dt/1e3; % nJ
noise_normalization = sqrt(noise_energy./current_noise_energy);
noise = noise_normalization.*noise;

%% Output
output.fields = pulse_gaussian.fields + noise;
output.dt = dt;

%{
% Compare energy
fprintf('noise energy: %6.4f(nJ)\n',noise_energy);
fprintf('pulse energy: %6.4f(nJ)\n',pulse_energy);
fprintf('total energy: %6.4f(nJ)\n',sum(trapz(abs(output.fields).^2)*dt/1e3));
%}
end