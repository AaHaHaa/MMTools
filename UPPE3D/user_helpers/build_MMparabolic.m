function output = build_MMparabolic(tfwhm, time_window, total_energy, num_modes, Nt, varargin)
%BUILD_MMPARABOLIC Build a multimode temporally-parabolic pulse using the following parameters:
%
% tfwhm - full width at half maximum of pulse, in ps
% time_window - width of entire time window, in ps
% total_energy - total energy of the pulse in all modes, in nJ
% num_modes - number of modes
% Nt - number of time grid points
%
% Optional inputs (varargin):
%   frequency_shift - a cell with two elements:
%                     the amount of shifted frequency (THz) (default is 0)
%                     and
%                     Fourier Transform type: 'fft' or 'ifft' (default is 'ifft')
%	coeffs - the normalized complex amplitude coefficients of the different modes (default is equal across all modes)
%   t_center - temporal position of the pulse in the time window (default is 0)
%
% Note:
%   Nonlinear Fiber Optics by Agrawal defines 'ifft' for Fourier Transform.
%   This convention is used as a default here.

%% Default optional input arguments
% Accept only 4 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 4
    error('build_MMparabolic:TooManyInputs', ...
          'It takes only at most 3 optional inputs');
end

% Set defaults for optional inputs
frequency_shift = {'ifft',0};
coeffs = ones(1,num_modes); % the energy is spreaded into all the modes equally
optargs = {frequency_shift,coeffs,0};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[frequency_shift, coeffs, t_center] = optargs{:};

% "coeffs" needs to be a row vector
if size(coeffs,2) == 1
    coeffs = coeffs.';
end
coeffs = coeffs./sqrt(sum(abs(coeffs).^2)); % normalization

%% Parabolic fields
% pulse power = B - A(t-t_center)^2
C = tfwhm/sqrt(2); % C = sqrt(B/A)
B = (total_energy*1e3)/(4/3*C); % "*1e3" is make it into pJ
A = B/C^2;

dt = time_window/Nt;  % ps
t = (-floor(Nt/2):floor((Nt-1)/2))'*dt; % ps

% Construct a single parabolic electric field envelope, in W^0.5
time_profile = sqrt(B - A*(t-t_center).^2);
time_profile(t<t_center-C | t>t_center+C) = 0;

% Apply the frequency shift
switch frequency_shift{1}
    case 'ifft'
        time_profile = time_profile.*exp(-1i*(2*pi*frequency_shift{2}).*t);
    case 'fft'
        time_profile = time_profile.*exp(1i*(2*pi*frequency_shift{2}).*t);
    otherwise
        error('build_MMparabolic:frequency_shiftError',...
              'The type of the Fourier Transform can only be ''ifft'' or ''fft''.');
end

% Apply this time profile to each mode using the coefficients
field = coeffs.*time_profile;

% Output as a struct
output = struct('fields',field,'dt',dt);

end