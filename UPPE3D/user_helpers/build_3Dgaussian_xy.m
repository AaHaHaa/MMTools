function output = build_3Dgaussian_xy(MFD, spatial_window, tfwhm, time_window, energy, Nt, Nx, varargin)
%BUILD_3DGAUSSIAN_XY Build a 3D temporally-superGaussian pulse with a Gaussian spatial profile using the following parameters:
%
% MFD - mode-field diameter (m)
% spatial_window - width of the spatial domain (m)
% Nx - number of spatial grid points
%
% tfwhm - full width at half maximum of pulse (ps)
% time_window - width of entire time window (ps)
% Nt - number of time grid points
%
% energy - pulse energy (nJ)
%
% Optional inputs (varargin):
%   frequency_shift - a cell with two elements:
%                     the amount of shifted frequency (THz) (default is 0)
%                     and
%                     Fourier Transform type: 'fft' or 'ifft' (default is 'ifft')
%   t_center - temporal position of the pulse in the time window (default is 0)
%   x_center - x spatial position (default is 0)
%   y_center - y spatial position (default is 0)
%   gaussexpo - supergaussian exponent (~exp(-t^(2*gaussexpo))) (default is 1)
%
% Note:
%   Nonlinear Fiber Optics by Agrawal defines 'ifft' for Fourier Transform.
%   This convention is used as a default here.

%% Default optional input arguments
% Accept only 5 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 5
    error('build_3Dgaussian:TooManyInputs', ...
        'It takes only at most 5 optional inputs');
end

% Set defaults for optional inputs
frequency_shift = {'ifft',0};
optargs = {frequency_shift,0,0,0,1};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[frequency_shift, t_center, x_center, y_center, gaussexpo] = optargs{:};

%% Gaussian temporal profile
temporal_profile = build_MMgaussian(tfwhm, time_window, energy, 1, Nt, frequency_shift, 1, t_center, gaussexpo);

%% Gausian spatial profile
spatial_profile = build_2Dgaussian_xy(MFD, spatial_window, Nx, x_center, y_center, gaussexpo);

%% Output as a struct
output = struct('field',temporal_profile.fields.*permute(spatial_profile.field,[3,1,2]),...
                'dt',temporal_profile.dt,'dx',spatial_profile.dx,...
                                         'dy',spatial_profile.dx);

end