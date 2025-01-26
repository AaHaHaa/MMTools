function output = build_MMsoliton(tfwhm, beta2, fiiii, lambda0, time_window, num_modes, Nt,varargin)
%BUILD_MMGAUSSIAN Build soliton pulses using the following parameters:
%
% tfwhm - full width at half maximum of pulse, in ps; (1,num_modes)
% beta2 - group velocity dispersion, ps^2/m; (1,num_modes)
% fiiii - overlap integral; fiiii=1/Aeff_i, in m^(-2); (1,num_modes)
% lambda0 - the central wavelength, in m; (1,num_modes)
% time_window - the width of entire time window, in ps
% num_modes - the number of modes
% Nt - the number of time grid points
%
% Optional inputs:
%   frequency_shift - a cell with two elements:
%                     Fourier Transform type: 'fft' or 'ifft' (default is 'ifft')
%                     and
%                     the amount of shifted frequency (THz) (default is 0)
%                     
%   N - soliton order (default is 1, fundamental soliton)
%   center - temporal position of the pulse in the time window (default is 0)
%   n2 - nonlinear refractive index, in m^2 W^(-1) (default is 2.3e-20)
%
% Note:
%   Nonlinear Fiber Optics by Agrawal defines 'ifft' for Fourier Transform.
%   This convention is used as a default here.

%% Default optional input arguments
% Accept only 4 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 4
    error('build_MMsoltion:TooManyInputs', ...
        'It takes only at most 4 optional inputs');
end

% Set defaults for optional inputs
frequency_shift = {'ifft',0};
optargs = {frequency_shift,1,0,2.3e-20};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[frequency_shift, N, center, n2] = optargs{:};

%% Check validity
% Check num_modes
check_num_modes(tfwhm,num_modes);
check_num_modes(beta2,num_modes);
check_num_modes(fiiii,num_modes);
check_num_modes(lambda0,num_modes);
check_num_modes(N,num_modes);
check_num_modes(frequency_shift{2},num_modes);
check_num_modes(center,num_modes);
check_num_modes(n2,num_modes);

% "tfwhm, beta2, fiiii, lambda0, N, center, n2" need to be scalars or row vectors
tfwhm = change_to_row(tfwhm);
beta2 = change_to_row(beta2);
fiiii = change_to_row(fiiii);

%% General parameters for solitons
lambda0 = change_to_row(lambda0);
N = change_to_row(N);
center = change_to_row(center);
n2 = change_to_row(n2);

t0 = tfwhm/(asech(1/sqrt(2))*2);    % ps; 2*sqrt(log(2))=1.665 is for Gaussian
                                    %     2*asech(1/sqrt(2))=1.7627 is for fundamental solitons
dt = time_window/Nt;  % ps
t = (-floor(Nt/2):floor((Nt-1)/2))'*dt; % ps

%% Soliton fields
% Construct a single sech2 electric field envelope, in W^0.5
LD = t0.^2./abs(beta2);
c = 299792458;
lambda0 = c./(c./lambda0 + frequency_shift{2}*1e12); % calibrate lambda0 by including the frequency shift
gamma = n2*2*pi./lambda0.*fiiii; % m/W
P0 = abs(N.^2./gamma./LD);
fields = sqrt(P0).*sech((t-center)./t0);

if size(fields,2)==1
    fields = repmat(fields,1,num_modes);
end

% Apply the frequency shift
switch frequency_shift{1}
    case 'ifft'
        fields = fields.*exp(-1i*(2*pi*frequency_shift{2}).*t);
    case 'fft'
        fields = fields.*exp(1i*(2*pi*frequency_shift{2}).*t);
    otherwise
        error('build_MMsoliton:frequency_shiftError',...
              'The type of the Fourier Transform can only be ''ifft'' or ''fft''.');
end

% Output as a struct
output = struct('fields',fields,'dt',dt);

end

%% Sub-functions
function v = change_to_row(v)

if size(v,1) ~= 1 && size(v,2) == 1
    v = v.';
end

end

function check_num_modes(v,num_modes)

if ~ismember(length(v),[1 num_modes])
    error('build_MMsoliton:num_modesError',...
        'Please check the input arguments. They need to match with "num_modes".');
end

end