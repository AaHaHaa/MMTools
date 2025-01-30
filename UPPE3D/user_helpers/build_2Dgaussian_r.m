function output = build_2Dgaussian_r(MFD, r, varargin)
%BUILD_2DGAUSSIAN_R Build a superGaussian radially-symmetric spatial profile using the following parameters:
%
% MFD - mode-field diameter (m)
% spatial_window - width of the spatial domain (um)
% r - (1,Nr); radial sampling positions (m)
%
% Optional inputs (varargin):
%   gaussexpo - supergaussian exponent (~exp(-t^(2*gaussexpo))) (default is 1)

%% Default optional input arguments
% Accept only 3 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('build_2Dgaussian:TooManyInputs', ...
          'It takes only at most 1 optional inputs');
end

% Set defaults for optional inputs
optargs = {1};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
gaussexpo = optargs{:};

%% Gausian spatial profile
MFR = MFD/2; % mode-field radius; m

gexpo = 2*gaussexpo;

spatial_profile = exp(-r.^gexpo/(MFR^gexpo));
spatial_profile = spatial_profile/sqrt(2*pi*trapz(r,abs(spatial_profile).^2.*r));

%% Output as a struct
output = struct('field',spatial_profile);

end