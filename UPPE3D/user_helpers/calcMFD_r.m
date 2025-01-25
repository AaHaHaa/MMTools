function MFD = calcMFD_r(field,r,varargin)
%CALCMFD_R This calculates the MFD of radially-symmetric spatial fields.
%
% It finds the size of the cleaned field by calculating its effective
% mode-field area. By assuming a circular profile, the mode-field radius is
% computed with
%   MFR = sqrt(Aeff/pi)
%   --> MFD = MFR*2
%
% field: a (Nr,num_fields) two-dimensional array; the spatial "field", not intensity (be careful!)
% r: (1,Nr); radial sampling positions (m)
%
% optional input arguments:
%   remove_noise_model:
%       0 (default): don't remove the noise
%       1: remove the noise that is too weak based on thresholding
%       2: remove the noise by calculating I.^C and then multiply the found MFD by sqrt(C) times
%          This calculation assumes a Gaussian spatial profile.

%% Default optional input arguments
% Accept only 1 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('calcMFD:TooManyInputs', ...
          'It takes only at most 1 optional input');
end

% Set defaults for optional inputs
optargs = {0};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
remove_noise_model = optargs{:};

%%
if size(r,1) == 1
    r = r'; % r is a column vector in this code
end
num_fields = size(field,2);

% Intensity
I = abs(field).^2;
% Removen noise for the mode-field calculation
switch remove_noise_model
    case 1
        for j = 1:num_fields
            idx = I(:,j) < max(I(:,j),[],1)/30; % noise region
            Itmp = I(:,j);
            max_noise = max(Itmp(idx),[],1);
            min_noise = min(Itmp(idx),[],1);
            Itmp(idx) = Itmp(idx).*((Itmp(idx)-min_noise)./(max_noise-min_noise)).^2; % remove noise by incrementally reducing their strength (weaker one becomes weaker)
            I(:,j) = Itmp;
        end
    case 2
        multiplication_ratio = 3;
        I = I./max(I,[],1); % make I from 0-1
        I = I.^multiplication_ratio;
end

Aeff = (2*pi*trapz(r,I.*r)).^2./(2*pi*trapz(r,I.^2.*r)); % effective mode-field area (m^2)
MFD = sqrt(Aeff/pi)*2;

if remove_noise_model == 2
    MFD = MFD*sqrt(multiplication_ratio);
end

end