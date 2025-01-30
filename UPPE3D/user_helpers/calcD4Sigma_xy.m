function [D4SigmaX, D4SigmaY] = calcD4Sigma_xy(field,spatial_window,varargin)
%CALCD4SIGMA_XY This calculates the D4Sgima of spatial fields.
%
% It finds the size of the cleaned central field by calculating its second
% moment of x and y axes, followed by a multiplicaiton of 4
%
% field: a (Nx,Nx,num_fields) (multi)dimensional array; the spatial "field", not intensity (be careful!)
% spatial_window: a scalar; the size of the spatial window
%
% optional input arguments:
%   remove_noise_model:
%       0 (default): don't remove the noise
%       1: remove the noise that is too weak based on thresholding
%       2: remove the noise by calculating I.^C and then multiply the found MFD by sqrt(C) times
%          This calculation assumes a Gaussian spatial profile.
%
% =========================================================================
% Reference: https://www.rp-photonics.com/beam_radius.html
%
% The recommended definition of the beam size is that of ISO Standard 11146,
% based on the second moment of the intensity distribution I(x,y).
%
% The common method is called D4Sigma because for the beam diameter, one
% obtains 4 times the standard deviation of the intensity distribution.
%
% =========================================================================

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

[Nx, Ny, num_fields] = size(field);
dx = spatial_window/Nx;
dy = dx;
x = (-floor(Nx/2):floor((Nx-1)/2))'*dx;
y = (-floor(Ny/2):floor((Ny-1)/2))*dy;

% Intensity
I = abs(field).^2;
% Removen noise for the mode-field calculation
switch remove_noise_model
    case 1
        for j = 1:num_fields
            idx = I(:,:,j) < max(max(I(:,:,j),[],1),[],2)/30; % noise region
            Itmp = I(:,:,j);
            max_noise = max(max(Itmp(idx),[],1),[],2);
            min_noise = min(min(Itmp(idx),[],1),[],2);
            Itmp(idx) = Itmp(idx).*((Itmp(idx)-min_noise)./(max_noise-min_noise)).^2; % remove noise by incrementally reducing their strength (weaker one becomes weaker)
            I(:,:,j) = Itmp;
        end
    case 2
        multiplication_ratio = 3;
        I = I./max(max(I,[],1),[],2); % make I from 0-1
        I = I.^multiplication_ratio;
end

denominator = sum(I,[1,2]);
xCenter = sum(x.*I,[1,2])./denominator;
yCenter = sum(y.*I,[1,2])./denominator;
D4SigmaX = 4*sqrt( sum((x-xCenter).^2.*I,[1,2])./denominator );
D4SigmaY = 4*sqrt( sum((y-yCenter).^2.*I,[1,2])./denominator );

if remove_noise_model == 2
    D4SigmaX = D4SigmaX*sqrt(multiplication_ratio);
    D4SigmaY = D4SigmaY*sqrt(multiplication_ratio);
end

end