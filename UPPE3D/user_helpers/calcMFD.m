function [D4SigmaX, D4SigmaY] = calcMFD(field,spatial_window)
%CALCMFD This calculates the D4Sgima of spatial fields.
%It finds the size of the cleaned central field by calculating its second
%moment of x and y axes, followed by a multiplicaiton of 4
%
% field: a (Nx,Nx,num_fields) (multi)dimensional array; the spatial "field", not intensity (be careful!)
% spatial_window: a scalar; the size of the spatial window
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

[Nx, Ny, num_fields] = size(field);
dx = spatial_window/Nx;
dy = dx;
x = (-Nx/2:Nx/2-1)'*dx;
y = (-Ny/2:Ny/2-1)*dy;

% Intensity
I = abs(field).^2;
for j = 1:num_fields
    idx = I(:,:,j) < max(max(I(:,:,j)))/10;
    Itmp = I(:,:,j);
    max_noise = max(max(Itmp(idx)));
    min_noise = min(min(Itmp(idx)));
    Itmp(idx) = Itmp(idx).*((Itmp(idx)-min_noise)./(max_noise-min_noise)).^2; % remove noise by incrementally reducing their strength (weaker one becomes weaker)
    I(:,:,j) = Itmp;
end

denominator = sum(I,[1,2]);
xCenter = sum(x.*I,[1,2])./denominator;
yCenter = sum(y.*I,[1,2])./denominator;
D4SigmaX = 4*sqrt( sum((x-xCenter).^2.*I,[1,2])./denominator );
D4SigmaY = 4*sqrt( sum((y-yCenter).^2.*I,[1,2])./denominator );

end