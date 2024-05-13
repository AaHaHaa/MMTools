function [D4SigmaX, D4SigmaY] = calcMFD(image,spatial_window)
%CALCMFD This calculates the mode field of an image.
%It finds the size of the cleaned central field by calculating its second
%moment of x and y axes, followed by a multiplicaiton of 4. Besides, the
%final cleaned image and the center of the mode field are also returned.
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
%
% image - (Nx,Ny,num_images); in current implementation, Nx=Ny

image = abs(image);

[Nx, Ny, ~] = size(image);
dx = spatial_window/Nx;
dy = dx;
x = (-Nx/2:Nx/2-1)'*dx;
y = (-Ny/2:Ny/2-1)*dy;

% Mode field info
denominator = sum(image.^2,[1,2]);
xCenter = sum(x.*image.^2,[1,2])./denominator;
yCenter = sum(y.*image.^2,[1,2])./denominator;
D4SigmaX = 4*sqrt( sum((x-xCenter).^2.*image.^2,[1,2])./denominator );
D4SigmaY = 4*sqrt( sum((y-yCenter).^2.*image.^2,[1,2])./denominator );

end