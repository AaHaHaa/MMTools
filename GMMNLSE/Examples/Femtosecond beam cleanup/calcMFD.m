function [D4SigmaX, D4SigmaY, x, y] = calcMFD(image)
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

[M, N] = size(image);
x = -N/2:N/2-1;
y = (-M/2:M/2-1)';

% Mode field info
denominator = sum(sum(image.^2));
xCenter = sum(sum(x.*image.^2))/denominator;
yCenter = sum(sum(y.*image.^2))/denominator;
D4SigmaX = 4*sqrt( sum(sum((x-xCenter).^2.*image.^2))/denominator );
D4SigmaY = 4*sqrt( sum(sum((y-yCenter).^2.*image.^2))/denominator );

x = x - xCenter;
y = y - yCenter;

end