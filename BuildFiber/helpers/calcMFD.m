function [D4SigmaX, D4SigmaY, clean_image, x, y] = calcMFD(image, threshold)
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
xCenter = sum(sum(x.*image))/sum(sum(image));
yCenter = sum(sum(y.*image))/sum(sum(image));
D4SigmaX = 4*sqrt( sum(sum((x-xCenter).^2.*image))/sum(sum(image)) );
D4SigmaY = 4*sqrt( sum(sum((y-yCenter).^2.*image))/sum(sum(image)) );

% Iterate the center and the size of the mode field by increasing the 
% threshold value until the size doesn't change significantly or the
% threshold is larger than the maximumm of the image.
flag = 1;
while flag > 0
    % Clean the image
    clean_image = image;
    clean_image(clean_image<threshold) = 0;
    
    % Find the new mode field info
    xCenter = sum(sum(x.*clean_image))/sum(sum(clean_image));
    yCenter = sum(sum(y.*clean_image))/sum(sum(clean_image));
    D4SigmaXTmp = 4*sqrt( sum(sum((x-xCenter).^2.*clean_image))/sum(sum(clean_image)) );
    D4SigmaYTmp = 4*sqrt( sum(sum((y-yCenter).^2.*clean_image))/sum(sum(clean_image)) );
    
    % Increase the threshold for the next iteration
    threshold = threshold*1.1;
    
    % Check if the conditions are fulfilled
    if D4SigmaXTmp/D4SigmaX > 0.98 || D4SigmaYTmp/D4SigmaY > 0.98 || threshold > max(image(:))
        flag = 0;
    end
    
    % Update the mode field size
    D4SigmaX = D4SigmaXTmp;
    D4SigmaY = D4SigmaYTmp;
end

x = x - xCenter;
y = y - yCenter;

end