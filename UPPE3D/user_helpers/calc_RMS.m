function [RMS,T1] = calc_RMS(x,y)
%CALC_RMS It calculates RMS width
%
%   x: a column or row vector
%   y: a multidimensional array composed of "column vectors".
%      y should be intensities of pulses or spectra, instead of complex-number fields

sx = size(x);
if length(sx)==2 && sx(1)==1
    x = x';
end

area = trapz(x,y);

T1 = trapz(x,x.*y)./area;
T2 = trapz(x,x.^2.*y)./area;

RMS = sqrt(T2-T1.^2);

end