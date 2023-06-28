function [ duration,bandwidth ] = calc_duration_bandwidth( t,wavelength_f,field,wavelength_f_type )
%CALC_DURATION_BANDWIDTH It calculates the RMS duration and bandwidth.
%
%   t: (N,1); time
%   wavelength_f: (N,1); wavelength (nm) or frequency (THz)
%   field: (N,......), a multidimensional array composed of columns of
%          fields to be calculated
%   wavelength_f_type: 'lambda' or 'f'

intensity = abs(field).^2;
intensity(intensity<max(intensity)/50) = 0;
spectrum = abs(fftshift(ifft(field),1)).^2;
spectrum(spectrum<max(spectrum)/50) = 0;

if isequal(wavelength_f_type,'lambda')
    c = 299792.458; % nm/ps
    factor = c./wavelength.^2; % change the spectrum from frequency domain into wavelength domain
    spectrum = spectrum.*factor;
end

duration = calc_RMS(t,intensity);
bandwidth = calc_RMS(wavelength_f,spectrum);

% fwhm
RMS2fwhm = 2*sqrt(log(2))*sqrt(2);
duration = duration*RMS2fwhm;
bandwidth = bandwidth*RMS2fwhm;

end