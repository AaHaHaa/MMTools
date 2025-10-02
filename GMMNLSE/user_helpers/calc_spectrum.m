function [spectrum,spectrum_wavelength] = calc_spectrum(fields,time_window,f)
%CALC_SPECTRUM This computes the spectrum in the frequency and wavelength
%domains
%
% Inputs:
%   fields: (Nt,num_modes,Nz); electric fields (sqrt(W))
%   time_window: time window (ps)
%   f: (Nt,1); frequency points
%
% Outputs:
%   spectrum: (Nt,num_modes,Nz); the spectrum in the frequency domain (nJ/THz)
%   spectrum_wavelength: (Nt,num_modes,Nz); the spectrum in the wavelength domain (nJ/nm)

c = 299792458*1e-12; % m/ps
lambda = c./f*1e9; % nm

factor_correct_unit = time_window^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                         % "/1e3" is to make pJ into nJ
spectrum = abs(fftshift(ifft(fields),1)).^2*factor_correct_unit; % in frequency domain

factor = (c*1e9)./lambda.^2; % change the spectrum from frequency domain into wavelength domain
spectrum_wavelength = spectrum.*factor; % nJ/nm

end

