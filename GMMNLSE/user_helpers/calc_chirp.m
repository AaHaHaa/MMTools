function func = calc_chirp
%CALC_CHIRP This is a caller for two functions, Gaussian() and General().
%Read below for how to use them.
%
%   The chirp is the C in "exp(i*C/2* (w-w0)^2)" under frequency domain.
%   Adding a quadratic phase to a spectrum is equivalent to adding a linear
%   chirp under the time domain.
%
%   Use:
%       func = calc_chirp;
%
%       chirp_sign = 1; % for positive chirp (normal dispersive delay)
%
%       duration = ?; % ps
%       time_window = ?; % ps
%       omega = ifftshift(2*pi*(-N/2:N/2-1)'/time_window,1); % 2*pi*THz
%       spectrum_amplitude = ifft(E); % E is the electric field
%       
%       For a Gaussian pulse,
%       [chirp,chirped_pulse] = func.Gaussian( duration,omega,spectrum_amplitude,chirp_sign );
%       
%       For a random pulse,
%       [chirp,chirped_pulse] = func.General( duration,omega,spectrum_amplitude,chirp_sign );

func.Gaussian = @Gaussian;
func.General  = @General;

end

function  [chirp,chirped_pulse] = Gaussian( duration,omega,spectrum_amplitude,chirp_sign )
%CALC_CHIRP It calculates the chirp and the chirped pulse of a Gaussian 
%pulse based on its duration and bandwidth.
%
%   Please check "Agrawal, Ch.3.2, Gaussian pulse" for details.
%   However, C in Agrawal doesn't directly correspond to the quadratic 
%   coefficient I define above.
%   Check its Eq(3.2.19), Nonlinear Fiber Optics, v5, and assume that the 
%   pulse starts without chirp (C=0 in that equation). beta2*z is exactly 
%   the quadratic phase coefficient that adds to the spectrum during 
%   propagation. Therefore, beta2*z is, in fact, the "C" we want here.
%   Check Eq(3.2.5) for details as well. We then obtain:
%
%       T1/T0 = sqrt( 1 + (C/T0^2) )
%
%   Input:
%       duration:  FWHM under time domain; ps
%       omega: the angular frequency column vector (before fftshift); 2*pi*THz
%       spectrum_amplitude: the spectrum amplitude column vector which is the "Fourier Transform of the electric field".
%                           It's a complex column vector. Maybe it's not correct to call it "amplitude".
%                           Don't take "fftshift" and "abs()^2" for this input argument.
%       chirp_sign: add a positive chirp or a negative chirp; 1 or -1
%
%   Output:
%       chirp: the chirp coefficient or the quadratic phase coefficient;
%              (2*pi*THz)^(-2), but we usually don't show 2*pi factor in the unit, so it becomes ps^2
%       chirped_pulse: the electric field column vector for the chirped pulse

duration = duration(1); % in case it's an array

% Spectrum information
fftshift_omega = fftshift(omega,1);
abs2_spectrum = abs(fftshift(spectrum_amplitude,1)).^2;
omega0 = trapz(fftshift_omega,fftshift_omega.*abs2_spectrum)./trapz(fftshift_omega,abs2_spectrum);

% Gaussian bandwidth
bandwidth = calc_RMS(fftshift_omega,abs2_spectrum)*sqrt(2)*(2*sqrt(log(2))); % fwhm

t1 = duration/(2*sqrt(log(2))); % ps
w0 = bandwidth/(2*sqrt(log(2))); % 2*pi*THz

t0 = 1/w0; % calculate t0 from the time-bandwidth product; ps
           % For Gaussian, its time-bandwidth product equals one.

% t1 = t0*sqrt( 1+(C/t0^2)^2 ), from "Agrawal, Ch.3.2, Gaussian pulse"
if t1 > t0
    chirp = sign(chirp_sign)*t0*sqrt(t1^2-t0^2); % ps^2
else
    chirp = 0;
end

chirped_pulse = fft(spectrum_amplitude.*exp(1i*chirp/2*(omega-omega0).^2));

end

function [chirp,chirped_pulse] = General( duration,omega,spectrum_amplitude,chirp_sign )
%GENERAL It calculates the chirp of a spectrum for the desired pulse
%duration.
%
%   This code uses an optimization function to find the chirp that creates 
%   the desired pulse duration for a specific spectrum.
%
%   Input:
%       duration: FWHM under time domain; ps
%       omega: the angular frequency column vector (before fftshift); 2*pi*THz
%       spectrum_amplitude: the spectrum amplitude column vector which is the "Fourier Transform of the electric field".
%                           It's a complex column vector. Maybe it's not correct to call it "amplitude".
%                           Don't take "fftshift" and "abs()^2" for this input argument.
%       chirp_sign: having positive chirp or a negative chirp eventually; 1, -1, or 0
%                   This is relative to the transform-limited pulse duration.
%                   1 means the pulse has a positive chirp to obtain desired pulse duration;
%                   -1 is related to the negative chirp.
%                   0 represents transform-limited pulse.
%
%   Output:
%       chirp: the chirp coefficient or the quadratic phase coefficient;
%              (2*pi*THz)^(-2), but we usually don't show 2*pi factor in the unit, so it becomes ps^2
%       chirped_pulse: the electric field column vector for the chirped pulse

Nt = length(omega);
dt = 2*pi/(max(omega)-min(omega)); % ps
time = (ceil(-Nt/2):ceil(Nt/2-1))*dt; % ps

% Find omega0
fftshift_omega = fftshift(omega,1);
abs2_spectrum = abs(fftshift(spectrum_amplitude,1)).^2;
omega0 = trapz(fftshift_omega,fftshift_omega.*abs2_spectrum)./trapz(fftshift_omega,abs2_spectrum);

% Kill the background spectral noise, particularly the high-frequency noise
spectrum_amplitude_thresholded = spectrum_amplitude;
spectrum_amplitude_thresholded(abs(spectrum_amplitude)<max(abs(spectrum_amplitude))/10) = 0;

% Guess the initial chirp by assuming a Gaussian shape
input_duration = calc_RMS(time,abs(fft(spectrum_amplitude)).^2)*sqrt(2)*(2*sqrt(log(2))); % fwhm
%input_bandwidth = calc_RMS(fftshift_omega,abs2_spectrum)*sqrt(2)*(2*sqrt(log(2))); % fwhm
guess_chirp0 = Gaussian( input_duration,omega,spectrum_amplitude,1 );
% Guess the chirp required to obtain the desired duration
guess_desired_chirp = Gaussian( duration,omega,spectrum_amplitude,chirp_sign );

% Dechirp the pulse first
options = optimset('TolX',1e-20);
%options = optimset('PlotFcns',@optimplotfval,'TolX',1e-20); % plot the process of optimization
find_dechirp = @(C) find_fwhm( C,time,0,omega,omega0,spectrum_amplitude_thresholded );
% Because I don't know the input pulse is initially positively or
% negatively chirped, I run the optimization twice with different initial
% guess and find the optimum.
[dechirp_pos,diff_pos] = fminsearch(find_dechirp, guess_chirp0,options);
[dechirp_neg,diff_neg] = fminsearch(find_dechirp,-guess_chirp0,options);
if diff_pos < diff_neg
    dechirp = dechirp_pos;
else
    dechirp = dechirp_neg;
end

% Add the chirp based on chirp_sign
find_optimal_chirp = @(C) find_fwhm( abs(C)*sign(chirp_sign)+dechirp,time,duration,omega,omega0,spectrum_amplitude_thresholded );
stretching_chirp = fminsearch(find_optimal_chirp,guess_desired_chirp,options);
chirp = abs(stretching_chirp)*sign(chirp_sign) + dechirp;
chirped_pulse = fft(spectrum_amplitude.*exp(1i*chirp/2*(omega-omega0).^2));

end

function difference = find_fwhm( C,time,duration,omega,omega0,spectrum_amplitude )
%FIND_FWHM

pulse = sum(abs(fft(spectrum_amplitude.*exp(1i*C/2*(omega-omega0).^2))).^2,2);

% Find the FWHM
pulse(pulse<max(pulse)/50) = 0; % It's important to kill the noise first; otherwise, noise will make calc_RMS return a large value
pulse_FWHM = calc_RMS(time,pulse)*sqrt(2)*(2*sqrt(log(2))); % transform from t0 to tfwhm

if isempty(pulse_FWHM)
    pulse_FWHM = 0;
else
    pulse_FWHM = max(pulse_FWHM);
end
difference = abs(duration - pulse_FWHM);

end

function RMS = calc_RMS(x,y)
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