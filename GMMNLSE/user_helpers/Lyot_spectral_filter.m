function [output,fig] = Lyot_spectral_filter(input, f0, target_wavelength, thickness, varargin)
%LYOT_SPECTRAL_FILTER Apply a spectral Lyot filter to a field
%
% Input:
%   input.fields - a (N, num_modes, m) matrix with each mode's field, in the time domain
%   input.dt - the time grid spacing (ps)
%   f0 - the center frequency of simulations (THz)
%   target_wavelength - the targeted center wavelength of the Lyot filter (nm)
%   thickness - birefringent plate thickness (mm)
%
%   Optional inputs (varargin):
%       birefringent_angle
%       polarizer_angle
%       verbose - true(1) or false(0); whether to plot the input and the output spectra (default: false)
%
% Output:
%   output.fields
%   output.rejected_fields
%   output.dt
%   fig - the figure handle of the figure if "verbose" is true

c = 299792.458; % in nm/ps

if ~isstruct(input)
    error('Unlike most auxiliary functions, Lyot_spectral_filter requires that the input be a struct with at least fields and dt');
end

optargs = {true,0.0088};
optargs(1:length(varargin)) = varargin;
[verbose, dn] = optargs{:};

input_field = input.fields;
N = size(input_field, 1);
dt = input.dt;

f = f0 + ifftshift(linspace(-N/2, N/2-1, N))'/(N*dt); % in THz, in the order that the fft gives

c = 299792458; % m/s
phase_retardation = 2*pi*dn*(f*1e12)/c*(thickness*1e-3);

birefringent_angle = -pi/4; % 45 degree

target_phase_retardation = mod(2*pi*dn/(target_wavelength*1e-9)*(thickness*1e-3),2*pi);
polarizer_angle = -mod(target_phase_retardation, 2*pi)/2;

% Calculate the filter profile in frequency space
birefringent_phase_matrix = @(theta,eta) exp(-1i*eta/2).*[cos(theta)^2+exp(1i*eta)*sin(theta)^2,(1-exp(1i*eta))*cos(theta)*sin(theta)];
linear_polarizer_matrix = @(theta,Ex,Ey) [Ex*cos(theta)^2+Ey*cos(theta)*sin(theta),Ex*cos(theta)*sin(theta)+Ey*sin(theta)^2];
transmission = @(Ex,Ey) abs(Ex).^2+abs(Ey).^2;
quarter_wave_matrix = @(theta,Ex,Ey) exp(-1i*pi/4)*[(cos(theta)^2+1i*sin(theta)^2)*Ex+((1-1i)*sin(theta)*cos(theta))*Ey,((1-1i)*sin(theta)*cos(theta))*Ex+(sin(theta)^2+1i*cos(theta)^2)*Ey];

E = birefringent_phase_matrix(birefringent_angle,phase_retardation);
E = quarter_wave_matrix(birefringent_angle+pi/4,E(:,1),E(:,2));
E = linear_polarizer_matrix(polarizer_angle,E(:,1),E(:,2));
mult_factor = transmission(E(:,1),E(:,2));

% Apply the filter in frequency space
output = struct('dt',input.dt,...
                'fields',fft(mult_factor.*ifft(input_field)),...
                'rejected_fields',fft((1-mult_factor).*ifft(input_field)));

if verbose
    fig = figure('Name','Filter');
    f = fftshift(f,1); % THz
    wavelength = 299792.458./f; % nm
    factor_correct_unit = (N*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                        % "/1e3" is to make pJ into nJ
    factor = c./wavelength.^2; % change the spectrum from frequency domain into wavelength domain
    spectrum = abs(fftshift(ifft(input_field(:,:,end)),1)).^2*factor_correct_unit.*factor;
    
    [wavelength_rms,wavelength_c] = calc_RMS(wavelength,spectrum);
    
    yyaxis left
    h1 = plot(wavelength,spectrum);
    hold on;
    h2 = plot(wavelength,spectrum.*fftshift(mult_factor,1).^2,'--');
    hold off;
    xlabel('Wavelength (nm)'); ylabel('Spectrum (nJ/nm)');
    xlim(wavelength_c+[-1,1]*wavelength_rms*6);
    title('Spectral filter');
    set(h1,'linewidth',2); set(h2,'linewidth',2);
    set(gca,'fontsize',14);
    yyaxis right
    h3 = plot(wavelength,fftshift(mult_factor,1).^2);
    ylabel('Filter');
    set(h3,'linewidth',2);
    ylim([0,1]);
    drawnow;
end

end

