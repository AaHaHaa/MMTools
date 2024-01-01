function fields = include_shot_noise(sim,omegas,Nt,dt,fields)
%INCLUDE_SHOT_NOISE It adds shot noise to the fields
%
% Input arguments:
%   sim.num_photon_noise_per_bin: the number of photon noise per spectral discretization bin
%   sim.f0: center frequency of the numerical frequency window (THz);
%           This is used to compute the "real_omegas" for the photon noise.
%   sim.gpu_yes: whether to use GPU or not
%   omegas: relative omegas of the frequency window (THz)
%   Nt: the number of numerical sampling points
%   dt: the temporal spacing (ps)
%   fields: the electric field (sqrt(W))
%
% -------------------------------------------------------------------------
% Unit explanation:
%   intensity = abs(field).^2;
%   energy = trapz(t,intensity) = trapz(intensity)*dt;       % pJ
%   
%   spectrum_unknown_unit = abs(fftshift(ifft(field),1)).^2;
%
%   Parseval's theorem: sum(intensity) = sum(spectrum_unknown_unit)*N;
%                       * Note that spectrum_unknown_unit is from "ifft".
%   therefore sum(intensity)*dt = sum(spectrum_unknown_unit)*N*dt
%                               = sum(spectrum_unknown_unit)*(N*dt)^2/(N*dt)
%                               = sum(spectrum_unknown_unit)*(N*dt)^2*df
%                               = sum(spectrum_f)*df
%
%   spectrum_f = spectrum_unknown_unit*(N*dt)^2;
%   energy = trapz(f,spectrum_f) = trapz(spectrum_f)*df      % pJ
%                                = trapz(spectrum_f)/(N*dt);
%                                = trapz(spectrum_unknown_unit)*(N*dt)
%
% -------------------------------------------------------------------------
%   c = 299792.458;     % nm/ps
%   wavelength = c./f;  % nm
%   spectrum_wavelength = spectrum_f.*(c./wavelength.^2);
%   energy = -trapz(wavelength,spectrum_wavelength);         % pJ
% -------------------------------------------------------------------------
%   The noise is added as "one photon per frequency band," so each
%   frequency band adds one photon energy, hbar*omega, to the total energy.
%   This corresponds to its "spectrum_unknown_unit" counterpart as
%   "hbar*omega/(N*dt)," whose field amplitude is
%   "sqrt(hbar*omega/(N*dt))."
%
%   If the frequency window is small, assume that all omegas=omegas0,
%   adding photon noise adds N*hbar*omegas0 to the total energy,
%   proportional to the number of number of numerical sampling points.
% -------------------------------------------------------------------------
%   Note that the implementation of having photon noise being
%   hbar*omegas/(N*dt) relies on having "ifft" as Fourier Transform. The
%   relation might be different when Fourier Transform becomes "fft"
%   because they have different constants to divide in the integral
%   relations.
% -------------------------------------------------------------------------

if sim.num_photon_noise_per_bin ~= 0
    hbar = 6.62607015e-34/2/pi*1e24; % pJ*ps
    photon_noise_intensity = hbar*(omegas+2*pi*sim.f0)/(Nt*dt)*sim.num_photon_noise_per_bin;
    % I use analytical-signal representation for solving GMMNLSE, so the field covers only the positive frequencies.
    photon_noise_intensity(photon_noise_intensity<0) = 0; % no noise at negative frequencies
    photon_noise_amplitude = sqrt(photon_noise_intensity);
    if sim.gpu_yes
        rand_value = exp(1i*2*pi*rand(Nt,size(fields,2),'gpuArray'));
    else
        rand_value = exp(1i*2*pi*rand(Nt,size(fields,2)));
    end
    fields = fft(ifft(fields) + photon_noise_amplitude.*rand_value);

    clear hbar photon_noise_intensity photon_noise_amplitude rand_value;
end

if sim.gpu_yes
    fields = gpuArray(fields);
end

end

