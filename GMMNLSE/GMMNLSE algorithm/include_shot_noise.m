function fields = include_shot_noise(sim,omegas,Nt,dt,fields)
%INCLUDE_SHOT_NOISE It adds shot noise to the fields
%
% Input arguments:
%   sim.num_photon_noise_per_band
%   sim.f0
%   sim.gpu_yes
%   omegas
%   Nt
%   dt
%   fields

if sim.num_photon_noise_per_band ~= 0
    hbar = 6.62607015e-34/2/pi*1e24; % pJ*ps
    photon_noise_intensity = hbar*(omegas+2*pi*sim.f0)/(Nt*dt)*sim.num_photon_noise_per_band;
    % I use analytical-signal representation for solving GMMNLSE, so the field covers only the positive frequencies.
    photon_noise_intensity(photon_noise_intensity<0) = 0; % no noise at negative frequencies
    photon_noise_amplitude = sqrt(photon_noise_intensity);
    if sim.gpu_yes
        rand_phase = exp(1i*2*pi*rand(Nt,size(fields,2),'gpuArray'));
    else
        rand_phase = exp(1i*2*pi*rand(Nt,size(fields,2)));
    end
    fields = fft(ifft(fields) + photon_noise_amplitude.*rand_phase);

    clear hbar photon_noise_intensity photon_noise_amplitude rand_phase;
end

if sim.gpu_yes
    fields = gpuArray(fields);
end

end

