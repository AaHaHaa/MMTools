function E_tr_noise = shot_noise(sim,...
                                 Nt,dt,...
                                 field_size,...
                                    dx,dy)
%SHOT_NOISE It computes the shot noise included in the governing equation

h = 6.62607015e-34; % J*s

time_window = Nt*dt; % ps
f = ifftshift((-Nt/2:Nt/2-1)'/time_window,1); % THz
real_f = (f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_f(real_f<0) = 0; % no noise at negative frequencies

% Find the spatial distribution of the spontaneous Raman scattering
spatial_distribution = 1/sqrt(dx*dy);

noise_amplitude = sqrt(h*real_f.*spatial_distribution.^2/(time_window*1e-12));

E_tr_noise = fft(noise_amplitude.*randn(field_size).*exp(1i*2*pi*rand(field_size)));

end