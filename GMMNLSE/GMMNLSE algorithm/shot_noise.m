function At_noise = shot_noise(Nt,dt,sim,num_modes)
%SHOT_NOISE It computes the shot noise included in the governing equation

h = 6.62607015e-34; % J*s

time_window = Nt*dt; % ps
f = ifftshift((-floor(Nt/2):floor((Nt-1)/2))'/time_window,1); % THz
real_f = (f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_f(real_f<0) = 0; % no noise at negative frequencies

noise_amplitude = sqrt(h.*real_f/(time_window*1e-12));

if isequal(sim.step_method,'RK4IP')
    At_noise = fft(noise_amplitude.*randn(Nt,num_modes).*exp(1i*2*pi*rand(Nt,num_modes)),[],1);
else % 'MPA'
    At_noise = repmat(fft(noise_amplitude.*randn(Nt,1,num_modes).*exp(1i*2*pi*rand(Nt,1,num_modes)),[],1),1,sim.MPA.M+1,1);
end

end