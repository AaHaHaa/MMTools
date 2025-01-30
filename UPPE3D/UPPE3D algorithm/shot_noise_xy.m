function E_tr_noise = shot_noise_xy(f0,...
                                    Nt,dt,...
                                    field,...
                                       dx,dy)
%SHOT_NOISE_XY It computes the shot noise included in the governing
% equation of 3D-UPPE.
%
% Because the shot noise approach inherently relies on
% mode-resolved formulation, the initial mode field is assumed as the input
% single mode. This is a strong assumption. For example, if the fiber
% supports multimode, shot noise should have higher energy = num_modes*h*nu.
% In addition, each mode will have a different spatial profile. A
% full-field approach requires recovery back to the mode picture for
% correct estimate of the shot-noise intensity (W/m^2) at each spatial
% position. Even for a single-mode propagation, the beam spatial profile
% can vary, leading to a changing shot-noise intensity. For example, a
% Gaussian beam, when propagating in free space, diverges.

field_size = size(field);
I = squeeze(sum(abs(field).^2,1)); % size: (Nx,Ny)
Aeff = (sum(I*dx*dy,[1,2]))^2/sum(I.^2*dx*dy,[1,2]); % effective mode-field area

h = 6.62607015e-34; % J*s

time_window = Nt*dt; % ps
f = ifftshift((-floor(Nt/2):floor((Nt-1)/2))'/time_window,1); % THz
real_f = (f+f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_f(real_f<0) = 0; % no noise at negative frequencies

noise_amplitude = sqrt(h*real_f/Aeff/(time_window*1e-12));

E_tr_noise = fft(noise_amplitude.*randn(field_size).*exp(1i*2*pi*rand(field_size)),[],1);

end