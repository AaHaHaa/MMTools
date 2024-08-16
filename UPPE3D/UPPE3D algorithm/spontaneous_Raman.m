function sponRS = spontaneous_Raman(sim,fr,...
                                    Nt,dt,...
                                       dx,dy)
%SPONTANEOUS_RAMAN It computes the spontaneous Raman term.
%
%   It creates the spontaneous-Raman counterpart of |E|^2, E: the pulse field
%
%   I consider one-photon noise per frequency bin.
%   spectral noise field = counterpart of E
%                        = sponRS_prefactor{1}.*randn(size(sponRS_prefactor{1})).*exp(1i*2*pi*rand(size(sponRS_prefactor{1})))
%   noise temporal intensity = abs( fft(spectral noise field) ).^2
%   Transform into the spectral domain for convolution = ifft( noise temporal intensity ).*sponRS_prefactor{2}
%   sponRS_prefactor{2} modifies the spectral noise according to the Bose-Einstein distribution and Stokes generation.
%   sponRS_Gamma = fft(haw.*sponRS) finishes the spontaneous-Raman convolution
%       such that
%   sponRS_Gamma*E is the spontaneous-Raman field in the 3D-GNLSE.
%
% -------------------------------------------------------------------------
%   They're summarized and implemented as the following the stepping functions:
%
%      sponRS = ifft(abs(fft(sponRS_prefactor{1}.*randn(size(sponRS_prefactor{1})).*exp(1i*2*pi*rand(size(sponRS_prefactor{1}))))).^2).*sponRS_prefactor{2};
%      sponRS_Gamma = fft(haw.*sponRS);

h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

time_window = Nt*dt; % ps
f = ifftshift((-Nt/2:Nt/2-1)'/time_window,1); % THz
real_omegas = 2*pi*(f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_omegas(real_omegas<0) = 0; % no noise at negative frequencies

% Find the spatial distribution of the spontaneous Raman scattering
spatial_distribution = 1/sqrt(dx*dy);

% Bose-Einstein distribution
temperature = 273.15 + 25; % 25 degree Celsius
nth = 1./(exp(h*abs(f*1e12)./k/temperature)-1); % phonon number based on Bose-Einstein distribution
nth(isinf(nth)) = 0; % if f=0
Heaviside = double(f<0); % Stokes wave

sponRS = {sqrt(fr*hbar*real_omegas.*spatial_distribution.^2/(time_window*1e-12)),...
          nth+Heaviside};
end