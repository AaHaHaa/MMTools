function sponRS = spontaneous_Raman(Nt,dt,sim)
%SPONTANEOUS_RAMAN It computes the spontaneous Raman term.
% Spontaneous Raman scattering is equivalent to scattering with a field
% with one photon per frequency band.

h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

time_window = Nt*dt; % ps
f = ifftshift((-Nt/2:Nt/2-1)'/time_window,1); % THz
real_omegas = 2*pi*(f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_omegas(real_omegas<0) = 0; % no noise at negative frequencies

% Bose-Einstein distribution
temperature = 273.15 + 25; % 25 degree Celsius
nth = 1./(exp(h*abs(f*1e12)./k/temperature)-1); % phonon number based on Bose-Einstein distribution
nth(isinf(nth)) = 0; % if f=0
Heaviside = double(f<0); % Stokes wave

sponRS = {sqrt(hbar.*real_omegas/(time_window*1e-12)),...
          nth+Heaviside};
end