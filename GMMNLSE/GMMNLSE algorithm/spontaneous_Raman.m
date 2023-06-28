function sponRS = spontaneous_Raman(Nt,dt,sim,fiber,num_modes)
%SPONTANEOUS_RAMAN It computes the spontaneous Raman term.

h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

time_window = Nt*dt; % ps
f = ifftshift((-Nt/2:Nt/2-1)'/time_window,1); % THz
real_omegas = 2*pi*(f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_omegas(real_omegas<0) = 0; % no noise at negative frequencies

if sim.gpu_yes
    SR = zeros(1,num_modes,'gpuArray');
else
    SR = zeros(1,num_modes);
end
for midx = 1:num_modes
    if sim.scalar
        SR(:,midx) = fiber.SR(midx,midx,midx,midx);
    else
        if mod(midx,2) == 0
            SR(:,midx) = SR(:,midx/2);
        else
            SR(:,midx) = fiber.SR(ceil(midx/2),ceil(midx/2),ceil(midx/2),ceil(midx/2));
        end
    end
end

% Bose-Einstein distribution
temperature = 273.15 + 25; % 25 degree Celsius
nth = 1./(exp(h*abs(f*1e12)./k/temperature)-1); % phonon number based on Bose-Einstein distribution
nth(isinf(nth)) = 0; % if f=0
Heaviside = double(f<0); % Stokes wave

sponRS = {sqrt(hbar*real_omegas.*SR/(time_window*1e-12)),...
          nth+Heaviside};
end