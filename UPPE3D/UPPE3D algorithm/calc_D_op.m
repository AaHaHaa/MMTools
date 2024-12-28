function [D_op,W_op,...
          kc,k0,...
          n,...
          sim] = calc_D_op(sim,...
                           n,...
                           Nt,dt,...
                           Nx,dx,...
                           Ny,dy,...
                           omegas,omegas_real,...
                           fields,...
                           F_op)
%D_OP It computes the dispersion operator used in 3D-GNLSE

kx = 2*pi*ifftshift(linspace(-floor(Nx/2), floor((Nx-1)/2), Nx),2)/(Nx*dx); % in 1/m, in the order that the fft gives
ky = 2*pi*ifftshift(linspace(-floor(Ny/2), floor((Ny-1)/2), Ny),2)/(Ny*dy); % in 1/m, in the order that the fft gives
ky = permute(ky,[1,3,2]);

% The dispersion term in the 3D-UPPE, in frequency space
n = ifftshift(n,1); % (Nt,Nx,Ny,1,Np); in the order that the fft gives
c = 299792458e-12; % m/ps
k0 = omegas_real/c; % 1/m
k = n.*k0;
% Refractive index is separated into the space-invariant and variant parts.
% The invariant part is taken at where the highest refractive index is.
% For example, for graded-index fibers, it's at the fiber center.
% The space-variant part becomes the waveguide term.
[~,max_n_idx] = max(feval(@(x)x(:),sum(real(n),[1,5]))); [max_n_idx1,max_n_idx2] = ind2sub([Nx,Ny],max_n_idx);
nc = n(:,max_n_idx1,max_n_idx2,:,:);% max(max(max(real(n),[],2),[],3),[],5);
kc = nc.*k0;

kW2 = k.^2 - kc.^2; % waveguide contribution
W_op = 1i*kW2/2./kc;

% Obtain the omega0 of the input pulse
fftshift_omegas = fftshift(omegas,1);
spectrum = sum(abs(fftshift(F_op.Fs(fields),1)).^2,[2,3,5]);
omega0 = sum(fftshift_omegas.*spectrum)/sum(spectrum); % 2*pi*THz; the pulse center frequency (under shifted omega)
%{
idx0 = find(fftshift_omegas > omega0,1);
if idx0 >= Nt/2
    idx0 = idx0 - Nt/2 +1;
else
    idx0 = idx0 + Nt/2 -1;
end

% Get the waveguide term only near the pulse frequency
W_op = W_op(idx0,:,:,:,:);
%}

if ~isfield(sim,'betas')
    if Nt == 1 || ~any(fields(:)) % CW case
        sim.betas = [0;0];
    else
        sim.betas = zeros(2,1,'gpuArray');

        % Obtain the betas of the input pulse
        omega_range = 1/dt; % 2*pi*THz
        omegas_idx_near_pulse = fftshift_omegas>omega0-omega_range/5 & fftshift_omegas<omega0+omega_range/5;% pick only the data near the pulse center frequency to find its beta0 and beta1
        clear spectrum omega0 omega_range;

        fftshift_kc = fftshift(kc,1);
        fit_order = max(2,min(7,sum(omegas_idx_near_pulse)-1)); % 2~7
        [betas_Taylor_coeff,~,mu] = polyfit(fftshift_omegas(omegas_idx_near_pulse),real(fftshift_kc(omegas_idx_near_pulse)),fit_order);
        sim.betas = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
        sim.betas = [sim.betas(1)-sim.betas(2)*mu(1)/mu(2);...
                     sim.betas(2)/mu(2)];
        clearvars fftshift_omegas fftshift_kc fit_order betas_Taylor_coeff mu
    end
end

Kz = sqrt(complex(kc.^2 - (kx.^2+ky.^2))); % the dispersion term; in k-space
D_op = 1i*(Kz-(sim.betas(1)+sim.betas(2)*omegas));

% Remove the high spatial frequencies and fill it with zeros
% It's lossy (negative real part) here, so making it zero is equivalent but
% it makes the adaptive step-size control simpler and faster.
D_op(repmat(kc,1,length(kx),length(ky),1,1) < sqrt(kx.^2+ky.^2)) = 0;

end