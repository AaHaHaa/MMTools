function [G,saturation_parameter] = Gaussian_gain(fiber,sim,omegas)
%GAUSSIAN_GAIN It computes the gain parameters for the Gaussian-gain model.
%
% Input arguments:
%   fiber.saturation_energy
%   fiber.saturation_intensity
%   fiber.fr
%   fiber.gain_coeff
%   fiber.gain_fwhm
%   sim.gain_model
%   sim.f0
%   sim.gpu_yes
%   omegas = 2*pi*ifftshift(linspace(-floor(Nt/2), floor((Nt-1)/2), Nt))'/(Nt*dt); % in 1/ps, in the order that the fft gives

c = 2.99792458e-4; % speed of ligth m/ps

G = [];
saturation_parameter = [];
if sim.gain_model == 1
    if isscalar(sim.midx) % single-mode
        saturation_parameter = fiber.saturation_energy;
    else
        saturation_parameter = fiber.saturation_intensity;
    end
    if sim.gpu_yes
        saturation_parameter = gpuArray(saturation_parameter);
    end
    
    w_fwhm = 2*pi*sim.f0^2/c*fiber.gain_fwhm;
    w_0 = w_fwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665
    G = fiber.gain_coeff/2.*exp(-omegas.^2/w_0^2);
    
    if sim.gpu_yes
        G = gpuArray(G);
    end
end

end

