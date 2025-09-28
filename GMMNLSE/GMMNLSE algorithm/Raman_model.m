function [fiber,haw,hbw] = Raman_model(fiber,sim,Nt,dt)
%RAMAN_MODEL It calculates the Raman response of several solid materials.
%   Input:
%       fiber.material: a string; the type of the fiber;
%                         Currently it supports only silica, chalcogenide(As2S3), and ZBLAN(fluoride) fibers
%       sim: a structure containing
%           sim.include_Raman
%           sim.scalar
%           sim.ellipticity
%           sim.gpu_yes
%       Nt: the number of time points in the simulation
%       dt: the time interval (Nt*dt=time window)
%
%   Output:
%       fiber: Besides "material", "fiber.fr" is included based on the material
%       haw: the isotropic Raman response for the convolution under the frequency domain
%       hbw: the anisotropic Raman response for the convolution under the frequency domain;
%            Currently anisotropic one is only implemented in silica fibers
%
%   For the details about ZBLAN fibers, please refer to
%   "Supercontinuum generation in ZBLAN fibers-detailed comparison between measurement and simulation" by Agger (2012)

haw = [];
hbw = []; hb = zeros(Nt,1);

if ~sim.include_Raman
    return;
end

t_shifted = dt*(0:Nt-1)'; % time starting at 0
if sim.gpu_yes
    t_shifted = gather(t_shifted);
end

switch fiber.material
    case 'silica'
        % only the isotropic Raman
        % Ch. 2.3, p.42, Nonlinear Fiber Optics (5th), Agrawal
        %fiber.fr = 0.18; % 0.18 is standard for silica fibers
        %t1 = 12.2e-3; % Raman parameter t1 [ps]
        %t2 = 32e-3; % Raman parameter t2 [ps]

        %ha = ((t1^2+t2^2)/(t1*t2^2)).*exp(-t_shifted/t2).*sin(t_shifted/t1);
        
        % isotropic and anisotropic Raman
        % "Ch. 2.3, p.43" and "Ch. 8.5, p.340" in Nonlinear Fiber Optics (5th), Agrawal
        % It includes the "boson peak".
        % For more details, please read "Raman response function for
        % silica fibers", by Q. Lin and Govind P. Agrawal (2006)
        fiber.fr = 0.245;
        fa = 0.75;
        fb = 0.21;  % 1-fb = fa+fc
        fc = 0.04;

        t1 = 12.2e-3; % Raman parameter t1 [ps]
        t2 = 32e-3;   % Raman parameter t2 [ps]
        tb = 96e-3; % ps

        ha_model = ((t1^2+t2^2)/(t1*t2^2)).*exp(-t_shifted/t2).*sin(t_shifted/t1); % isotropic Raman
        hb_model = ((2*tb-t_shifted)/tb^2).*exp(-t_shifted/tb); % anisotropic Raman

        ha = fa*ha_model;
        hb = fc*ha_model + fb*hb_model;

        include_anisotropic_Raman = true;
    
    case 'chalcogenide'
        fiber.fr = 0.115;
        t1 = 15.2e-3; % Raman parameter t1 [ps]
        t2 = 230.5e-3; % Raman parameter t2 [ps]

        ha = ((t1^2+t2^2)/(t1*t2^2)).*exp(-t_shifted/t2).*sin(t_shifted/t1);

        include_anisotropic_Raman = false;
    
    case 'ZBLAN'
        fiber.fr = 0.062;
        
        a1 = 0.54e-13; % m/W
        a2 = 0.25e-13; % m/W
        v1 = 17.4; % THz
        v2 = 12.4; % THz
        w1 = 0.68; % THz
        w2 = 3.5;  % THz
        
        max_Stokes_freq_shift = 200; % 2*pi*THz; check Fig.4 in the reference paper mentioned above
        N_Omega = 4000; % just choose some number of points
        Omega = linspace(0,max_Stokes_freq_shift,N_Omega);
        ga = a1*exp(-(Omega/2/pi-v1).^2/(2*w1^2)) + a2*exp(-(Omega/2/pi-v2).^2/(2*w2^2)); % Raman gain profile
        % Because this step is easy to blow up the RAM, I use for-loop.
        % integral_part = trapz(Omega,ga.*sin(Omega.*t_shifted),2);
        integral_part = zeros(Nt,1);
        for i = 1:Nt
            integral_part(i) = trapz(Omega,ga.*sin(Omega.*t_shifted(i)),2);
        end
        
        c = 299792458; % m/s
        pump_wavelength = 1060e-9; % m; the pump wavelength for measuring this Raman response
        
        % fiber.fr is chosen to normalized ha
        % ha = integral_part/trapz(t_shifted,ha);
        ha = c/(fiber.fr*pi*fiber.n2*(2*pi*c/pump_wavelength))*integral_part;

        include_anisotropic_Raman = false;

    otherwise
        error('Raman_model:fiberMaterialError',...
              'Raman model is implemented only for silica, chalcogenide, and ZBLAN for now.');
end

% For scalar fields, anisotropic Raman is incorporated into ha to faciliate computations.
if include_anisotropic_Raman
    if sim.scalar
        if sim.ellipticity == 0 % linear polarization
            ha = ha + hb;
        else % circular polarization: its SRb=SRa/2, so the factor 1/2 is included here
            ha = ha + hb/2;
        end

        include_anisotropic_Raman = false; % its computation is integrated in isotropic part, so there's no need for a separate computation
    end
end

haw = ifft(ha,[],1)*Nt*dt; % isotropic Raman; The factor of Nt is needed and used later in the convolution theorem because of how ifft is defined.
                           %                  The convolution theorem is PQ=F[(p*q)]/N, where (p*q) is the discrete-form convolution, for discrete Fourier transforms.
                           %                  Normal (integral-form) convolution is p*q=(p*q)*dt.
hbw = ifft(hb,[],1)*Nt*dt; % anisotropic Raman

% Incoporate fiber.fr into haw and hbw to save the computational time
haw = fiber.fr*haw;
hbw = fiber.fr*hbw;

% Put the data in GPU if needed
if sim.gpu_yes
    haw = gpuArray(haw);
    hbw = gpuArray(hbw);
end

% If there is no anisotropic Raman contribution, set hbw=[], which will be
% used to determine whether there is an anisotropic Raman effect through
% "isempty(hbw) = true or false" later.
if ~include_anisotropic_Raman
    hbw = [];
end

end