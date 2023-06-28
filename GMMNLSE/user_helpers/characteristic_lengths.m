function [ dispersion_length,nonlinear_length ] = characteristic_lengths( field_power,t,f0,beta2,Aeff,n2 )
%CHARACTERISTIC_LENGTHS Calculates the characteristic lengths of a pulse
%   
%   field_power: (N,num_modes); the power of the electric field of a pulse (W)
%   t: (N,1); temporal points (ps)
%   f0: (1,num_modes) or a scalar; central frequency (THz)
%   beta2: (1,num_modes); the group dispersion (ps^2/m)
%   Aeff: (1,num_modes); the effective area (m^2)
%   n2: a scalar; nonlinear coefficient (m^2/W)

if nargin < 6
    % the nonlinear coefficient for silica fibers
    n2 = 2.3e-20; % m^2/W
end

num_modes = size(field_power,2);
nonlinear_length = zeros(1,num_modes);
dispersion_length = zeros(1,num_modes);

for n = 1:num_modes

    nonzero_idx = find(field_power(:,n)~=0);
    if length(nonzero_idx) < 3 % almost all elements zero
        % Characteristic lengths
        nonlinear_length(n)  = NaN; % m
        dispersion_length(n) = NaN; % m
    else
        
        %% Extract only the central part of the field to save the computational time
        lower_bound_ratio = 2.2; % It just needs to be larger than 2. Choose 2.2 of no reason.
        field_power( abs(field_power(:,n))<max(field_power(:,n))/lower_bound_ratio ,n) = 0;
        first_nonzero_idx = find(field_power(:,n),1);
        last_nonzero_idx = find(field_power(:,n),1,'last');
        intensity_n = field_power(first_nonzero_idx:last_nonzero_idx,n);
        t_n = t(first_nonzero_idx:last_nonzero_idx);

        if size(t_n,1) == 1
            t_n = t_n';
        end
        
        %% Some parameters
        c = 2.99792458e-4; % speed of ligth m/ps
        if length(f0) ~= 1
            w0 = 2*pi*f0(n);
        else
            w0 = 2*pi*f0; % angular frequency (THz)
        end
        nonlin_const = n2*w0/c; % m/W
        gamma = nonlin_const/Aeff(n); % 1/W/m

        %% Pulse parameters
        P0 = max(intensity_n);
        
        threshold = max(intensity_n)/1.001;
        [~,~,T_fwhm,~] = findpeaks(intensity_n,t_n,'MinPeakHeight',threshold,'WidthReference','halfheight');
        T_fwhm = sum(T_fwhm); % if the field has many small peaks, the T_fwhm is approximately the sum of all T_fwhm's because of the threshold being close to the maximum of the field
        
        T0 = T_fwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665

        % Characteristic lengths
        nonlinear_length(n) = 1/(gamma*P0);     % m
        dispersion_length(n) = T0^2/abs(beta2(n)); % m
        
    end
end

end

