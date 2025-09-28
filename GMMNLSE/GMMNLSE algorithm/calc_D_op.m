function [D_op,sim] = calc_D_op(fiber,sim,Nt,dt,Omega,fields)
%D_OP It computes the dispersion operator used in GMMNLSE

% If the code applies narrowband transformation to the coherent fields,
% fields are downsampled in GMMNLSE_propagate() already. If
% fiber.betas is a function of frequency, it will then contain too many
% points. Its downsampling operation is applied below.
if sim.cs.cs > 1
    if any(size(fiber.betas) >= Nt)
        fiber.betas = fiber.betas(1:sim.cs.cs:end,:);
    end
end

% The dispersion term in the GMMNLSE, in frequency space
if any(size(fiber.betas) == Nt) % the betas is given over different frequencies,
                                % instead of the coefficients of the Taylor series expansion over the center frequency
    if size(fiber.betas,2) == Nt % betas should be a column vector
        fiber.betas = fiber.betas.';
    end
    if ~isfield(sim,'betas')
        if ~any(fields) % fields is all-zero, which typically happens when you want to compute ASE or inversion under rate-eqn gain modeling
            sim.betas = [0;0];
        else
            % Obtain the betas of the input pulse
            fftshift_Omega = fftshift(Omega,1);
            spectrum = sum(abs(fftshift(ifft(fields,[],1),1)).^2,2);
            omega0 = sum(fftshift_Omega.*spectrum)/sum(spectrum); % 2*pi*THz; the pulse center frequency (under shifted Omega)
            Omega_range = 2*pi/dt; % 2*pi*THz
            Omega_idx_near_pulse = fftshift_Omega>omega0-Omega_range/5 & fftshift_Omega<omega0+Omega_range/5;% pick only the data near the pulse center frequency to find its beta0 and beta1
            clear spectrum omega0 Omega_range;
    
            fit_order = max(2,min(7,sum(Omega_idx_near_pulse)-1)); % 2~7
            [betas_Taylor_coeff,~,mu] = polyfit(fftshift_Omega(Omega_idx_near_pulse),real(fiber.betas(Omega_idx_near_pulse,1)),fit_order);
            sim.betas = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
            sim.betas = [sim.betas(1)-sim.betas(2)*mu(1)/mu(2);...
                         sim.betas(2)/mu(2)];
            clear fftshift_Omega fit_order betas_Taylor_coeff mu;
            %{
            % The following code finds sim.betas that minimizes the range and
            % values of betas.
            % However, it can pick a sim.betas that is far away from the input
            % pulse if the center frequency of the frequency window is too far away
            % from the pulse frequency. The result can deviate from the physics if
            % the step size isn't small enough.
            sim.betas = zeros(2,1,'gpuArray');
    
            % Obtain the betas of the input pulse
            fftshift_Omega = fftshift(Omega,1);
            fit_order = 7;
            [betas_Taylor_coeff,~,mu] = polyfit(fftshift_Omega,real(fiber.betas(:,1)),fit_order);
            sim.betas = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
            new_betas = real(fiber.betas(:,1))-(sim.betas(1)+sim.betas(2)*(fftshift_Omega-mu(1))/mu(2));
            sim.betas = [sim.betas(1)-sim.betas(2)*mu(1)/mu(2) + (max(new_betas) + min(new_betas))/2;...
                         sim.betas(2)/mu(2)];
            clear fftshift_Omega fit_order betas_Taylor_coeff mu new_betas;
            %}
        end
    end
    
    D_op = 1i*(ifftshift(fiber.betas,1)-(sim.betas(1)+sim.betas(2)*Omega));
elseif size(fiber.betas,1) < 20
    % D0_op = sum( i*beta_n/n!*omega^n ,n)
    if ~isfield(sim,'betas')
        sim.betas = real(fiber.betas([1 2],1));
    end
    betas = fiber.betas;
    betas([1 2],:) = betas([1 2],:) - sim.betas; % beta0 and beta1 are set relative to sim.betas, or the fundamental mode if sim.betas doesn't exist
    % D_op = sum( 1i*beta_n/n!*Omega^n ,n,0,size(fiber.betas,1)-1 )
    taylor_n = permute(0:size(betas,1)-1,[1 3 2]); % calculation starting here is under the dimension (N,num_modes,order_betas)
    taylor_power = Omega.^taylor_n;
    D_op = sum( 1i./factorial(taylor_n).*permute(betas,[3 2 1]).*taylor_power,3); % sum(...,3) sums over different "n" adding all Taylor series expansion terms
    clear taylor_n taylor_power
else
    error('GMMNLSE_propagate:betasError',...
          '"fiber.betas" can only have the size of (Nt,num_modes) or (num_Taylor_series_coeff,num_modes),\nwhere num_Taylor_series_coeff<20.');
end

% Scaled dispersion due to the narrowband transformation (scaled Fourier transform)
if sim.cs.cs > 1
    D_op = D_op/sim.cs.cs;
end

end