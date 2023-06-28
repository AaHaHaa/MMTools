function [D_op,sim] = calc_D_op(fiber,sim,Nt,dt,omegas,fields)
%D_OP It computes the dispersion operator used in GMMNLSE

% The dispersion term in the GMMNLSE, in frequency space
if any(size(fiber.betas) == Nt) % the betas is given over different frequencies,
                                % instead of the coefficients of the Taylor series expansion over the center frequency
    if size(fiber.betas,2) == Nt % betas should be a column vector
        fiber.betas = fiber.betas.';
    end
    if ~isfield(sim,'betas')
        % Obtain the betas of the input pulse
        fftshift_omegas = fftshift(omegas,1);
        spectrum = sum(abs(fftshift(ifft(fields),1)).^2,2);
        omega0 = sum(fftshift_omegas.*spectrum)/sum(spectrum); % 2*pi*THz; the pulse center frequency (under shifted omega)
        omega_range = 1/dt; % 2*pi*THz
        omegas_idx_near_pulse = fftshift_omegas>omega0-omega_range/5 & fftshift_omegas<omega0+omega_range/5;% pick only the data near the pulse center frequency to find its beta0 and beta1
        clear spectrum omega0 omega_range;

        fit_order = 7;
        [betas_Taylor_coeff,~,mu] = polyfit(fftshift_omegas(omegas_idx_near_pulse),real(fiber.betas(omegas_idx_near_pulse,1)),fit_order);
        sim.betas = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
        sim.betas = [sim.betas(1)-sim.betas(2)*mu(1)/mu(2);...
                     sim.betas(2)/mu(2)];
        clear fftshift_omegas fit_order betas_Taylor_coeff mu;
        %{
        % The following code finds sim.betas that minimizes the range and
        % values of betas.
        % However, it can pick a sim.betas that is far away from the input
        % pulse if the center frequency of the frequency window is too far away
        % from the pulse frequency. The result can deviate from the physics if
        % the step size isn't small enough.
        sim.betas = zeros(2,1,'gpuArray');

        % Obtain the betas of the input pulse
        fftshift_omegas = fftshift(omegas,1);
        fit_order = 7;
        [betas_Taylor_coeff,~,mu] = polyfit(fftshift_omegas,real(fiber.betas(:,1)),fit_order);
        sim.betas = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
        new_betas = real(fiber.betas(:,1))-(sim.betas(1)+sim.betas(2)*(fftshift_omegas-mu(1))/mu(2));
        sim.betas = [sim.betas(1)-sim.betas(2)*mu(1)/mu(2) + (max(new_betas) + min(new_betas))/2;...
                     sim.betas(2)/mu(2)];
        clear fftshift_omegas fit_order betas_Taylor_coeff mu new_betas;
        %}
    end
    
    D_op = 1i*(ifftshift(fiber.betas,1)-(sim.betas(1)+sim.betas(2)*omegas));
else
    % D0_op = sum( i*beta_n/n!*omega^n ,n)
    if ~isfield(sim,'betas')
        sim.betas = real(fiber.betas([1 2],1));
    end
    betas = fiber.betas;
    betas([1 2],:) = betas([1 2],:) - sim.betas; % beta0 and beta1 are set relative to sim.betas, or the fundamental mode if sim.betas doesn't exist
    % D_op = sum( 1i*beta_n/n!*omegas^n ,n,0,size(fiber.betas,1)-1 )
    taylor_n = permute(0:size(betas,1)-1,[1 3 2]); % calculation starting here is under the dimension (N,num_modes,order_betas)
    taylor_power = omegas.^taylor_n;
    D_op = sum( 1i./factorial(taylor_n).*permute(betas,[3 2 1]).*taylor_power,3); % sum(...,3) sums over different "n" adding all Taylor series expansion terms
    clear taylor_n taylor_power
end

end

