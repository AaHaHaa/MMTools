function [A1w, dummy_a5,...
          Power_pump_forward, Power_pump_backward,...
          opt_dz, success,...
          N] = stepping_MPA_rategain_adaptive(A0w, N0,...
                                              dt,...
                                              sim, gain_rate_eqn,...
                                              SK_info, SRa_info, SRb_info,...
                                              haw, hbw,...
                                              sponRS_prefactor,...
                                              n2_prefactor, omegas, D_op,...
                                              Power_pump_forward, Power_pump_backward, dummy_a5_1,...
                                              dummy_var)
%STEPPING_MPA_RATEGAIN_ADAPTIVE Take one step with MPA with a gain model 
%solved from rate equations. It's only used for multimode cases and the 
%gain term is treated as a nonlinear term, instead of a dispersion term. 
%Because of this, the step sim.dz needs to be small enough for MPA to 
%converge.
%
% Input:
%    A0w - initial forward-propagating field in the frequency domain (Nt, num_modes); sqrt(W)
%    dt - time grid point spacing; ps
%
%    sim.small_dz - small step size; m
%
%    sim.MPA.M - parallel extent, 1 is no parallelization
%    sim.MPA.n_tot_max - maximum number of iterations
%    sim.MPA.n_tot_min - minimum number of iterations
%    sim.MPA.tol - tolerance for convergence at each iteration
%
%    sim.scalar - scalar or polarized fields
%    sim.gpu_yes - true = GPU, false = CPU
%
%    sim.cuda_SRSK - the cuda for computing SR and SK values
%
%    sim.include_Raman - whether Raman is included
%
%    gain_rate_eqn - container of rate-eqn-gain parameters
%
%    haw - isotropic Raman response in the frequency domain
%    hbw - anisotropic Raman response in the frequency domain
%
%    sponRS_prefactor - prefactor for the spontaneous Raman scattering
%
%    SRa_info.SRa - SRa tensor; m^-2
%    SRa_info.nonzero_midx1234s - required SRa indices in total
%    SRa_info.nonzero_midx34s - required (SRa) indices for partial Raman term (only for CPU computation)
%    SRb_info.SRb - SRb tensor; m^-2
%    SRb_info.nonzero_midx1234s - required SRb indices in total
%    SRb_info.nonzero_midx34s - required (SRb) indices for partial Raman term (only for CPU computation)
%    SK_info.SK - SK tensor; m^2 (unempty if considering polarizaton modes)
%    SK_info.nonzero_midx1234s - required SK indices in total (unempty if considering polarizaton modes)
%
%    n2_prefactor - 1i*n2*omega/c; m/W
%    omegas - angular frequencies in 1/ps, in the fft ordering
%    D_op - dispersion term D (Nt, num_modes)
%
%    cross_sections_pump
%    cross_sections
%    overlap_factor - no unit for single-mode and 1/um^2 for multimode
%    N_total - (Nx,Nx); the doped ion density; in "1/um^3"
%    FmFnN - the integral2(overlap_factor*N_total) for the signal and ASE
%    GammaN - the integral2(overlap_factor*N_total) for the pump
%
%    Power_pump_forward - scalar; the power of the co-propagating pump
%    Power_pump_backward - scalar; the power of the counter-propagating pump
%
%    dummy_var - unused variable
%
% Output:
%    A1w - the field (in the frequency domain) after one step size (Nt, num_modes)
%    a5 - the RK4 term that can be reused in the next step
%    Power_pump_forward - scalar; the power of the co-propagating pump
%    Power_pump_backward - scalar; the power of the counter-propagating pump
%    opt_dz - recommended step size
%    success - whether the current step size is sufficiently small for the required tolerance
%    N - (Nx,Nx); the ion density/population of each energy level

sim.small_dz = sim.dz/sim.MPA.M; % the step size between each parallelization plane in MPA

[Nt,num_modes] = size(A0w);

% Spontaneous Raman scattering
if sim.gpu_yes
    At_noise = fft(sponRS_prefactor{1}.*sqrt(abs(randn(Nt,sim.MPA.M+1,num_modes,'gpuArray'))).*exp(1i*2*pi*rand(Nt,sim.MPA.M+1,num_modes,'gpuArray')));
else
    At_noise = fft(sponRS_prefactor{1}.*sqrt(abs(randn(Nt,sim.MPA.M+1,num_modes))).*exp(1i*2*pi*rand(Nt,sim.MPA.M+1,num_modes)));
end

% We can pre-compute exp(D_op*z) and exp(-D_op*z) for all z.
z = permute((0:sim.MPA.M)',[2,3,1]);% to make "Dz" into (Nt,num_modes,M+1)
Dz = D_op*sim.small_dz.*z;
D = struct('pos', exp(Dz), 'neg', exp(-Dz));

% Set initial values for psi
psi = repmat(A0w, 1, 1, sim.MPA.M+1); % M copies of psi(w,z) = A(w,z), in the frequency domain!
current_psi_half_parallelization = psi(:,:,sim.MPA.M+1);

for n_it = 1:sim.MPA.n_tot_max
    % Calculate A(w,z) at all z:
    Aw = D.pos.*psi;

    % Calculate A(t,z) at all z
    At = permute(fft(Aw),[1 3 2]); % (Nt, M+1, num_modes)
    At_wNoise = At + At_noise;

    % Set up matrices for the following Kerr, Ra, and Rb computations
    if sim.gpu_yes
        Kerr = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Ra = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Rb = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));

        Ra_sponRS = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes, 'gpuArray')); % spontaneous isotropic Raman scattering
        Rb_sponRS = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes, 'gpuArray')); % spontaneous anisotropic Raman scattering
    else
        Kerr = complex(zeros(Nt, sim.MPA.M+1, num_modes));
        Ra = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes));
        Rb = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes));

        Ra_sponRS = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes)); % spontaneous isotropic Raman scattering
        Rb_sponRS = complex(zeros(Nt, sim.MPA.M+1, num_modes, num_modes)); % spontaneous anisotropic Raman scattering
    end
    
    % Calculate large num_modes^4 Kerr, Ra, and Rb terms.
    % If not using the GPU, we will precompute Ra_mn and Rb_mn before the num_modes^4 sum
    if sim.gpu_yes
        % If using the GPU, do the computation with fast CUDA code
        if sim.scalar % scalar fields
            [Kerr,...
             Ra,...
             Ra_sponRS] = feval(sim.cuda_SRSK,...
                                Kerr, Ra, Ra_sponRS,...
                                complex(At), At_noise,...
                                SK_info.SK, SRa_info.SRa,...
                                SRa_info.nonzero_midx1234s,...
                                SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                sim.include_Raman,...
                                int32(Nt), sim.MPA.M+1,...
                                num_modes,...
                                sim.cuda_num_operations_SRSK);
        else % polarized fields
            [Kerr,...
             Ra, Rb,...
             Ra_sponRS, Rb_sponRS] = feval(sim.cuda_SRSK,...
                                           Kerr, Ra, Rb, Ra_sponRS, Rb_sponRS,...
                                           complex(At), At_noise,...
                                           SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                                           SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                           SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                                           sim.include_Raman, ~isempty(hbw),...
                                           int32(Nt), sim.MPA.M+1,...
                                           num_modes,...
                                           sim.cuda_num_operations_SRSK);
        end
        Kerr = sum(Kerr,4);
    else
        % If using the CPU, first precompute Ra_mn and Rb_mn.
        if sim.include_Raman
            if n_it == 1
                midx34s_sub2ind = @(x)...
                    cellfun(@(xx)...
                        feval(@(sub) sub2ind(num_modes*ones(1,2),sub{:}), num2cell(xx)),... % this extra "feval" is to get "xx", which is of the size 2x1, into the input arguments of "sub2ind", so transforming "xx" into a 2x1 cell, each containing an integer, and using {:} expansion is necessary
                    mat2cell(x,2,ones(1,size(x,2)))); % transform (2,num_nonzero34) midx34s into linear indices of a num_modes-by-num_modes matrix
                    % What "midx34s_sub2ind" does (e.g.):
                    %
                    %   x = [1 3;
                    %        5 4]
                    %
                    %   After "mat2cell": {[1;  {[3;  (2x1 cells, each having 2x1 array)
                    %                       5]}   4]}
                    %
                    %   First,
                    %
                    %   xx = {[1;  , then after "num2cell": {{1}; (1 cell with 2x1 cell)
                    %          5]}                           {5}}
                    %
                    %   The purpose of separating 1 and 5 into cells is to use
                    %   index expansion, {:}, to put them into the input
                    %   arguments of "sub2ind" function.
                    %
                    %   For 6 modes and thus for 6x6 matrix, sub2ind([6 6],1,5) = 25
                    %
                    %   Do the same for xx = {[3;  and get sub2ind([6 6],3,4) = 21
                    %                          4]}
                    %   Finally, midx34s_sub2ind = [25 21] (1x2 array)

                SRa_nonzero_midx34s = midx34s_sub2ind(SRa_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
                if ~isempty(hbw)
                    SRb_nonzero_midx34s = midx34s_sub2ind(SRb_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
                end
            end
            Ra_mn = At(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(At(:, :, SRa_info.nonzero_midx34s(2,:))); % (Nt,M+1,num_nonzero34)
            if ~isempty(hbw)
                Rb_mn = At(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(At(:, :, SRb_info.nonzero_midx34s(2,:))); % (Nt,M+1,num_nonzero34)
            end
            
            % spontaneous Raman scattering
            Ra_mn_sponRS =      At(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(At_noise(:, :, SRa_info.nonzero_midx34s(2,:))) +... % (Nt,M+1,num_nonzero34)
                           At_noise(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(     At(:, :, SRa_info.nonzero_midx34s(2,:))) +...
                           At_noise(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(At_noise(:, :, SRa_info.nonzero_midx34s(2,:)));
            if ~isempty(hbw)
                Rb_mn_sponRS =      At(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(At_noise(:, :, SRb_info.nonzero_midx34s(2,:))) +... % (Nt,M+1,num_nonzero34)
                               At_noise(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(     At(:, :, SRb_info.nonzero_midx34s(2,:))) +...
                               At_noise(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(At_noise(:, :, SRb_info.nonzero_midx34s(2,:)));
            end
        end

        % Then calculate Kerr,Ra,Rb.
        for midx1 = 1:num_modes
            % Kerr
            nz_midx1 = find( SK_info.nonzero_midx1234s(1,:)==midx1 );
            midx2 = SK_info.nonzero_midx1234s(2,nz_midx1);
            midx3 = SK_info.nonzero_midx1234s(3,nz_midx1);
            midx4 = SK_info.nonzero_midx1234s(4,nz_midx1);
            Kerr(:,:,midx1) = sum(permute(SK_info.SK(nz_midx1),[3 2 1]).*At_wNoise(:, :, midx2).*At_wNoise(:, :, midx3).*conj(At_wNoise(:, :, midx4)),3);
            if sim.include_Raman
                % Ra
                for midx2 = 1:num_modes
                    nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                    nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                    midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                    midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                    idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                    idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of T_mn
                    Ra(:, :, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[3 2 1]).*Ra_mn(:, :, idx),3);
                    Ra_sponRS(:, :, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[3 2 1]).*Ra_mn_sponRS(:, :, idx),3);
                end
                % Rb
                if ~isempty(hbw)
                    for midx2 = 1:num_modes
                        nz_midx1 = find( SRb_info.nonzero_midx1234s(1,:)==midx1 );
                        nz_midx = nz_midx1( SRb_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                        midx3 = SRb_info.nonzero_midx1234s(3,nz_midx);
                        midx4 = SRb_info.nonzero_midx1234s(4,nz_midx);
                        idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                        idx = arrayfun(@(i) find(SRb_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of T_mn
                        Rb(:, :, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[3 2 1]).*Rb_mn(:, :, idx),3);
                        Rb_sponRS(:, :, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[3 2 1]).*Rb_mn_sponRS(:, :, idx),3);
                    end
                end
            end
        end
        if sim.include_Raman
            clearvars Ra_mn Ra_mn_sponRS
            if ~isempty(hbw)
                clearvars Rb_mn Rb_mn_sponRS;
            end
        end
    end
    
    % Calculate the gain term
    if gain_rate_eqn.counterpump_power == 0 % copumping
        [Power_pump_forward,~,~,...
         gAwdz,N] = solve_gain_rate_eqn('forward',...
                                        sim,gain_rate_eqn,...
                                        N0,...
                                        Aw,dummy_var,...
                                        Power_pump_forward,dummy_var,...
                                        dummy_var,dummy_var,...
                                        omegas,dt,...
                                        false);
        Power_pump_backward = 0;
    else % bi-pumping or counterpumping
        [Power_pump_forward,Power_pump_backward,~,~,...
         gAwdz,N] = solve_gain_rate_eqn_linear_oscillator(sim,gain_rate_eqn,...
                                                          N0,...
                                                          Aw,dummy_var,...
                                                          Power_pump_forward,Power_pump_backward,...
                                                          dummy_var,dummy_var,...
                                                          omegas,dt);
    end
    
    % Apply the convolution for each part of the Raman sum:
    % The convolution using Fourier Transform is faster if both arrays are
    % large. If one of the array is small, "conv" can be faster.
    % Please refer to
    % "https://blogs.mathworks.com/steve/2009/11/03/the-conv-function-and-implementation-tradeoffs/"
    % for more information.
    if sim.include_Raman
        Ra = fft(haw.*ifft(Ra));
        Ra_sponRS = fft(haw.*ifft(Ra_sponRS).*sponRS_prefactor{2});
        
        if ~isempty(hbw) % polarized fields with an anisotropic Raman
            Rb = fft(hbw.*ifft(Rb));
            Rb_sponRS = fft(hbw.*ifft(Rb_sponRS).*sponRS_prefactor{2});
        end
        
        if isempty(hbw)
            nonlinear = Kerr + sum((Ra+Ra_sponRS).*permute(At,[1 2 4 3]),4);
        else % polarized fields with an anisotropic Raman
            nonlinear = Kerr + sum((Ra+Rb+Ra_sponRS+Rb_sponRS).*permute(At,[1 2 4 3]),4);
        end
    else
        nonlinear = Kerr;
    end
    
    % Multiply the nonlinear factor
    nonlinear = n2_prefactor.*permute(ifft(nonlinear),[1 3 2]); % (Nt,num_modes,M+1)
    
    % Incorporate dz and D.neg term for the integration
    %zero_idx = (max(abs(A0))==0);
    %gain_term(:,zero_idx,:) = 0; % clear effects from the gain model by making them all zeros
    nonlinear = D.neg.*(sim.small_dz*nonlinear + gAwdz); % (Nt, num_modes, M+1)
	
    % Save the previous psi at the right end, then compute the new psi's
    last_psi_half_parallelization = current_psi_half_parallelization;
    % ---------------------------------------------------------------------
    % "psi" is calculated from trapezoidal integrals, which results in
    % (1/2,1,1,......,1,1/2) coefficients.
    % Note that for each z plane, it has a (1/2, 1, 1...,1/2) factor in "nonlinear".
    %nonlinear(:,:,1) = nonlinear(:,:,1)/2;
    %psi(:,:,2:end) = psi(:,:,1) + cumsum(nonlinear(:,:,1:sim.MPA.M),3) + nonlinear(:,:,2:end)/2;
    % ---------------------------------------------------------------------
    % "psi" is calculated with a multistep Adams-Moulton method for each parallelization plane.
    % A M-step Adams-Moulton method (with the summation of M's derivatives) has a (M+1)th-order accuracy.
    % With it, the latter psi will not have a bigger error.
    % The error builds up for latter parallization planes in the typical trapezoidal quadrature.
    % The highest-order accuracy is 2nd order O(h^2) from
    %   y(n+1)=y(n)+1/2*h*( f(t(n+1),y(n+1)) + f(t(n),y(n)) ),
    % so the error estimate of adaptive-step-size control uses cubic root: error = O(h^3).
    if sim.gpu_yes
        psi = feval(sim.cuda_MPA_psi_update,psi,...
                                            nonlinear,...
                                            sim.MPA.coeff,...
                                            Nt,sim.MPA.M+1,num_modes);
    else
        for Midx = 1:sim.MPA.M
            psi(:,:,Midx+1) = psi(:,:,1) + sum(nonlinear(:,:,1:Midx+1).*sim.MPA.coeff(Midx,:,1:Midx+1),3);
        end
    end

    % Calculate the average NRMSE = take the RMSE between the previous psi
    % and the current psi at the right edge, normalize by the absolute max,
    % and average over all modes
    current_psi_half_parallelization = psi(:,:,1) + 2*sum(nonlinear(:,:,1:2:end).*sim.MPA.coeff(sim.MPA.M/2,:,1:sim.MPA.M/2+1),3);
    energy_current_psi_half_parallelization = sum(abs(current_psi_half_parallelization).^2);
    weight = energy_current_psi_half_parallelization/sum(energy_current_psi_half_parallelization);
    NRMSE_p = sqrt(sum(abs(current_psi_half_parallelization-last_psi_half_parallelization).^2)./energy_current_psi_half_parallelization).*weight;
    NRMSE_p(isnan(NRMSE_p)) = 0; % exclude modes with all zero fields
    avg_NRMSE = sum(NRMSE_p);
    
    % If it has converged to within tol, then quit
    if avg_NRMSE < sim.MPA.tol && n_it >= sim.MPA.n_tot_min
        break
    end
    
    if n_it == sim.MPA.n_tot_max
        error('Error in GMMNLSE_MPA_step: The step did not converge after %d iterations, aborting.', sim.MPA.n_tot_max);
    end
end

% Get back to A from psi
psi_full_parallelization = psi(:,:,sim.MPA.M+1);
A1w = D.pos(:,:,sim.MPA.M+1).*psi_full_parallelization;

% Output powers and doped ion density at the final parallelization plane (M+1)
Power_pump_forward = Power_pump_forward(:,:,:,:,:,end);
N = N(:,:,:,:,:,end,:,:);

% Stepsize control
err = sqrt(sum(abs(current_psi_half_parallelization-psi_full_parallelization).^2)./energy_current_psi_half_parallelization).*weight;
err(isnan(err)) = 0; % exclude modes with all zero fields
err = sum(err);
if isnan(err) % the computation is just so wrong, so we reduce the step size and do it again
    opt_dz = 0.5*sim.dz;
    success = false;
else
    opt_dz = max(0.5,min(2,0.8*(sim.adaptive_dz.threshold/err)^(1/3)))*sim.dz; % limit the optimal step size to 0.5~2 times as the current step size

    success = err < sim.adaptive_dz.threshold;
end

dummy_a5 = dummy_a5_1; % unused variable

end