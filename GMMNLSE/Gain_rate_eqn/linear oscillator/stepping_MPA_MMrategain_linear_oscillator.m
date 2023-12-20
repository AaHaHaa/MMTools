function [A1,...
          Power_pump_forward, Power_pump_backward, Power_ASE_forward,...
          N2] = stepping_MPA_MMrategain_linear_oscillator(A0, A0_backward,...
                                                          Power_pump_forward, Power_pump_backward,...
                                                          Power_ASE_forward, Power_ASE_backward,...
                                                          dt, sim, prefactor,...
                                                          SRa_info, SRb_info, SK_info,...
                                                          omegas, D,...
                                                          haw, hbw, sponRS_prefactor,...
                                                          gain_rate_eqn)
%STEPPING_MPA_MMRATEGAIN_LINEAR_OSCILLATOR Take one step with MPA with a
%gain model solved from rate equations. It's only used for multimode
%cases and the gain term is treated as a nonlinear term, instead of a 
%dispersion term in MPA. Because of this, the step sim.deltaZ needs to
%be small enough for MPA to converge.
%
% Input:
%    A0 - initial forward-propagating field in the frequency domain (N, num_modes); sqrt(W)
%    A0_backward - initial backward-propagating field in the frequency domain (N, num_modes); sqrt(W)
%
%    Power_pump_forward - scalar; the power of the co-propagating pump
%    Power_pump_backward - scalar; the power of the counter-propagating pump
%    Power_ASE_forward - (N,num_modes); the power of the co-propagating ASE
%    Power_ASE_backward - (N,num_modes); the power of the counter-propagating ASE
%
%    dt - time grid point spacing; ps
%
%    sim.small_deltaZ - small step size; m
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
%    sim.Raman_model - which Raman model is used
%    sim.Raman_sponRS - consider spontaneous Raman or not
%
%    prefactor - 1i*n2*omega/c; m/W
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
%    omegas - angular frequencies in 1/ps, in the fft ordering
%    D.pos - dispersion term exp(Dz) (N, num_modes, M+1)
%    D.neg - dispersion term exp(-Dz) (N, num_modes, M+1)
%
%    haw - isotropic Raman response in the frequency domain
%    hbw - anisotropic Raman response in the frequency domain
%
%    sponRS_prefactor - prefactor for the spontaneous Raman scattering
%
% *Please refer to "gain_info.m" for details.
%    gain_rate_eqn
%    cross_sections_pump
%    cross_sections
%    overlap_factor - no unit for single-mode and 1/um^2 for multimode
%    N_total - (Nx,Nx); the doped ion density; in "1/um^3"
%    FmFnN - the integral2(overlap_factor*N_total) for the signal and ASE
%    GammaN - the integral2(overlap_factor*N_total) for the pump
%
%    first_backward_before_iterations - 1(true) or 0(false);
%       For bi/counter-pumping cases, the first backward propagation doesn't consider the signal fields.
%
% Output:
%    A1 - the field (in the frequency domain) after one step size (N, num_modes)
%    Power_pump_forward - scalar; the power of the co-propagating pump
%    Power_pump_backward - scalar; the power of the counter-propagating pump
%    Power_ASE_forward - (N,num_modes); the power of the co-propagating ASE
%    N2 - (Nx,Nx); the ion density of the upper state

Power_ASE_forward0 = Power_ASE_forward;

[N,num_modes] = size(A0);

anisotropic_Raman_included = ~sim.scalar & sim.Raman_model==2;

% Spontaneous Raman scattering
if sim.Raman_model ~= 0 && sim.Raman_sponRS
    sponRS = ifft(abs(fft(sponRS_prefactor{1}.*randn(size(sponRS_prefactor{1})).*exp(1i*2*pi*rand(size(sponRS_prefactor{1}))))).^2).*sponRS_prefactor{2};
    sponRS_Gamma = permute(fft(haw.*sponRS),[1,3,2]);
else
    sponRS_Gamma = 0;
end

% Set initial values for psi
psi = repmat(A0, 1, 1, sim.MPA.M+1); % M copies of psi(w,z) = A(w,z), in the frequency domain!

for n_it = 1:sim.MPA.n_tot_max
    % Calculate A(w,z) at all z:
    A_w = D.pos.*psi;

    % Calculate A(t,z) at all z
    A_t = permute(fft(A_w),[1 3 2]); % (N, M+1, num_modes)
    
    % Set up matrices for the following Kerr, Ra, and Rb computations
    if sim.gpu_yes
        Kerr = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Ra = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Rb = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
    else
        Kerr = complex(zeros(N, sim.MPA.M+1, num_modes));
        Ra = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes));
        Rb = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes));
    end
    
    % Calculate large num_modes^4 Kerr, Ra, and Rb terms.
    % If not using the GPU, we will precompute Ra_mn and Rb_mn before the num_modes^4 sum
    if sim.gpu_yes
        % If using the GPU, do the computation with fast CUDA code
        if sim.scalar % scalar fields
            [Kerr,...
             Ra] = feval(sim.cuda_SRSK,...
                         Kerr, Ra,...
                         complex(A_t),...
                         SK_info.SK, SRa_info.SRa,...
                         SRa_info.nonzero_midx1234s,...
                         SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                         sim.Raman_model~=0,...
                         N, sim.MPA.M+1,...
                         num_modes,...
                         sim.cuda_num_operations);
        else % polarized fields
            [Kerr,...
             Ra, Rb] = feval(sim.cuda_SRSK,...
                             Kerr, Ra, Rb,...
                             complex(A_t),...
                             SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                             SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                             SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                             sim.Raman_model~=0, sim.Raman_model==2,...
                             N, sim.MPA.M+1,...
                             num_modes,...
                             sim.cuda_num_operations);
        end
        Kerr = sum(Kerr,4);
    else
        % If using the CPU, first precompute Ra_mn and Rb_mn.
        if sim.Raman_model ~= 0
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
                if anisotropic_Raman_included
                    SRb_nonzero_midx34s = midx34s_sub2ind(SRb_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
                end
            end
            Ra_mn = A_t(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(A_t(:, :, SRa_info.nonzero_midx34s(2,:))); % (N,M+1,num_nonzero34)
            if anisotropic_Raman_included
                Rb_mn = A_t(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(A_t(:, :, SRb_info.nonzero_midx34s(2,:))); % (N,M+1,num_nonzero34)
            end
        end

        % Then calculate Kerr, Ra, and Rb.
        for midx1 = 1:num_modes
            % Kerr
            nz_midx1 = find( SK_info.nonzero_midx1234s(1,:)==midx1 );
            midx2 = SK_info.nonzero_midx1234s(2,nz_midx1);
            midx3 = SK_info.nonzero_midx1234s(3,nz_midx1);
            midx4 = SK_info.nonzero_midx1234s(4,nz_midx1);
            Kerr(:,:,midx1) = sum(permute(SK_info.SK(nz_midx1),[3 2 1]).*A_t(:, :, midx2).*A_t(:, :, midx3).*conj(A_t(:, :, midx4)),3);
            if sim.Raman_model ~= 0
                % Ra
                for midx2 = 1:num_modes
                    nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                    nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                    midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                    midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                    idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                    idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of T_mn
                    Ra(:, :, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[3 2 1]).*Ra_mn(:, :, idx),3);
                end
                % Rb
                if anisotropic_Raman_included
                    for midx2 = 1:num_modes
                        nz_midx1 = find( SRb_info.nonzero_midx1234s(1,:)==midx1 );
                        nz_midx = nz_midx1( SRb_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                        midx3 = SRb_info.nonzero_midx1234s(3,nz_midx);
                        midx4 = SRb_info.nonzero_midx1234s(4,nz_midx);
                        idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                        idx = arrayfun(@(i) find(SRb_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of T_mn
                        Rb(:, :, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[3 2 1]).*Rb_mn(:, :, idx),3);
                    end
                end
            end
        end
        if anisotropic_Raman_included
            clear Ra_mn Rb_mn
        elseif sim.Raman_model ~= 0
            clear Ra_mn
        end
    end
    
    % Calculate the gain term
    [Power_pump_forward,Power_pump_backward,Power_ASE_forward,...
     gain_term,N2] = solve_gain_rate_eqn_linear_oscillator(sim,gain_rate_eqn,...
                                                           A_w,A0_backward,...
                                                           Power_pump_forward,Power_pump_backward,...
                                                           Power_ASE_forward,Power_ASE_backward,...
                                                           omegas,dt);
    
    % Apply the convolution for each part of the Raman sum:
    % The convolution using Fourier Transform is faster if both arrays are
    % large. If one of the array is small, "conv" can be faster.
    % Please refer to
    % "https://blogs.mathworks.com/steve/2009/11/03/the-conv-function-and-implementation-tradeoffs/"
    % for more information.
    if sim.Raman_model ~= 0
        Ra = fft(haw.*ifft(Ra));

        if ~anisotropic_Raman_included
            nonlinear = Kerr + sum(Ra.*permute(A_t,[1 2 4 3]),4) + sponRS_Gamma.*A_t;
        else % polarized fields with an anisotropic Raman
            Rb = fft(hbw.*ifft(Rb));

            nonlinear = Kerr + sum((Ra+Rb).*permute(A_t,[1 2 4 3]),4) + sponRS_Gamma.*A_t;
        end
    else
        nonlinear = Kerr;
    end
    
    % Multiply the nonlinear factor
    nonlinear = prefactor.*permute(ifft(nonlinear),[1 3 2]); % (N,num_modes,M+1)

    % Incorporate deltaZ and D.neg term for the integration
    %zero_idx = (max(abs(A0))==0);
    %gain_term(:,zero_idx,:) = 0; % clear effects from the gain model by making them all zeros
    nonlinear = D.neg.*(sim.small_deltaZ*nonlinear+gain_term); % (N, num_modes, M+1)
	
    % Save the previous psi at the right end, then compute the new psi's
    last_psi = psi(:,:,sim.MPA.M+1);
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
    if sim.gpu_yes
        psi = feval(sim.cuda_MPA_psi_update,psi,...
                                            nonlinear,...
                                            sim.MPA.coeff,...
                                            N,sim.MPA.M+1,num_modes);
    else
        for Midx = 1:sim.MPA.M
            psi(:,:,Midx+1) = psi(:,:,1) + sum(nonlinear(:,:,1:Midx+1).*sim.MPA.coeff(Midx,:,1:Midx+1),3);
        end
    end

    % Calculate the average NRMSE = take the RMSE between the previous psi
    % and the current psi at the right edge, normalize by the absolute max,
    % and average over all modes
    current_psi = psi(:,:,sim.MPA.M+1);
    energy_current_psi = sum(abs(current_psi).^2);
    weight = energy_current_psi/sum(energy_current_psi);
    NRMSE_p = sqrt(sum(abs(current_psi-last_psi).^2)./energy_current_psi).*weight;
    NRMSE_p(isnan(NRMSE_p)) = 0; % exclude modes with all zero fields
    avg_NRMSE = sum(NRMSE_p);
    
    % If it has converged to within tol, then quit
    if (avg_NRMSE < sim.MPA.tol && n_it >= sim.MPA.n_tot_min)
        break
    end
    
    if n_it == sim.MPA.n_tot_max
        error('Error in GMMNLSE_MPA_step: The step did not converge after %d iterations, aborting.', sim.MPA.n_tot_max);
    end
end

% Get back to A from psi
A1 = D.pos(:,:,sim.MPA.M+1).*psi(:,:,sim.MPA.M+1);

% Output powers and doped ion density at the final parallelization plane (M+1)
Power_pump_forward = Power_pump_forward(:,:,:,:,:,end);
Power_pump_backward = Power_pump_backward(:,:,:,:,:,end);
Power_ASE_forward = Power_ASE_forward(:,:,:,:,:,end);
N2 = N2(:,:,:,:,:,end);

% Include ASE power to the signal field
% Because the time window is small, on the order of picoseconds, inclusion
% of ASE power to the field doesn't compute the ASE effect to the gain
% twice besides the Power_ASE_forward term.
% It's continuous light, so the inclusion needs to add a random spectral
% phase.
% If the signal field is zero, we don't need to compute ASE effect to the
% field. All we need is then the Power_ASE_forward.
%
% Unit:
%   Power_ASE_forward: W/THz
%       ASE spectral energy = Power_ASE_forward*(N*dt*1e12): pJ/THz
%   Spectral energy = abs(A1)^2*(N*dt)^2: pJ/THz
if gain_rate_eqn.include_ASE && any(abs(A1(:)))
    dP_ASE = permute(abs(Power_ASE_forward-Power_ASE_forward0),[5,3,1,2,4])/(N*dt);
    A1 = A1 + sqrt(dP_ASE).*exp(1i*2*pi*rand(N,num_modes));
end

end
