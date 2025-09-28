function [A1w, a5,...
          Power_pump_forward, Power_pump_backward,...
          opt_dz, success,...
          N] = stepping_RK4IP_rategain_adaptive(A0w, N0,...
                                                dt,...
                                                sim, gain_rate_eqn,...
                                                SK_info, SRa_info, SRb_info,...
                                                haw, hbw,...
                                                At_noise,...
                                                n2_prefactor, Omega, D_op,...
                                                Power_pump_forward, Power_pump_backward, a5_1,...
                                                dummy_var)
%STEPPING_RK4IP_RATEGAIN_ADAPTIVE Take one step with RK4IP with a gain model
%solved from rate equations. The gain term is treated as a dispersion term,
%instead of a nonlinear term.
%
% Input:
%    A0w - initial forward-propagating field in the frequency domain (Nt, num_modes); sqrt(W)
%    dt - time grid point spacing; ps
%
%    sim.dz - small step size; m
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
%    Omega - offset angular frequencies; in 1/ps, in the fft ordering
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
%
% For adaptive-step implementation, check http://www.sciencedirect.com/science/article/pii/S0010465512004262
%
%    Balac and Mahe, Embedded Runge-Kutta scheme for step-size control in 
%    the interaction picture method, "Comput. Phys. Commun. 184(4), 1211-
%    1219 (2013)

[Nt,num_modes] = size(A0w);

% Set up matrices for the following Kerr, Ra, and Rb computations
if sim.gpu_yes
    Kerr = complex(zeros(Nt, num_modes, num_modes, 'gpuArray'));
    Ra = complex(zeros(Nt, num_modes, num_modes, 'gpuArray'));
    Rb = complex(zeros(Nt, num_modes, num_modes, 'gpuArray'));
else
    Kerr = complex(zeros(Nt, num_modes));
    Ra = complex(zeros(Nt, num_modes, num_modes));
    Rb = complex(zeros(Nt, num_modes, num_modes));
end

% Represented under the interaction picture (dispersion + gain)
Dz = D_op*sim.dz/2;
if gain_rate_eqn.counterpump_power == 0 % copumping
    [Power_pump_forward,~,~,...
     G,N] = solve_gain_rate_eqn('forward',...
                                sim,gain_rate_eqn,...
                                N0,...
                                A0w,dummy_var,...
                                Power_pump_forward,dummy_var,...
                                dummy_var,dummy_var,...
                                Omega,dt,...
                                false);
    Power_pump_backward = 0;
else % bi-pumping or counterpumping
    [Power_pump_forward,Power_pump_backward,~,~,...
     G,N] = solve_gain_rate_eqn_linear_oscillator(sim,gain_rate_eqn,...
                                                  N0,...
                                                  A0w,dummy_var,...
                                                  Power_pump_forward,Power_pump_backward,...
                                                  dummy_var,dummy_var,...
                                                  Omega,dt);
end
gz_over_2 = log(G)/2; % E = G*E0 = exp(g*dz)*E0. The gz_over_2 here has already been multiplied by "dz/2" (symmetrized SS).
expDG = exp(Dz + gz_over_2);

A_IP = expDG.*A0w;

% Propagate through the nonlinearity
if isempty(a5_1)
    a5_1 = N_op(        A0w,...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                At_noise,...
                n2_prefactor,...
                Nt, num_modes);
end
a1 = expDG.*a5_1;
a2 =       N_op(        A_IP+a1*(sim.dz/2),...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                At_noise,...
                n2_prefactor,...
                Nt, num_modes);
a3 =       N_op(        A_IP+a2*(sim.dz/2),...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                At_noise,...
                n2_prefactor,...
                Nt, num_modes);
a4 =       N_op(expDG.*(A_IP+a3*(sim.dz)),...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                At_noise,...
                n2_prefactor,...
                Nt, num_modes);

A1w = expDG.*(A_IP + (a1+2*a2+2*a3)*(sim.dz/6)) + a4*(sim.dz/6);

% Local error estimate
a5 =       N_op(        A1w,...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                At_noise,...
                n2_prefactor,...
                Nt, num_modes);
err = sum(abs((a4-a5)*(sim.dz/10)).^2,1);

% Stepsize control
normA = sum(abs(A1w).^2,1);
err = sqrt(err./normA);
err = max(err(normA~=0));
if normA == 0 % all-zero field; this will make err empty, so this condition needs to be determined first
    opt_dz = 2*sim.dz;
    success = true;
elseif isnan(err) % the computation is just so wrong, so we reduce the step size and do it again
    opt_dz = 0.5*sim.dz;
    success = false;
else
    opt_dz = max(0.5,min(2,0.8*(sim.adaptive_dz.threshold/err)^(1/4)))*sim.dz;

    success = err < sim.adaptive_dz.threshold;
end

end

function dAdz = N_op(Aw,...
                     sim,...
                     SK_info, SRa_info, SRb_info,...
                     Kerr, Ra, Rb,...
                     haw, hbw,...
                     At_noise,...
                     n2_prefactor,...
                     Nt, num_modes)
%N_op Calculate dAdz

At = fft(Aw,[],1);
At_wNoise = At + At_noise;

% Calculate large num_modes^4 Kerr, Ra, and Rb terms.
% If not using the GPU, we will precompute Ra_mn and Rb_mn before the num_modes^4 sum
if sim.gpu_yes
    % If using the GPU, do the computation with fast CUDA code
    if sim.scalar % scalar fields
        [Kerr,...
         Ra] = feval(sim.cuda_SRSK,...
                     Kerr, Ra,...
                     complex(At_wNoise),...
                     SK_info.SK, SRa_info.SRa,...
                     SRa_info.nonzero_midx1234s,...
                     SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                     sim.include_Raman,...
                     int32(Nt), 1,...
                     num_modes,...
                     sim.cuda_num_operations_SRSK);
    else % polarized fields
        [Kerr,...
         Ra, Rb] = feval(sim.cuda_SRSK,...
                         Kerr, Ra, Rb,...
                         complex(At_wNoise),...
                         SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                         SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                         SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                         sim.include_Raman, ~isempty(hbw),...
                         int32(Nt), 1,...
                         num_modes,...
                         sim.cuda_num_operations_SRSK);
    end
    Kerr = sum(Kerr,3);
else
    % If using the CPU, first precompute Ra_mn and Rb_mn.
    if sim.include_Raman
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
        Ra_mn = At_wNoise(:, SRa_info.nonzero_midx34s(1,:)).*conj(At_wNoise(:, SRa_info.nonzero_midx34s(2,:))); % (Nt,num_nonzero34)
        if ~isempty(hbw)
            SRb_nonzero_midx34s = midx34s_sub2ind(SRb_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
            Rb_mn = At_wNoise(:, SRb_info.nonzero_midx34s(1,:)).*conj(At_wNoise(:, SRb_info.nonzero_midx34s(2,:))); % (Nt,num_nonzero34)
        end
    end
    
    % Then calculate Kerr, Ra, and Rb.
    for midx1 = 1:num_modes
        % Kerr
        nz_midx1 = find( SK_info.nonzero_midx1234s(1,:)==midx1 );
        midx2 = SK_info.nonzero_midx1234s(2,nz_midx1);
        midx3 = SK_info.nonzero_midx1234s(3,nz_midx1);
        midx4 = SK_info.nonzero_midx1234s(4,nz_midx1);
        Kerr(:,midx1) = sum(permute(SK_info.SK(nz_midx1),[2 1]).*At_wNoise(:, midx2).*At_wNoise(:, midx3).*conj(At_wNoise(:, midx4)),2);
        if sim.include_Raman
            % Ra
            for midx2 = 1:num_modes
                nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 2nd-dimensional "num_nonzero34" of Ra_mn
                Ra(:, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[2 1]).*Ra_mn(:, idx),2);
            end
            % Rb
            if ~isempty(hbw)
                for midx2 = 1:num_modes
                    nz_midx1 = find( SRb_info.nonzero_midx1234s(1,:)==midx1 );
                    nz_midx = nz_midx1( SRb_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                    midx3 = SRb_info.nonzero_midx1234s(3,nz_midx);
                    midx4 = SRb_info.nonzero_midx1234s(4,nz_midx);
                    idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                    idx = arrayfun(@(i) find(SRb_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of Rb_mn
                    Rb(:, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[2 1]).*Rb_mn(:, idx),2);
                end
            end
        end
    end
    if sim.include_Raman
        clearvars Ra_mn;
        if ~isempty(hbw)
            clearvars Rb_mn;
        end
    end
end

% Calculate h*Ra as F-1(h F(Ra))
% The convolution using Fourier transform is faster if both arrays are
% large. If one of the array is small, "conv" can be faster.
% Please refer to
% "https://blogs.mathworks.com/steve/2009/11/03/the-conv-function-and-implementation-tradeoffs/"
% for more information.
if sim.include_Raman
    Ra = fft(haw.*ifft(Ra,[],1),[],1);
    
    if ~isempty(hbw) % polarized fields with an anisotropic Raman
        Rb = fft(hbw.*ifft(Rb,[],1),[],1);
    end
    
    if isempty(hbw)
        nonlinear = Kerr + sum(Ra.*permute(At_wNoise,[1 3 2]),3);
    else % polarized fields with an anisotropic Raman
        nonlinear = Kerr + sum((Ra+Rb).*permute(At_wNoise,[1 3 2]),3);
    end
else
    nonlinear = Kerr;
end

% Now everything has been summed into "nonlinear", so transform into the
% frequency domain for the n2_prefactor
dAdz = n2_prefactor.*ifft(nonlinear,[],1);

end