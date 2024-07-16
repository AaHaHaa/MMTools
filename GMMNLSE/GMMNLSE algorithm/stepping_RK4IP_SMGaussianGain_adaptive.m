function [A1,a5,...
          opt_dz,success] = stepping_RK4IP_SMGaussianGain_adaptive(A0, dt,...
                                                                   sim, prefactor,...
                                                                   SK_info, SRa_info, SRb_info,...
                                                                   D_op,...
                                                                   haw, hbw,...
                                                                   haw_sponRS, hbw_sponRS, sponRS_prefactor,...
                                                                   a5_1,...
                                                                   G, saturation_energy)
%STEPPING_RK4IP_SMAGUSSIANGAIN_ADAPTIVE Take one step with RK4IP with a 
%total-energy saturating gain model
%
% Input:
%    A0 - initial field in the frequency domain (N, num_modes); sqrt(W)
%    dummy_dt - time grid point spacing; ps; but this isn't used in this function
%
%    sim.scalar - scalar or polarized fields
%    sim.gpu_yes - true = GPU, false = CPU
%
%    sim.cuda_SRSK - the cuda for computing SR and SK values
%
%    sim.Raman_model - which Raman model is used
%    sim.Raman_sponRS - consider spontaneous Raman or not
%
%    sim.adaptive_dz.threshold - the accuracy used to determined whether to increase or decrease the step size
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
%    D_op - dispersion term Dz/2 (N, num_modes)
%    haw - isotropic Raman response in the frequency domain
%    hbw - anisotropic Raman response in the frequency domain
%
%    sponRS_prefactor - prefactor for the spontaneous Raman scattering
%
%    a5_1 - the previous RK4 term that can be reused
%
%    dummy_G - gain term (N, 1); but this isn't used in this function
%    dummy_saturation_intensity - saturation intensity; J/m^2; but this isn't used in this function
%    dummy_fr - Raman fraction; but this isn't used in this function
%
% Output:
%    A1 - the field (in the frequency domain) after one step size (N, num_modes)
%    a5 - the RK4 term that can be reused in the next step
%    opt_dz - recommended step size
%    success - whether the current step size is sufficiently small for the required tolerance

[N,num_modes] = size(A0);

anisotropic_Raman_included = ~sim.scalar & sim.Raman_model==2;

% Set up matrices for the following Kerr, Ra, and Rb computations
if sim.gpu_yes
    Kerr = complex(zeros(N, num_modes, num_modes, 'gpuArray'));
    Ra = complex(zeros(N, num_modes, num_modes, 'gpuArray'));
    Rb = complex(zeros(N, num_modes, num_modes, 'gpuArray'));
    if sim.include_sponRS % spontaneous Raman scattering
        Ra_sponRS = complex(zeros(N, num_modes, num_modes, 'gpuArray'));
        Rb_sponRS = complex(zeros(N, num_modes, num_modes, 'gpuArray'));
        
        A_sponRS = fft(sponRS_prefactor{1}.*sqrt(abs(randn(N,num_modes,'gpuArray'))).*exp(1i*2*pi*rand(N,num_modes,'gpuArray')));
    else
        Ra_sponRS = [];
        Rb_sponRS = [];
        A_sponRS = [];
    end
else
    Kerr = complex(zeros(N, num_modes));
    Ra = complex(zeros(N, num_modes, num_modes));
    Rb = complex(zeros(N, num_modes, num_modes));
    if sim.include_sponRS % spontaneous Raman scattering
        Ra_sponRS = complex(zeros(N, num_modes, num_modes));
        Rb_sponRS = complex(zeros(N, num_modes, num_modes));
        
        A_sponRS = fft(sponRS_prefactor{1}.*sqrt(abs(randn(N,num_modes))).*exp(1i*2*pi*rand(N,num_modes)));
    else
        Ra_sponRS = [];
        Rb_sponRS = [];
        A_sponRS = [];
    end
end

% Represented under the interaction picture (dispersion + gain)
total_energy = calc_total_energy(A0,N,dt); % in nJ
expDG = exp((D_op + G/(1+total_energy/saturation_energy))*sim.dz/2);
A_IP = expDG.*A0;

% Propagate through the nonlinearity
if isempty(a5_1)
    a5_1 = N_op(       A0,...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                haw_sponRS, hbw_sponRS,...
                A_sponRS, sponRS_prefactor,...
                Ra_sponRS, Rb_sponRS,...
                anisotropic_Raman_included,...
                prefactor,...
                N, num_modes);
end
a1 = expDG.*a5_1;
a2 =       N_op(        A_IP+a1*(sim.dz/2),...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                haw_sponRS, hbw_sponRS,...
                A_sponRS, sponRS_prefactor,...
                Ra_sponRS, Rb_sponRS,...
                anisotropic_Raman_included,...
                prefactor,...
                N, num_modes);
a3 =       N_op(        A_IP+a2*(sim.dz/2),...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                haw_sponRS, hbw_sponRS,...
                A_sponRS, sponRS_prefactor,...
                Ra_sponRS, Rb_sponRS,...
                anisotropic_Raman_included,...
                prefactor,...
                N, num_modes);
a4 =       N_op(expDG.*(A_IP+a3*(sim.dz)),...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                haw_sponRS, hbw_sponRS,...
                A_sponRS, sponRS_prefactor,...
                Ra_sponRS, Rb_sponRS,...
                anisotropic_Raman_included,...
                prefactor,...
                N, num_modes);

A1 = expDG.*(A_IP + (a1+2*a2+2*a3)*(sim.dz/6)) + a4*(sim.dz/6);

% Local error estimate
a5 =       N_op(        A1,...
                sim,...
                SK_info, SRa_info, SRb_info,...
                Kerr, Ra, Rb,...
                haw, hbw,...
                haw_sponRS, hbw_sponRS,...
                A_sponRS, sponRS_prefactor,...
                Ra_sponRS, Rb_sponRS,...
                anisotropic_Raman_included,...
                prefactor,...
                N, num_modes);
err = sum(abs((a4-a5)*(sim.dz/10)).^2,1);

% Stepsize control
normA = sum(abs(A1).^2,1);
err = sqrt(err./normA);
err = max(err(normA~=0));
if isnan(err) % the computation is just so wrong, so we reduce the step size and do it again
    opt_dz = 0.5*sim.dz;
    success = false;
else
    opt_dz = max(0.5,min(2,0.8*(sim.adaptive_dz.threshold/err)^(1/4)))*sim.dz;

    success = err < sim.adaptive_dz.threshold;
end

end

function dAdz = N_op(A_w,...
                     sim,...
                     SK_info, SRa_info, SRb_info,...
                     Kerr, Ra, Rb,...
                     haw, hbw,...
                     haw_sponRS, hbw_sponRS,...
                     A_sponRS, sponRS_prefactor,...
                     Ra_sponRS, Rb_sponRS,...
                     anisotropic_Raman_included,...
                     prefactor,...
                     N, num_modes)
%N_op Calculate dAdz

A_t = fft(A_w);

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
                     int32(N), 1,...
                     num_modes,...
                     sim.cuda_num_operations_SRSK);
        if sim.include_sponRS % spontaneous Raman scattering
            Ra_sponRS = feval(sim.cuda_sponRS,...
                              Ra_sponRS,...
                              complex(A_t), A_sponRS,...
                              SRa_info.SRa,...
                              SRa_info.nonzero_midx1234s,...
                              SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                              int32(N), 1,...
                              num_modes);
        end
    else % polarized fields
        [Kerr,...
         Ra, Rb] = feval(sim.cuda_SRSK,...
                         Kerr, Ra, Rb,...
                         complex(A_t),...
                         SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                         SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                         SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                         sim.Raman_model~=0, sim.Raman_model==2,...
                         int32(N), 1,...
                         num_modes,...
                         sim.cuda_num_operations_SRSK);
        if sim.include_sponRS % spontaneous Raman scattering
            [Ra_sponRS, Rb_sponRS] = feval(sim.cuda_sponRS,...
                                           Ra_sponRS, Rb_sponRS,...
                                           complex(A_t), A_sponRS,...
                                           SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                           SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                                           sim.Raman_model==2,...
                                           int32(N), 1,...
                                           num_modes,...
                                           sim.cuda_num_operations_sponRS);
        end
    end
    Kerr = sum(Kerr,3);
else
    % If using the CPU, first precompute Ra_mn and Rb_mn.
    if sim.Raman_model ~= 0
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
        Ra_mn = A_t(:, SRa_info.nonzero_midx34s(1,:)).*conj(A_t(:, SRa_info.nonzero_midx34s(2,:))); % (N,num_nonzero34)
        if anisotropic_Raman_included
            SRb_nonzero_midx34s = midx34s_sub2ind(SRb_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
            Rb_mn = A_t(:, SRb_info.nonzero_midx34s(1,:)).*conj(A_t(:, SRb_info.nonzero_midx34s(2,:))); % (N,num_nonzero34)
        end
        if sim.include_sponRS % spontaneous Raman scattering
            Ra_mn_sponRS =      A_t(:, SRa_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, SRa_info.nonzero_midx34s(2,:))) +... % (N,num_nonzero34)
                           A_sponRS(:, SRa_info.nonzero_midx34s(1,:)).*conj(     A_t(:, SRa_info.nonzero_midx34s(2,:))) +...
                           A_sponRS(:, SRa_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, SRa_info.nonzero_midx34s(2,:)));
            if anisotropic_Raman_included
                Rb_mn_sponRS =      A_t(:, SRb_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, SRb_info.nonzero_midx34s(2,:))) +... % (N,num_nonzero34)
                               A_sponRS(:, SRb_info.nonzero_midx34s(1,:)).*conj(     A_t(:, SRb_info.nonzero_midx34s(2,:))) +...
                               A_sponRS(:, SRb_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, SRb_info.nonzero_midx34s(2,:)));
            end
        end
    end
    
    % Then calculate Kerr, Ra, and Rb.
    for midx1 = 1:num_modes
        % Kerr
        nz_midx1 = find( SK_info.nonzero_midx1234s(1,:)==midx1 );
        midx2 = SK_info.nonzero_midx1234s(2,nz_midx1);
        midx3 = SK_info.nonzero_midx1234s(3,nz_midx1);
        midx4 = SK_info.nonzero_midx1234s(4,nz_midx1);
        Kerr(:,midx1) = sum(permute(SK_info.SK(nz_midx1),[2 1]).*A_t(:, midx2).*A_t(:, midx3).*conj(A_t(:, midx4)),2);
        if sim.Raman_model ~= 0
            % Ra
            for midx2 = 1:num_modes
                nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 2nd-dimensional "num_nonzero34" of Ra_mn
                Ra(:, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[2 1]).*Ra_mn(:, idx),2);
                if sim.include_sponRS % spontaneous Raman scattering
                    Ra_sponRS(:, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[2 1]).*Ra_mn_sponRS(:, idx),2);
                end
            end
            % Rb
            if anisotropic_Raman_included
                for midx2 = 1:num_modes
                    nz_midx1 = find( SRb_info.nonzero_midx1234s(1,:)==midx1 );
                    nz_midx = nz_midx1( SRb_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                    midx3 = SRb_info.nonzero_midx1234s(3,nz_midx);
                    midx4 = SRb_info.nonzero_midx1234s(4,nz_midx);
                    idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                    idx = arrayfun(@(i) find(SRb_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of Rb_mn
                    Rb(:, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[2 1]).*Rb_mn(:, idx),2);
                    if sim.include_sponRS % spontaneous Raman scattering
                        Rb_sponRS(:, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[2 1]).*Rb_mn_sponRS(:, idx),2);
                    end
                end
            end
        end
    end
    if anisotropic_Raman_included
        clear Ra_mn Rb_mn;
        if sim.include_sponRS % spontaneous Raman scattering
            clear Ra_mn_sponRS Rb_mn_sponRS;
        end
    elseif sim.Raman_model~=0
        clear Ra_mn;
        if sim.include_sponRS % spontaneous Raman scattering
            clear Ra_mn_sponRS;
        end
    end
end

% Calculate h*Ra as F-1(h F(Ra))
% The convolution using Fourier Transform is faster if both arrays are
% large. If one of the array is small, "conv" can be faster.
% Please refer to
% "https://blogs.mathworks.com/steve/2009/11/03/the-conv-function-and-implementation-tradeoffs/"
% for more information.
if sim.Raman_model ~= 0
    Ra = fft(haw.*ifft(Ra));

    if sim.include_sponRS % spontaneous Raman scattering
        Ra_sponRS = fft(haw_sponRS.*ifft(Ra_sponRS).*sponRS_prefactor{2});
    end
    
    if anisotropic_Raman_included % polarized fields with an anisotropic Raman
        Rb = fft(hbw.*ifft(Rb));
        
        if sim.include_sponRS % spontaneous Raman scattering
            Rb_sponRS = fft(hbw_sponRS.*ifft(Rb_sponRS).*sponRS_prefactor{2});
        end
    end
    
    if ~anisotropic_Raman_included
        if sim.include_sponRS % spontaneous Raman scattering
            nonlinear = Kerr + sum((Ra+Ra_sponRS).*permute(A_t,[1 3 2]),3);
        else
            nonlinear = Kerr + sum(Ra.*permute(A_t,[1 3 2]),3);
        end
    else % polarized fields with an anisotropic Raman
        if sim.include_sponRS % spontaneous Raman scattering
            nonlinear = Kerr + sum((Ra+Rb+Ra_sponRS+Rb_sponRS).*permute(A_t,[1 3 2]),3);
        else
            nonlinear = Kerr + sum((Ra+Rb).*permute(A_t,[1 3 2]),3);
        end
    end
else
    nonlinear = Kerr;
end

% Now everything has been summed into "nonlinear", so transform into the
% frequency domain for the prefactor
dAdz = prefactor.*ifft(nonlinear);

end

%% helper function
function total_energy = calc_total_energy(A,N,dt)
%CALC_TOTAL_ENERGY
% total_energy = sum(trapz(abs(fftshift(A,1)).^2),2)*N*dt/1e3; % nJ
%
% "A" is the spectral intensity of the pulse under the convention of
% "ifft", so it needs "fftshift" to reorder it correctly for "trapz".
%
% In practice, because "fftshift" is slow, I use "sum" instead.
% The difference between "trapz" and "sum" is half the sum of the first and
% last elements.

center_idx = ceil(size(A,1)/2);
A2 = abs(A).^2;

total_energy = ( sum(A2(:))-sum(sum(A2([center_idx,center_idx+1],:)))/2 )*N*dt/1e3; % nJ;

end