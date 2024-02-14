function A1 = stepping_MPA_MMGaussianGain_rmc(A0, dt,...
                                              sim, prefactor,...
                                              SRa_info, SRb_info, SK_info,...
                                              D, rmc_D,...
                                              haw, hbw,...
                                              haw_sponRS, hbw_sponRS, sponRS_prefactor,...
                                              G, saturation_intensity)
%STEPPING_MPA_MMGAUSSIANGAIN_RMC Take one step with MPA with a spatially 
%saturating gain model
%
% Input:
%    A0 - initial field in the frequency domain (N, num_modes); sqrt(W)
%    dt - time grid point spacing; ps
%
%    sim.deltaZ - step size; m
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
%    sim.adaptive_deltaZ.threshold - the accuracy used to determined whether to increase or decrease the step size
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
%    D.pos - dispersion term exp(Dz) (N, num_modes, M+1)
%    D.neg - dispersion term exp(-Dz) (N, num_modes, M+1)
%    haw - isotropic Raman response in the frequency domain
%    hbw - anisotropic Raman response in the frequency domain
%
%    sponRS_prefactor - prefactor for the spontaneous Raman scattering
%
%    dummy_a5_1 - unused variable
%
%    G - gain term (N, 1); m^-1
%    saturation_intensity - saturation intensity; J/m^2
%    fr - Raman fraction
%
% Output:
%    A1 - the field (in the frequency domain) after one step size (N, num_modes)
%    dummy_a5 - unused variable
%    opt_deltaZ - recommended step size
%    success - whether the current step size is sufficiently small for the required tolerance

anisotropic_Raman_included = ~sim.scalar & sim.Raman_model==2;

[N,num_modes] = size(A0);

% Explanation about the computation of the transfer matrix when considering
% polarization modes:
%   transfer_matrix and Bmn are in (M+1,num_spatial_modes,num_spatial_modes,2)
%   "2" represents two polarization modes.
%   I put the odd modes(e.g., x-polarization) in the 1st one, the even(e.g., y-polarization) in the 2nd one.
%   Now there are two sub-arrays, each related to different polarizations,
%   then I can do the multiplication for each polarization and matrix division later.
if sim.scalar
    num_spatial_modes = num_modes;
    polar = 1;
else
    num_spatial_modes = num_modes/2;
    polar = 2;
end

% Spontaneous Raman scattering
if sim.include_sponRS
    if sim.gpu_yes
        Ra_sponRS = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Rb_sponRS = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));

        A_sponRS = fft(sponRS_prefactor{1}.*sqrt(abs(randn(N,sim.MPA.M+1,num_modes,'gpuArray'))).*exp(1i*2*pi*rand(N,sim.MPA.M+1,num_modes,'gpuArray')));
    else
        Ra_sponRS = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes));
        Rb_sponRS = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes));

        A_sponRS = fft(sponRS_prefactor{1}.*sqrt(abs(randn(N,sim.MPA.M+1,num_modes))).*exp(1i*2*pi*rand(N,sim.MPA.M+1,num_modes)));
    end
else
    Ra_sponRS = [];
    Rb_sponRS = [];
    A_sponRS = [];
end

% Set initial values for psi
psi = repmat(A0, 1, 1, sim.MPA.M+1); % M copies of psi(w,z) = A(w,z), in the frequency domain!
current_psi_half_parallelization = psi(:,:,sim.MPA.M+1);

for n_it = 1:sim.MPA.n_tot_max
    % Calculate A(w,z) at all z:
    %    A_w = D.pos*psi
    A_w = add_D_and_rmc(true,...
                        D,rmc_D,...
                        psi,...
                        N,num_modes,sim.MPA.M+1,...
                        sim.rmc.cuda_mypagemtimes);
    % Calculate A(t,z) at all z
    A_t = permute(fft(A_w),[1 3 2]); % (N, M+1, num_modes)

    % Set up matrices for the following Kerr, Ra, Rb, transfer_matrix computations
    if sim.gpu_yes
        Kerr = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Ra = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Rb = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));

        transfer_matrix = complex(zeros(sim.MPA.M+1, num_spatial_modes, num_spatial_modes, polar, 'gpuArray'));
    else
        Kerr = complex(zeros(N, sim.MPA.M+1, num_modes));
        Ra = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes));
        Rb = complex(zeros(N, sim.MPA.M+1, num_modes, num_modes));

        transfer_matrix = complex(zeros(sim.MPA.M+1, num_spatial_modes, num_spatial_modes, polar));
    end

    % Calculate Bmn for all l and m, and M
    % This is faster to do in MATLAB and only scales like num_modes^2
    %    Bmn(:, midx3, midx4,p) = squeeze(dt/1e12*sum(A_t(:, :, midx3,p).*conj(A_t(:, :, midx4,p)), 1))/saturation_intensity;
    if sim.scalar
        Bmn = dt/1e12*squeeze(sum(A_t.*conj(permute(A_t,[1 2 4 3]))))/saturation_intensity;
    else
        A_t_for_Bmn = cat(5,A_t(:,:,1:2:num_modes-1),A_t(:,:,2:2:num_modes)); % separate the polarization modes
        Bmn = shiftdim(dt/1e12*sum(A_t_for_Bmn.*conj(permute(A_t_for_Bmn,[1 2 4 3 5])))/saturation_intensity, 1);
    end
    
    % Calculate large num_modes^4 Kerr, Ra, and Rb terms, as well as the 
    % transfer matrix for the multimode Gaussian (spatially saturating) gain.
    % If not using the GPU, we will precompute Ra_mn and Rb_mn before the num_modes^4 sum
    if sim.gpu_yes
        % If using the GPU, do the computation with fast CUDA code
        if sim.scalar % scalar fields
            [Kerr,...
             Ra,...
             transfer_matrix] = feval(sim.cuda_SRSK,...
                                      Kerr, Ra, transfer_matrix,...
                                      complex(A_t),...
                                      complex(Bmn),...
                                      SK_info.SK, SRa_info.SRa,...
                                      SRa_info.nonzero_midx1234s,...
                                      SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                      sim.Raman_model~=0,...
                                      int32(N), sim.MPA.M+1,...
                                      num_modes,...
                                      sim.cuda_num_operations_SRSK);
            if sim.include_sponRS % spontaneous Raman scattering
                Ra_sponRS = feval(sim.cuda_sponRS,...
                                  Ra_sponRS,...
                                  complex(A_t), A_sponRS,...
                                  SRa_info.SRa,...
                                  SRa_info.nonzero_midx1234s,...
                                  SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                  int32(N), sim.MPA.M+1,...
                                  num_modes);
            end
        else % polarized fields
            [Kerr,...
             Ra, Rb,...
             transfer_matrix] = feval(sim.cuda_SRSK,...
                                      Kerr, Ra, Rb, transfer_matrix,...
                                      complex(A_t),...
                                      complex(Bmn),...
                                      SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                                      SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                      SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                                      sim.Raman_model~=0, sim.Raman_model==2,...
                                      int32(N), sim.MPA.M+1,...
                                      num_modes,...
                                      sim.cuda_num_operations_SRSK);
            if sim.include_sponRS % spontaneous Raman scattering
                [Ra_sponRS, Rb_sponRS] = feval(sim.cuda_sponRS,...
                                               Ra_sponRS, Rb_sponRS,...
                                               complex(A_t), A_sponRS,...
                                               SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                               SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                                               sim.Raman_model == 2,...
                                               int32(N), sim.MPA.M+1,...
                                               num_modes,...
                                               sim.cuda_num_operations_sponRS);
            end
        end
        Kerr = sum(Kerr,4);
    else
        % If using the CPU, first precompute Ra_mn and Rb_mn.
        if sim.Raman_model ~= 0
            if n_it == 1
                midx34s_sub2ind = @(x)...
                    cellfun(@(xx)...
                        feval(@(sub) sub2ind(num_modes*ones(1,2),sub{:}), num2cell(xx)),... % this extra "feval" is to get "xx", which is of the size 2x1, into the input arguments of "sub2ind", so transforming "xx" into a 2x1 cell, each containing an integer, and using {:} expansion is necessary
                    mat2cell(x,2,ones(1,size(x,2)))); %#ok transform (2,num_nonzero34) midx34s into linear indices of a num_modes-by-num_modes matrix
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
            if sim.include_sponRS % spontaneous Raman scattering
                Ra_mn_sponRS =      A_t(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, :, SRa_info.nonzero_midx34s(2,:))) +... % (N,M+1,num_nonzero34)
                               A_sponRS(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(     A_t(:, :, SRa_info.nonzero_midx34s(2,:))) +...
                               A_sponRS(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, :, SRa_info.nonzero_midx34s(2,:)));
                if ~sim.scalar
                    Rb_mn_sponRS =      A_t(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, :, SRb_info.nonzero_midx34s(2,:))) +... % (N,M+1,num_nonzero34)
                                   A_sponRS(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(     A_t(:, :, SRb_info.nonzero_midx34s(2,:))) +...
                                   A_sponRS(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, :, SRb_info.nonzero_midx34s(2,:)));
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
            Kerr(:,:,midx1) = sum(permute(SK_info.SK(nz_midx1),[3 2 1]).*A_t(:, :, midx2).*A_t(:, :, midx3).*conj(A_t(:, :, midx4)),3);
            if sim.Raman_model ~= 0
                % Ra
                for midx2 = 1:num_modes
                    nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                    nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                    midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                    midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                    idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                    idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of Ra_mn
                    Ra(:, :, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[3 2 1]).*Ra_mn(:, :, idx),3);
                    if sim.include_sponRS % spontaneous Raman scattering
                        Ra_sponRS(:, :, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[3 2 1]).*Ra_mn_sponRS(:, :, idx),3);
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
                        Rb(:, :, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[3 2 1]).*Rb_mn(:, :, idx),3);
                        if sim.include_sponRS % spontaneous Raman scattering
                            Rb_sponRS(:, :, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[3 2 1]).*Rb_mn_sponRS(:, :, idx),3);
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
        
        % Transfer matrix:
        % Explanation for the computation below:
        %   transfer_matrix: T, Bmn: B
        %   S_plmn: (p,l) relates to the indices of T
        %           (m,n) relates to the indices of B
        %   I calculate first the linear indices of (p,l) and (m,n), use
        %   these two indices to do the multiplication: T=S*B
        if sim.scalar
            lidx_T = sub2ind(num_spatial_modes*ones(1,2),SRa_info.nonzero_midx1234s(1,:)',SRa_info.nonzero_midx1234s(2,:)');
            lidx_B = sub2ind(num_spatial_modes*ones(1,2),SRa_info.nonzero_midx1234s(3,:)',SRa_info.nonzero_midx1234s(4,:)');
            msize_T = num_spatial_modes^2;
        else
            % odd, even represents polarizations, which gives an extra
            % addition of "num_spatial_modes^2" for the linear indices of
            % the even modes.
            odd_or_even = @(x) rem(x+1,2);
            % Recovery back to spatial mode indices to calculate the
            % linear indices related to T and B.
            recover_spatial_midx = @(x) ceil(x/2);
            
            oddeven_T = double(odd_or_even(SRa_info.nonzero_midx1234s(1,:)'));
            oddeven_B = double(odd_or_even(SRa_info.nonzero_midx1234s(3,:)'));
            lidx = recover_spatial_midx(SRa_info.nonzero_midx1234s);
            
            lidx_T = sub2ind(num_spatial_modes*ones(1,2),lidx(1,:)',lidx(2,:)') + oddeven_T*num_spatial_modes^2;
            lidx_B = sub2ind(num_spatial_modes*ones(1,2),lidx(3,:)',lidx(4,:)') + oddeven_B*num_spatial_modes^2;
            
            msize_T = num_spatial_modes^2*2;
        end
        for midx = 1:msize_T
            idx_T = lidx_T==midx;
            idx_B = lidx_B(idx_T);
            transfer_matrix(:,midx) = sum(SRa_info.SRa(idx_T)'.*squeeze(Bmn(:,idx_B)),2);
        end
    end
    
    % Finish the transfer_matrix by adding the diagonals and reshaping
    transfer_matrix = permute(transfer_matrix, [2, 3, 1, 4]) + eye(num_spatial_modes); % now it goes (num_spatial_modes, num_spatial_modes, M+1)
    A_w = permute(A_w, [1, 3, 2]); % A_w goes (N, num_modes, M+1)

    % Invert and apply the transfer function for each mode
    gain_rhs = permute(G.*A_w,[2 1 3]); % (num_modes, N, M+1) matrix, each column is a single w value
    if ~sim.scalar
        gain_rhs = cat(4,gain_rhs(1:2:num_modes-1,:,:),gain_rhs(2:2:num_modes,:,:));
        % "recovery_idx" is used to put the separated polarization modes
        % back into the same dimension of array.
        recovery_idx = [(1:num_spatial_modes);(1:num_spatial_modes)+num_spatial_modes];
        recovery_idx = recovery_idx(:);
    end
    if sim.gpu_yes
        gain_term = permute(pagefun(@mldivide, transfer_matrix,gain_rhs),[2 1 3 4]);
    else % CPU
        gain_term = permute(multbslash(transfer_matrix,gain_rhs),[2 1 3 4]);
    end
    % Put the separated polarization modes back into the same dimension of array.
    if ~sim.scalar
        gain_term = cat(2,gain_term(:,:,:,1),gain_term(:,:,:,2));
        gain_term = gain_term(:,recovery_idx,:);
    end
    gain_term = permute(gain_term,[1,3,2]); % (N, M+1, num_modes)
    
    % Apply the convolution for the isotropic Raman:
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
                nonlinear = Kerr + sum((Ra+Ra_sponRS).*permute(A_t,[1 2 4 3]),4);
            else
                nonlinear = Kerr + sum(Ra.*permute(A_t,[1 2 4 3]),4);
            end
        else % polarized fields with an anisotropic Raman
            if sim.include_sponRS % spontaneous Raman scattering
                nonlinear = Kerr + sum((Ra+Rb+Ra_sponRS+Rb_sponRS).*permute(A_t,[1 2 4 3]),4);
            else
                nonlinear = Kerr + sum((Ra+Rb).*permute(A_t,[1 2 4 3]),4);
            end
        end
    else
        nonlinear = Kerr;
    end
    
    % Multiply the nonlinear factor
    nonlinear = prefactor.*ifft(nonlinear);
    
    % Incorporate deltaZ and D.neg term for the integration.
    nonlinear = permute(nonlinear + gain_term,[1,3,2]); % (N,num_modes,M+1)
    nonlinear = add_D_and_rmc(false,...
                              D,rmc_D,...
                              sim.small_deltaZ*nonlinear,...
                              N,num_modes,sim.MPA.M+1,...
                              sim.rmc.cuda_mypagemtimes); % (N, num_modes, M+1)

    % Save the previous psi at the right end, then compute the new psis
    last_psi_half_parallelization = current_psi_half_parallelization;
    % ---------------------------------------------------------------------
    % "psi" is calculated from trapezoidal quadrature, which results in
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
                                            N,sim.MPA.M+1,num_modes);
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
A1 = add_D_and_rmc(true,...
                   D,rmc_D,...
                   psi(:,:,end),...
                   N,num_modes,1,...
                   sim.rmc.cuda_mypagemtimes);

end