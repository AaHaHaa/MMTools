function [A_out,...
          save_z,save_dz,...
          T_delay_out] = SteppingCaller_adaptive(sim,...
                                                 G, saturation_parameter,...
                                                 save_z, save_points,...
                                                 initial_condition,...
                                                 n2_prefactor,...
                                                 D_op,...
                                                 SK_info, SRa_info, SRb_info,...
                                                 haw, hbw,...
                                                 At_noise)
%STEPPINGCALLER_ADAPTIVE It starts the pulse propagation.

Nt = size(initial_condition.fields,1);
num_modes = size(initial_condition.fields,2);
dt = initial_condition.dt;

save_dz = zeros(save_points,1);
T_delay_out = zeros(save_points,1);

% Pulse centering based on the moment of its intensity
if sim.pulse_centering
    % Center the pulse
    temporal_profile = abs(initial_condition.fields).^2;
    temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
    TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,[1,2])/sum(temporal_profile,[1,2]));
    % Because circshift is slow on GPU, I discard it.
    %last_result = ifft(circshift(initial_condition.fields,-tCenter),[],1);
    if ~isnan(TCenter) && TCenter ~= 0
        if TCenter > 0
            initial_condition.fields = [initial_condition.fields(1+TCenter:end,:);initial_condition.fields(1:TCenter,:)];
        elseif TCenter < 0
            initial_condition.fields = [initial_condition.fields(end+1+TCenter:end,:);initial_condition.fields(1:end+TCenter,:)];
        end

        if sim.gpu_yes
            TCenter = gather(TCenter);
        end
        T_delay = TCenter*initial_condition.dt;
    else
        T_delay = 0;
    end
else
    T_delay = 0;
end
T_delay_out(1) = T_delay;

A_out = zeros(Nt, num_modes, save_points);

% Start by saving the initial condition
if sim.gpu_yes
    A_out(:,:,1) = gather(initial_condition.fields);
else
    A_out(:,:,1) = initial_condition.fields;
end

last_A = ifft(initial_condition.fields,[],1);

% Create a progress bar first
if sim.progress_bar
    if ~isfield(sim,'progress_bar_name')
        sim.progress_bar_name = '';
    elseif ~ischar(sim.progress_bar_name)
        error('GMMNLSE_propagate:ProgressBarNameError',...
            '"sim.progress_bar_name" should be a string.');
    end
    h_progress_bar = waitbar(0,sprintf('%s   0.0%%',sim.progress_bar_name),...
        'Name',sprintf('Running GMMNLSE: %s...',sim.progress_bar_name),...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h_progress_bar,'canceling',0);

    % Create the cleanup object
    cleanupObj = onCleanup(@()cleanMeUp(h_progress_bar));

    % Use this to control the number of updated time for the progress bar below 1000 times.
    num_progress_updates = 1000;
    progress_bar_z = (1:num_progress_updates)*save_z(end)/num_progress_updates;
    progress_bar_i = 1;
end

z = 0;
save_i = 2; % the 1st one is the initial field
a5 = [];
if ~isfield(sim.adaptive_dz,'max_dz')
    sim.adaptive_dz.max_dz = sim.save_period/10;
end
sim.dz = min(1e-6,sim.adaptive_dz.max_dz); % m; start with a small value to avoid initial blowup
save_dz(1) = sim.dz;
sim.last_dz = 1; % randomly put a number, 1, for initialization

switch sim.gain_model
    case 0
        gain_str = 'nogain';
    case 1
        if isscalar(sim.midx) % single-mode
            gain_str = 'SMGaussianGain'; % with RK4IP
        else % multimode
            gain_str = 'MMGaussianGain'; % with MPA
        end
end
GMMNLSE_func = str2func(['stepping_' sim.step_method '_', gain_str, '_adaptive']);

% Then start the propagation
while z+eps(z) < save_z(end) % eps(z) here is necessary due to the numerical error
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('SteppingCaller_adaptive:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end

    ever_fail = false;
    previous_A = last_A;
    previous_a5 = a5;

    success = false;
    while ~success
        if ever_fail
            last_A = previous_A;
            a5 = previous_a5;
        end

        [last_A, a5,...
         opt_dz, success] = GMMNLSE_func(last_A, dt,...
                                         sim, n2_prefactor,...
                                         SK_info, SRa_info, SRb_info,...
                                         D_op,...
                                         haw, hbw,...
                                         At_noise,...
                                         a5,...
                                         G, saturation_parameter);

        if ~success
            if opt_dz < 1e-10
                error('SteppingCaller_adaptive:adaptiveRK4IPError',...
                      'Adaptive RK4IP continues to fail.\nCheck simulation parameters.');
            end

            ever_fail = true;

            sim.dz = opt_dz;
        end
    end
    sim.last_dz = sim.dz; % previous dz
    
    % Check for any NaN elements
    if any(any(isnan(last_A))) %any(isnan(last_A),'all')
        error('SteppingCaller_adaptive:NaNError',...
              'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
    end

    % Center the pulse in the time window
    % Important:
    % In the modified shot-noise approach, the noise cannot be changed, so it needs to be translated as well.
    % This took me more than 12 hours of debugging to realize it.
    % Otherwise, the output field, if it has a strong frequency shift and shifts a lot in time relative to the time window, 
    % the noise without translation cannot match with the translated field,
    % resulting in a different noise field overlapped with the coherent pulse.
    % This will artificially create a noisy output.
    if sim.pulse_centering
        last_A_in_time = fft(last_A,[],1);
        temporal_profile = abs(last_A_in_time).^2;
        temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
        TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,[1,2])/sum(temporal_profile,[1,2]));
        % Because circshift is slow on GPU, I discard it.
        %last_result = ifft(circshift(last_A_in_time,-tCenter),[],1);
        if ~isnan(TCenter) && TCenter ~= 0
            if ~isempty(a5) % RK4IP reuses a5 from the previous step
                a5 = fft(a5,[],1);
            end
            if TCenter > 0
                if ~isempty(a5) % RK4IP reuses a5 from the previous step
                    a5 = ifft([a5(1+TCenter:end,:,:,:);a5(1:TCenter,:,:,:)],[],1);
                end
                last_A = ifft([last_A_in_time(1+TCenter:end,:);last_A_in_time(1:TCenter,:)],[],1);
                At_noise = cat(1,At_noise(1+TCenter:end,:,:),At_noise(1:TCenter,:,:));
            elseif TCenter < 0
                if ~isempty(a5) % RK4IP reuses a5 from the previous step
                    a5 = ifft([a5(end+1+TCenter:end,:,:,:);a5(1:end+TCenter,:,:,:)],[],1);
                end
                last_A = ifft([last_A_in_time(end+1+TCenter:end,:);last_A_in_time(1:end+TCenter,:)],[],1);
                At_noise = cat(1,At_noise(end+1+TCenter:end,:,:),At_noise(1:end+TCenter,:,:));
            end
            if sim.gpu_yes
                TCenter = gather(TCenter);
            end
            T_delay = T_delay + TCenter*dt;
        end
    end

    % Update z
    z = z + sim.dz;
    % Because the adaptive-step algorithm determines the step size by 
    % checking the error of the spectral intensities from RK3 and RK4, it
    % might ignore changes at the weakest part of the spectrum. This
    % happens in cases of noise-seeded four-wave mixing and noise-seeded
    % Raman scattering far from the pulse central frequency.
    %
    % To account for this effect, I limit the step size to be 4x smaller 
    % than the effective maximum beat length which is
    % 2*pi/(max(eff_betas)-min(eff_betas)).
    % eff_betas is from betas, propagation constants, throughout the time 
    % window but scaled according to the spectral intensity to prevent 
    % taking into account betas where there's no light at all or where 
    % there's some light starting to grow.
    eff_range_D = find_range_D(abs(last_A).^2,imag(D_op));
    min_beat_length = 2*pi/eff_range_D;
    if strcmp(sim.step_method,'MPA')
        dz_resolve_beat_length = min_beat_length/4*sim.MPA.M;
    else
        dz_resolve_beat_length = min_beat_length/4;
    end

    sim.dz = min([opt_dz,save_z(end)-z,sim.adaptive_dz.max_dz,dz_resolve_beat_length]);

    % If it's time to save, get the result from the GPU if necessary,
    % transform to the time domain, and save it
    if z >= save_z(save_i)-eps(z)
        A_out_ii = fft(last_A,[],1);
        if sim.gpu_yes
            save_dz(save_i) = gather(sim.last_dz);
            save_z(save_i) = gather(z);
            A_out(:, :, save_i) = gather(A_out_ii);
        else
            save_dz(save_i) = sim.last_dz;
            save_z(save_i) = z;
            A_out(:, :, save_i) = A_out_ii;
        end

        T_delay_out(save_i) = T_delay;

        save_i = save_i + 1;
    end

    % Report current status in the progress bar's message field
    if sim.progress_bar
        if z >= progress_bar_z(progress_bar_i)
            waitbar(gather(z/save_z(end)),h_progress_bar,sprintf('%s%6.1f%%',sim.progress_bar_name,z/save_z(end)*100));
            progress_bar_i = find(z<progress_bar_z,1);
        end
    end
end

end

%% Helper functions
function eff_range_D = find_range_D(spectrum,D)
%FIND_RANGE_D
%
% For an adaptive-dz method, the maximum dz is also limited by the 
% range of the propagation constant, beta0.
% If the FWM, Raman, or anything else happens for multiple frequencies 
% where dz can't resolve their beta0 difference, the outcome can be 
% wrong.
% Here, I multiply the (intensity)^(1/2) of the spectrum to the beta0 to 
% consider beta0 difference of the pulse and exclude those without the 
% pulse but within the frequency window. (1/2) is to maximize the 
% contribution of the weak part of the spectrum.

nonzero_field = max(spectrum)~=0;
spectrum = spectrum./max(spectrum(:));

eff_D = D(:,nonzero_field).*(spectrum(:,nonzero_field)).^(1/2); % I use ^(1/2) to emphasize the weak part
eff_range_D = max(eff_D(:)) - min(eff_D(:));

end
% =========================================================================
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end