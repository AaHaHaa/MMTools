function [signal_fields_out,Power_out,...
          save_z,save_dz,...
          T_delay_out,...
          N] = SteppingCaller_adaptive_rategain(sim, gain_rate_eqn,...
                                                save_z0, save_points,...
                                                initial_condition,...
                                                prefactor,...
                                                omegas, D_op,...
                                                SK_info, SRa_info, SRb_info,...
                                                haw, hbw,...
                                                haw_sponRS, hbw_sponRS, sponRS_prefactor)
%STEPPINGCALLER_ADAPTIVE_RATEGAIN It attains the field after propagation 
%inside the gain medium solved by the rate equations with the adaptive-step
%RK4IP algorithm.
%   
% The computation of this code is based on
%   1. Lindberg et al., "Accurate modeling of high-repetition rate ultrashort pulse amplification in optical fibers", Scientific Reports (2016)
%   2. Chen et al., "Optimization of femtosecond Yb-doped fiber amplifiers for high-quality pulse compression", Opt. Experss (2012)
%   3. Gong et al., "Numerical modeling of transverse mode competition in strongly pumped multimode fiber lasers and amplifiers", Opt. Express (2007)
%
%   Please go to "gain_info.m" file to find the information about some input arguments.
%   The info of some other input arguments are inside "GMMNLSE_propagate.m"
%
%
% Implementation of adaptive-step algorithm:
%   The counterpumping power is found at the pulse input end of the fiber 
%   with a binary search. The solver does only pulse-forward propagation, 
%   so the counterpumping power increases during the propagation until it 
%   reaches the maximum value at the pulse-output end, which is the input 
%   end of  the counterpumping power. If this maximum value is higher than 
%   the user-specified counterpumping power, we need to reduce the
%   counterpumping power at the pulse input end, and vice versa.

%% Error check
if gain_rate_eqn.linear_oscillator
    error('GMMNLSE_adaptive_rategain:settingsError',...
          'Adaptive-step method doesn''t support "linear_oscillator".');
end
if gain_rate_eqn.include_ASE
    error('GMMNLSE_adaptive_rategain:propagationDirectionError',...
          'Adaptive scheme works only for cases without ASE.')
end
if gain_rate_eqn.reuse_data
    error('GMMNLSE_adaptive_rategain:propagationDirectionError',...
          'Adaptive scheme works only for cases without reusing the data (for an oscillator to converge faster).')
end

%%
Nt = size(initial_condition.fields,1);

%% Pump direction
if gain_rate_eqn.copump_power == 0
    if gain_rate_eqn.counterpump_power == 0
        gain_rate_eqn.pump_direction = 'co'; % use 'co' for zero pump power
    else
        gain_rate_eqn.pump_direction = 'counter';
    end
else
    if gain_rate_eqn.counterpump_power == 0
        gain_rate_eqn.pump_direction = 'co';
    else
        gain_rate_eqn.pump_direction = 'bi';
    end
end

%% Propagations
% If the pumping includes counterpumping, do backward propagation once
% without the electric field to first build up the inversion and obtain an
% initial guess of the counterpumping power at the pulse-input end.
% This value will be used as a starting point for the binary search.
if ismember(gain_rate_eqn.pump_direction,{'counter','bi'})
    % Load the folder
    if ispc
        sep_char = '\';
    else % unix
        sep_char = '/';
    end
    current_path = mfilename('fullpath');
    sep_pos = strfind(current_path,sep_char);
    current_folder = current_path(1:sep_pos(end));
    addpath([current_folder 'linear oscillator']);
    
    max_Power_pump_backward = gain_propagate('backward',...
                                             sim,gain_rate_eqn,...
                                             save_points,save_z0,...
                                             Nt,prefactor,omegas,D_op,...
                                             SK_info,SRa_info,SRb_info,haw,hbw,...
                                             haw_sponRS, hbw_sponRS, sponRS_prefactor,...
                                             initial_condition,...
                                             []);
else
    max_Power_pump_backward = 0;
end
guess_Power_pump_backward = max_Power_pump_backward;

% Start the pulse propagation:
% If it includes counterpumping, we need to find the counterpumping power
% at the pulse-input end that gives the incident counterpumping power (at  
% the pulse-output end).
Power_pump_backward = {0};
binary_L = 0; binary_R = 0;
iterations = 1; % the number of iterations under counterpumping and bi-pumping
% The following counterpump_ratio_L/R restrict the boundaries of the
% modified binary search within 0.5~2 times the predefined
% counterpump_power.
counterpump_ratio = Power_pump_backward{end}/gain_rate_eqn.counterpump_power;
if isnan(counterpump_ratio)
    counterpump_ratio = 0;
end
counterpump_ratio_L = counterpump_ratio;
counterpump_ratio_R = counterpump_ratio;
while (binary_L == 0 || binary_R == 0) || ...
      (counterpump_ratio_L < 0.5 || counterpump_ratio_R > 2)
    [signal_fields,Power_pump_forward,Power_pump_backward,...
     save_z,save_dz,...
     T_delay_out,...
     N]          = gain_propagate('forward',...
                                   sim,gain_rate_eqn,...
                                   save_points,save_z0,...
                                   Nt,prefactor,omegas,D_op,...
                                   SK_info,SRa_info,SRb_info,haw,hbw,...
                                   haw_sponRS, hbw_sponRS, sponRS_prefactor,...
                                   initial_condition,...
                                   guess_Power_pump_backward);

    if isequal(gain_rate_eqn.pump_direction,'co')
        break;
    end
    % We need to find two counterpump powers at the pulse-input end, which
    % give pulse-output-end counterpump powers smaller and larger than the
    % actual counterpump power. Two boundary values are required for the 
    % binary search afterwards.
    extra_ratio = 0.5; % added multiplication ratio for setting the boundaries.
    if Power_pump_backward{end} < gain_rate_eqn.counterpump_power
        binary_L = guess_Power_pump_backward;
        binary_L_counterpump_power = Power_pump_backward{end};
        
        counterpump_ratio_L = binary_L_counterpump_power/gain_rate_eqn.counterpump_power;
        
        guess_ratio = max((1+extra_ratio)/counterpump_ratio_L,1+extra_ratio);
        guess_Power_pump_backward = guess_Power_pump_backward*guess_ratio;
        
        if gain_rate_eqn.verbose
            fprintf('Gain rate equation, iteration %u: counterpump power (at seed output end) = %7.6g(W)\n',iterations,Power_pump_backward{end});
            iterations = iterations + 1;
        end
    elseif Power_pump_backward{end} > gain_rate_eqn.counterpump_power
        binary_R = guess_Power_pump_backward;
        binary_R_counterpump_power = Power_pump_backward{end};
        
        counterpump_ratio_R = binary_R_counterpump_power/gain_rate_eqn.counterpump_power;
        
        guess_ratio = min(1/counterpump_ratio_R/(1+extra_ratio),1/(1+extra_ratio));
        guess_Power_pump_backward = guess_Power_pump_backward*guess_ratio;
        
        if gain_rate_eqn.verbose
            fprintf('Gain rate equation, iteration %u: counterpump power (at seed output end) = %7.6g(W)\n',iterations,Power_pump_backward{end});
            iterations = iterations + 1;
        end
    end
end

% Use the modified binary search to find the counterpump power at the input 
% end.
% The middle point is determined based on the computed counterpumping power 
% at both ends of the interval.
% Each data point contains counterpumping power at the pulse-input and 
% output end. We aim to find the counterpumping power at the pulse-input 
% end such that it gives the user-specified counterpumping power at the 
% pulse-output end.
if ~isequal(gain_rate_eqn.pump_direction,'co')
    while abs(Power_pump_backward{end}-gain_rate_eqn.counterpump_power)/gain_rate_eqn.counterpump_power > gain_rate_eqn.tol && iterations <= gain_rate_eqn.max_iterations
        % I call this "modified" binary search because I don't use the
        % middle point as the next upper or lower boundary value; instead,
        % I use the linear interpolation to find it. It's just faster!
        % Traditional binary serarch uses
        %    binary_m = (binary_L + binary_R)/2;
        binary_m = binary_L*(binary_R_counterpump_power-gain_rate_eqn.counterpump_power)/(binary_R_counterpump_power-binary_L_counterpump_power) + binary_R*(gain_rate_eqn.counterpump_power-binary_L_counterpump_power)/(binary_R_counterpump_power-binary_L_counterpump_power);
        
        [signal_fields,Power_pump_forward,Power_pump_backward,...
         save_z,save_dz,...
         T_delay_out,...
         N]          = gain_propagate('forward',...
                                       sim,gain_rate_eqn,...
                                       save_points,save_z,...
                                       Nt,prefactor,omegas,D_op,...
                                       SK_info,SRa_info,SRb_info,haw,hbw,...
                                       haw_sponRS, hbw_sponRS, sponRS_prefactor,...
                                       initial_condition,...
                                       binary_m);
        
        if Power_pump_backward{end} > gain_rate_eqn.counterpump_power
            binary_R = binary_m;
            binary_R_counterpump_power = Power_pump_backward{end};
        else
            binary_L = binary_m;
            binary_L_counterpump_power = Power_pump_backward{end};
        end
        
        if gain_rate_eqn.verbose
            fprintf('Gain rate equation, iteration %u: counterpump power (at seed output end) = %7.6g(W)\n',iterations,Power_pump_backward{end});
        end
        iterations = iterations + 1;
    end
end

%% Output:
% Power_out is in frequency domain while field_out is transformed into time domain.

% Transform the current "signal_fields" and "Power..." into arrays.
if sim.gpu_yes
    [signal_fields,Power_pump_forward,Power_pump_backward] = mygather(signal_fields,Power_pump_forward,Power_pump_backward);
    [signal_fields_out,Power_pump_forward_out,Power_pump_backward_out] = deal(signal_fields,Power_pump_forward,Power_pump_backward);
else
    [signal_fields_out,Power_pump_forward_out,Power_pump_backward_out] = deal(signal_fields,Power_pump_forward,Power_pump_backward);
end
signal_fields_out      = cell2mat(signal_fields_out);
Power_pump_forward_out = cell2mat(Power_pump_forward_out);
Power_pump_backward_out = cell2mat(Power_pump_backward_out);

Power_out.pump = struct('forward', Power_pump_forward_out,...
                        'backward',Power_pump_backward_out);
signal_fields_out = fft(signal_fields_out);

end

%%
function varargout = gain_propagate(direction,...
                                    sim,gain_rate_eqn,...
                                    save_points,save_z,...
                                    Nt,prefactor,omegas,D_op,...
                                    SK_info,SRa_info,SRb_info,haw,hbw,...
                                    haw_sponRS,hbw_sponRS,sponRS_prefactor,...
                                    initial_condition,...
                                    input_Power_pump_backward)
%GAIN_PROPAGATE Runs the corresponding propagation method based on "direction".

if sim.gpu_yes
    dummy_var = zeros(size(initial_condition.fields),'gpuArray');
else
    dummy_var = zeros(size(initial_condition.fields));
end

dt = initial_condition.dt;
num_modes = size(initial_condition.fields,2);

save_dz = zeros(save_points,1);
T_delay_out = zeros(save_points,1);

% Pulse centering based on the moment of its intensity
if sim.pulse_centering
    % Center the pulse
    TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(initial_condition.fields).^2,[1,2])/sum(abs(initial_condition.fields).^2,[1,2]));
    % Because circshift is slow on GPU, I discard it.
    %last_result = ifft(circshift(initial_condition.fields,-tCenter));
    if TCenter ~= 0
        if TCenter > 0
            initial_condition.fields = [initial_condition.fields(1+TCenter:end,:);initial_condition.fields(1:TCenter,:)];
        elseif TCenter < 0
            initial_condition.fields = [initial_condition.fields(end+1+TCenter:end,:);initial_condition.fields(1:end+TCenter,:)];
        end

        if sim.gpu_yes
            TCenter = gather(TCenter);
        end
        T_delay = TCenter*dt;
    else
        T_delay = 0;
    end
else
    T_delay = 0;
end
T_delay_out(1) = T_delay;
initial_condition.fields = ifft(initial_condition.fields);

switch direction
    case 'forward'
        [signal_fields,...
         Power_pump_forward,Power_pump_backward] = initialization(sim,gain_rate_eqn,...
                                                                  Nt,num_modes,...
                                                                  save_points,...
                                                                  initial_condition,input_Power_pump_backward);
        last_Power_pump_forward = Power_pump_forward{1};
        last_Power_pump_backward = Power_pump_backward{1};
        last_signal_fields      = signal_fields     {1}; % = initial_condition.fields
    case 'backward'
        last_Power_pump_backward = gain_rate_eqn.counterpump_power;
end

% Initialize N to be exported, the ion density of the upper state
N = zeros([size(gain_rate_eqn.N_total) save_points,length(gain_rate_eqn.energy_levels)-1]);

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

% max dz
if ~isfield(sim.adaptive_dz,'max_dz')
    sim.adaptive_dz.max_dz = sim.save_period/10;
end

sim.dz = 1e-6; % m; start with a small value to avoid initial blowup
save_dz(1) = sim.dz;

% Then start the propagation
z = 0;
save_i = 2; % the 1st one is the initial field
a5 = []; % the temporary variable in the forward propagation
p5 = []; % the temporary variable in the backward propagation
last_N = cat(8, ones([size(gain_rate_eqn.N_total),1,1,1,1,1,1]).*gain_rate_eqn.N_total,...
               zeros([size(gain_rate_eqn.N_total),1,1,1,1,1,length(gain_rate_eqn.energy_levels)-2])); % initial guess for solving the population during propagation
sim.last_dz = 1; % randomly put a number, 1, for initialization
GMMNLSE_rategain_func = str2func(['stepping_',sim.step_method,'_rategain_adaptive']);
while z+eps(z) < save_z(end) % eps(z) here is necessary due to the numerical error
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('GMMNLSE_propagate:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end
    
    switch direction
        case 'forward'
            ever_fail = false;
            previous_signal_fields = last_signal_fields;
            previous_a5 = a5;

            success = false;
            while ~success
                if ever_fail
                    last_signal_fields = previous_signal_fields;
                    a5 = previous_a5;
                end

                [last_signal_fields, a5,...
                 last_Power_pump_forward,last_Power_pump_backward,...
                 opt_dz,success,...
                 last_N] = GMMNLSE_rategain_func(last_signal_fields,last_N,...
                                                 dt,...
                                                 sim,gain_rate_eqn,...
                                                 SK_info,SRa_info,SRb_info,...
                                                 haw,hbw,...
                                                 haw_sponRS,hbw_sponRS,sponRS_prefactor,...
                                                 prefactor,omegas,D_op,...
                                                 last_Power_pump_forward,last_Power_pump_backward,a5,...
                                                 dummy_var);

                if ~success
                    ever_fail = true;

                    sim.dz = opt_dz;
                end
            end
            sim.last_dz = sim.dz; % previous dz
            
            % Apply the damped frequency window
            last_signal_fields = last_signal_fields.*sim.damped_freq_window;
            
            % Check for any NaN elements
            if any(any(isnan(last_signal_fields))) %any(isnan(last_signal_fields),'all')
                error('SteppingCaller_adaptive_rategain:NaNError',...
                      'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
            end

            % Center the pulse
            if sim.pulse_centering
                last_signal_fields_in_time = fft(last_signal_fields);
                TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(last_signal_fields_in_time).^2,[1,2])/sum(abs(last_signal_fields_in_time).^2,[1,2]));
                if ~isnan(TCenter) && TCenter ~= 0 % all-zero fields; for calculating ASE power only
                    % Because circshift is slow on GPU, I discard it.
                    %last_signal_fields = ifft(circshift(last_signal_fields_in_time,-tCenter));
                    if ~isempty(a5) % RK4IP reuses a5 from the previous step
                        a5 = fft(a5);
                    end
                    if TCenter > 0
                        if ~isempty(a5) % RK4IP reuses a5 from the previous step
                            a5 = ifft([a5(1+TCenter:end,:);a5(1:TCenter,:)]);
                        end
                        last_signal_fields = ifft([last_signal_fields_in_time(1+TCenter:end,:);last_signal_fields_in_time(1:TCenter,:)]);
                    elseif TCenter < 0
                        if ~isempty(a5) % RK4IP reuses a5 from the previous step
                            a5 = ifft([a5(end+1+TCenter:end,:);a5(1:end+TCenter,:)]);
                        end
                        last_signal_fields = ifft([last_signal_fields_in_time(end+1+TCenter:end,:);last_signal_fields_in_time(1:end+TCenter,:)]);
                    end
                    if sim.gpu_yes
                        TCenter = gather(TCenter);
                    end
                    T_delay = T_delay + TCenter*initial_condition.dt;
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
            % To account for this effect, I limit the step size to be 10x smaller 
            % than the effective maximum beat length which is
            % 2*pi/(max(eff_betas)-min(eff_betas)).
            % eff_betas is from betas, propagation constants, throughout the time 
            % window but scaled according to the spectral intensity to prevent 
            % taking into account betas where there's no light at all or where 
            % there's some light starting to grow.
            eff_range_D = find_range_D(abs(last_signal_fields).^2,imag(D_op));
            min_beat_length = 2*pi/eff_range_D;
            if strcmp(sim.step_method,'MPA')
                dz_resolve_beat_length = min_beat_length/4*sim.MPA.M;
            else
                dz_resolve_beat_length = min_beat_length/4;
            end
            % Because I use the approximation, sqrt(1+x)=1+x/2 if x is small, in
            % calculating signal fields with MPA, the code will give error here if
            % this approximation is bad.
            if isequal(sim.step_method,'MPA')
                relative_N1 = max(gain_rate_eqn.N_total-last_N(:,:,:,:,:,:,:,1));
                switch gain_rate_eqn.gain_medium
                    case {'Yb','Er','Nd'}
                        emi10 = gain_rate_eqn.cross_sections(:,:,:,:,:,:,:,2);
                        GSA01 = gain_rate_eqn.cross_sections(:,:,:,:,:,:,:,1);
                    case 'Tm'
                        emi10 = gain_rate_eqn.cross_sections(:,:,:,:,:,:,:,end-1);
                        GSA01 = gain_rate_eqn.cross_sections(:,:,:,:,:,:,:,1);
                end
                relative_gain = relative_N1.*emi10 - (max(gain_rate_eqn.N_total(:))-relative_N1).*GSA01;
                relative_gain(relative_gain<0) = 0; % I think (?) it only neds to resolve the "gain" correctly

                tol_approximation = 1e-4; % I've found that 1e-3 is not enough
                approx_error = @(x)abs((sqrt(1+x)-(1+x/2))./sqrt(1+x));
                while approx_error( 2*(opt_dz/sim.MPA.M*1e6)*max(relative_gain(:)) ) > tol_approximation
                    opt_dz = opt_dz/2;
                end
            end
            sim.dz = min([opt_dz,save_z(end)-z,sim.adaptive_dz.max_dz,dz_resolve_beat_length]);

            % If it's time to save, get the result from the GPU if necessary,
            % transform to the time domain, and save it
            if z == sim.last_dz
                ready_save_N = true;
            end
            if z >= save_z(end)-eps(z) % save the last N
                sim.small_dz = sim.dz/sim.MPA.M; % dummy variable for the following solve_gain_rate_eqn() to run correctly
                if gain_rate_eqn.counterpump_power == 0 % copumping
                    [~,~,~,...
                     ~,last_N] = solve_gain_rate_eqn('forward',...
                                                     sim,gain_rate_eqn,...
                                                     last_N,...
                                                     last_signal_fields,dummy_var,...
                                                     last_Power_pump_forward,dummy_var,...
                                                     dummy_var,dummy_var,...
                                                     omegas,dt,...
                                                     false );
                else % bi-pumping or counterpumping
                    [~,~,~,~,...
                     ~,last_N] = solve_gain_rate_eqn_linear_oscillator(sim,gain_rate_eqn,...
                                                                       last_N,...
                                                                       last_signal_fields,dummy_var,...
                                                                       last_Power_pump_forward,last_Power_pump_backward,...
                                                                       dummy_var,dummy_var,...
                                                                       omegas,dt);
                end

                if sim.gpu_yes
                    N(:,:,end,:) = permute(gather(last_N),[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
                else
                    N(:,:,end,:) = permute(last_N,[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
                end
                
                N = N/gather(max(gain_rate_eqn.N_total(:)));
            end
            % Below saves N, the population of each energy level:
            % Current N is computed by the next propagating step, so 
            % there is "ready_save_N" is to delay the saving process 
            % by one propagating step with the "save_i-1" command.
            % Coincidentally, this command also helps save the first
            % N.
            if ready_save_N
                if sim.gpu_yes
                    N(:,:,save_i-1,:) = permute(gather(last_N),[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
                else
                    N(:,:,save_i-1,:) = permute(last_N,[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
                end
                ready_save_N = false; % finish saving, then reset it back to false
            end
            if z >= save_z(save_i)-eps(z)
                if sim.gpu_yes
                    save_z(save_i) = gather(z);
                    save_dz(save_i) = gather(sim.last_dz);
                    Power_pump_forward{save_i} = gather(last_Power_pump_forward);
                    Power_pump_backward{save_i} = gather(last_Power_pump_backward);
                    signal_fields{save_i} = gather(last_signal_fields);
                else
                    save_z(save_i) = z;
                    save_dz(save_i) = sim.last_dz;
                    Power_pump_forward{save_i} = last_Power_pump_forward;
                    Power_pump_backward{save_i} = last_Power_pump_backward;
                    signal_fields{save_i} = last_signal_fields;
                end
                ready_save_N = true;

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
            
            % Force to break if the counterpumping power is too high;
            % otherwise, it'll run extremely slowly due to the adaptive
            % step. In such a situation, we know our counterpumping power
            % at the pulse-input end is too large, so there's no point
            % finishing this propagation accurately. I'll guess the
            % counterpumping power at the pulse-output end with the
            % following methods.
            if last_Power_pump_backward > gain_rate_eqn.counterpump_power*1.5
                % Approximate the final counterpumping power by assuming that it grows exponentially to the end
                %
                % Ballpark estimate with only two points
                % This can be inaccurate if the coefficient of the exponent
                % is far from 1.
                % P=exp(a*z)+b
                %    b=P(0)-1
                %    a=log(P(z)-b)/z
                b = input_Power_pump_backward - 1;
                a = log(last_Power_pump_backward-b)/z;
                if save_i < 5
                    Power_pump_backward{end} = exp(a*save_z(end))+b;
                    
                    warning('Please increase the number of save points so that the code can find the counterpumping power at the input with MATLAB "fit()" faster and more accurate.');
                else
                    % If the number of saved points is large enough,
                    % they're used to approximate the final backward pump power with fitting.
                    if sim.gpu_yes
                        pump_backward = squeeze(cell2mat(mygather(Power_pump_backward(1:save_i-1))));
                        b = gather(b); a = gather(b);
                    else
                        pump_backward = squeeze(cell2mat(Power_pump_backward(1:save_i-1)));
                    end
                    
                    % Use MATLAB "fit()" to fit the backward pump power to
                    % c0*exp(a0*d)+b0.
                    exp_fun = @(c0,a0,b0,d) c0*exp(a0*d)+b0;
                    fitopt = fittype(exp_fun,'coefficients',{'c0','a0','b0'},'independent',{'d'},'dependent',{'power'});
                    fitcurve = fit(save_z(1:save_i-1),pump_backward,fitopt,'Robust','Bisquare','StartPoint',[1,abs(a),b],'Lower',[0,-Inf,-Inf],'Display','off');
                    %figure;plot(save_z(1:save_i-1),pump_backward);hold on;plot(save_z,fitcurve(save_z));hold off;drawnow;
                    Power_pump_backward{end} = fitcurve(save_z(end));
                end
                break;
            end
            
        case 'backward'
            ever_fail = false;
            previous_Power_pump_backward = last_Power_pump_backward;
            previous_p5 = p5;
            
            success = false;
            while ~success
                if ever_fail
                    last_Power_pump_backward = previous_Power_pump_backward;
                    p5 = previous_p5;
                end

                [last_Power_pump_backward, p5,...
                 last_N,...
                 opt_dz,success] = pumpStepping_RK4IP_rategain_adaptive(sim, gain_rate_eqn,...
                                                                        last_Power_pump_backward, p5,...
                                                                        last_N,...
                                                                        omegas, dt,...
                                                                        dummy_var);

                if ~success
                    ever_fail = true;

                    sim.dz = opt_dz;
                end
            end
            sim.last_dz = sim.dz; % previous dz
            
            % Update z
            z = z + sim.dz;
            sim.dz = min([opt_dz,save_z(end)-z,sim.adaptive_dz.max_dz]);

            save_i = save_i + 1;

            % Report current status in the progress bar's message field
            if sim.progress_bar
                if z >= progress_bar_z(progress_bar_i)
                    waitbar(gather(z/save_z(end)),h_progress_bar,sprintf('%s%6.1f%%',sim.progress_bar_name,z/save_z(end)*100));
                    progress_bar_i = find(z<progress_bar_z,1);
                end
            end
    end
end

% Output
switch direction
    case 'forward'
        varargout = {signal_fields,Power_pump_forward,Power_pump_backward,...
                     save_z,save_dz,...
                     T_delay_out,...
                     N};
    case 'backward'
        varargout = {last_Power_pump_backward};
end

end

%% Initialization
function [signal_fields,...
          Power_pump_forward,Power_pump_backward] = initialization(sim,gain_rate_eqn,...
                                                                   Nt,num_modes,...
                                                                   save_points,...
                                                                   initial_condition,input_Power_pump_backward)
%INITIALIZATION initializes "signal_fields" and "Power" based on
%"segment_idx/num_segment".
%
%   If segment_idx = 1, they need to include copump_power, initial_fields, and initial forward ASE.
%   If segment_idx = num_segment(the last segment), "Power" needs to include counterpump_power if it's nonzero.
%                                                   They also need to include initial backward ASE.
%   Otherwise, they're just a zero matrix and a structure with zero matrices.
%
%   The reason I use cell arrays instead of a matrix (N,num_modes,z_points) for signal_fields and Power:
%       It's faster!
%       e.g. "signal_fields(:,:,zi) = signal_fields_next" is very slow.

    function output = initialize_zeros(mat_size)
        output = cell(1,1,save_points);
        output(:) = {zeros(mat_size)};
    end

% Pump power
Power_pump_forward  = initialize_zeros(1);
Power_pump_backward = initialize_zeros(1);

% Put in the necessary information.
Power_pump_forward{1} = gain_rate_eqn.copump_power;
Power_pump_backward{1} = input_Power_pump_backward;

% -------------------------------------------------------------------------
% Signal field
% "cell2mat" doesn't support "gpuArray" in a cell array, which affects the process when getting the output matrix. 
signal_fields = initialize_zeros([Nt,num_modes]);
signal_fields{1} = initial_condition.fields;

% -------------------------------------------------------------------------
% GPU
if sim.gpu_yes
    [signal_fields,Power_pump_forward,Power_pump_backward] = mygpuArray(signal_fields,Power_pump_forward,Power_pump_backward);
end

end

%% MYGPUARRAY
function varargout = mygpuArray(varargin)
%MYGPUARRAY

varargout = cell(1,nargin);
for i = 1:nargin
    varargout{i} = cellfun(@gpuArray,varargin{i},'UniformOutput',false);
end

end

%% MYGATHER
function varargout = mygather(varargin)
%MYGATHER

varargout = cell(1,nargin);
for i = 1:nargin
    varargout{i} = cellfun(@gather,varargin{i},'UniformOutput',false);
end

end
%% EFF_RANGE_D
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

%% CLEANMEUP
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end