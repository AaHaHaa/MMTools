function [signal_fields_out,T_delay_out,...
          Power_out,N,...
          saved_data] = SteppingCaller_rategain(sim,gain_rate_eqn,...
                                                num_zPoints,save_points,num_zPoints_persave,...
                                                initial_condition,...
                                                n2_prefactor,...
                                                SK_info, SRa_info, SRb_info,...
                                                omegas, D_op,...
                                                haw, hbw,...
                                                sponRS_prefactor,...
                                                At_noise,...
                                                saved_data)
%STEPPINGCALLER_RATEGAIN It attains the field after propagation inside the  
%gain medium solved by the rate equations.
%   
% The computation of this code is based on
%   1. Lindberg et al., "Accurate modeling of high-repetition rate ultrashort pulse amplification in optical fibers", Scientific Reports (2016)
%   2. Chen et al., "Optimization of femtosecond Yb-doped fiber amplifiers for high-quality pulse compression", Opt. Experss (2012)
%   3. Gong et al., "Numerical modeling of transverse mode competition in strongly pumped multimode fiber lasers and amplifiers", Opt. Express (2007)
%
%   Please go to "gain_info.m" file to find the information about some input arguments.
%   The information of some other input arguments are inside "GMMNLSE_propagate.m"
%

Nt = size(initial_condition.fields,1);
num_modes = size(initial_condition.fields,2);

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

%% Segment the entire propagation
% In situations where iterations are needed, such as counterpumping, all
% the information is saved during pulse propagation.
% However, if GPU is used, this might overload its memory. In such a
% situation, the propagation is segmented into several pieces. Only one
% piece is computed by GPU at each time instance while the rest is kept
% into RAM.

if ~sim.gpu_yes % use CPU
    num_segments = 1;
    segments = num_zPoints; % the number of z points of all segments; [num_segment1,num_segment2...]
    zPoints_each_segment = num_zPoints;
else % use GPU
    % Below calculates how to segment the data based on the predicted size (memory usage) of the total data.
    precision = 8; % "double" precision
    mem_complex_number = precision*2;
    % The size of the variables:
    variable_size.signal_fields = 2*Nt*num_modes*num_zPoints; % forward and backward
    variable_size.Power_pump = 2*num_zPoints; % forward and backward
    variable_size.Power_ASE  = 2*Nt*num_modes*num_zPoints; % forward and backward
    variable_size.signal_out_in_solve_gain_rate_eqn = Nt*num_modes^2*sim.MPA.M; % Because it sometimes blows up the GPU memory here, so I added it.
    variable_size.cross_sections = 2*Nt;
    variable_size.overlap_factor = numel(gain_rate_eqn.overlap_factor);
    variable_size.N_total = numel(gain_rate_eqn.N_total);
    variable_size.FmFnN = numel(gain_rate_eqn.FmFnN);
    variable_size.GammaN = numel(gain_rate_eqn.GammaN);
    var_field = fieldnames(variable_size);
    used_memory = 0;
    for i = 1:length(var_field)
        used_memory = used_memory + variable_size.(var_field{i});
    end
    used_memory = used_memory*mem_complex_number;

    num_segments = ceil(used_memory/gain_rate_eqn.memory_limit);
    if num_segments == 1
        segments = num_zPoints; % the number of z points of all segments; [num_segment1,num_segment2...]
    else
        zPoints_each_segment = ceil(num_zPoints/num_segments);
        num_segments = ceil(num_zPoints/zPoints_each_segment); % this operation is super important because a small increase of zPoints_each_segment (due to the ceiling function) can reduce num_segments by a few
        segments = [zPoints_each_segment*ones(1,num_segments-1) num_zPoints-zPoints_each_segment*(num_segments-1)]; % the number of z points of all segments; [num_segment1,num_segment2...]
        if segments(end) == 0
            segments = segments(1:end-1);
        end
    end
end

%% Initialization
% They're both in the frequency domain.
if isequal(gain_rate_eqn.pump_direction,'co')
    segment_idx = 1;
else % 'counter', 'bi'
    segment_idx = num_segments;
end

% Initialize everything in the specified segment, including both propagating forward and backward
% Everything is initialized because it's not sure which will be used.
% If forward propagation is done first, we need the backward data since in this code, both directions are, by default, considered.
% If backward propagation is done first, we need the forward data since in this code, both directions are, by default, considered.
[signal_fields,signal_fields_backward,...
 Power_pump_forward,Power_pump_backward,...
 Power_ASE_forward,Power_ASE_backward] = initialization('both',...
                                                        gain_rate_eqn,...
                                                        segment_idx,num_segments,...
                                                        Nt,num_modes,...
                                                        num_zPoints,...
                                                        initial_condition);

%% Propagations
% Load the saved data if "reuse_data=true" and they exist.
% This occurs when it's an oscillator where the data of each roundtrip is 
% saved for the next roundtrip to speed up the convergence.
%
% However, if this is the first time solving for the gain fiber in an 
% oscillator, there is no saved data and thus no need to load saved data.
need_to_load_saved_data = false;
first_iteration = true;
if gain_rate_eqn.reuse_data
    first_iteration = isempty(saved_data);
    if ~first_iteration % ready to load the saved data
        need_to_load_saved_data = true;
        
        Power_pump_forward     = saved_data.Power_pump_forward;
        Power_pump_backward    = saved_data.Power_pump_backward;
        Power_ASE_forward      = saved_data.Power_ASE_forward;
        Power_ASE_backward     = saved_data.Power_ASE_backward;
        signal_fields          = saved_data.signal_fields;
        signal_fields_backward = saved_data.signal_fields_backward;
    end
end

% =========================================================================
% Start the first pulse propagation
% =========================================================================
% If counter/bi-pumping, backward-propagate first without the signal to set up the population-inversion level.
if ismember(gain_rate_eqn.pump_direction,{'counter','bi'})
    need_to_load_saved_data = true;
    
    if first_iteration % when there is no saved data to load
                  gain_propagate('backward',...
                                 sim,gain_rate_eqn,...
                                 num_zPoints,segments,...
                                 need_to_load_saved_data,...
                                 save_points,num_zPoints_persave,...
                                 Nt,n2_prefactor,...
                                 SK_info,SRa_info,SRb_info,...
                                 omegas,D_op,...
                                 haw,hbw,...
                                 sponRS_prefactor,...
                                 At_noise,...
                                 initial_condition);
    end
end

% Start the forward pulse propagation
[T_delay_out,N] = gain_propagate('forward',...
                                 sim,gain_rate_eqn,...
                                 num_zPoints,segments,...
                                 need_to_load_saved_data,...
                                 save_points,num_zPoints_persave,...
                                 Nt,n2_prefactor,...
                                 SK_info,SRa_info,SRb_info,...
                                 omegas,D_op,...
                                 haw,hbw,...
                                 sponRS_prefactor,...
                                 At_noise,...
                                 initial_condition);

% -------------------------------------------------------------------------
% Initialize some variables for the following iterations
% -------------------------------------------------------------------------
if sim.gpu_yes
    pulse_energy = zeros(1,gain_rate_eqn.max_iterations,'gpuArray');
else
    pulse_energy = zeros(1,gain_rate_eqn.max_iterations);
end
if gain_rate_eqn.include_ASE
    if sim.gpu_yes
        energy_ASE_forward  = zeros(1,gain_rate_eqn.max_iterations,'gpuArray');
        energy_ASE_backward = zeros(1,gain_rate_eqn.max_iterations,'gpuArray');
    else
        energy_ASE_forward  = zeros(1,gain_rate_eqn.max_iterations);
        energy_ASE_backward = zeros(1,gain_rate_eqn.max_iterations);
    end
end

% =========================================================================
% Start the iterations (the 2nd, 3rd, ... pulse propagations) if necessary
% =========================================================================
% There's no need to compute backward propagation if copumping and ignoring ASE,
% so the previous first forward propagation solves everything.
if ~isequal(gain_rate_eqn.pump_direction,'co') || gain_rate_eqn.include_ASE
    need_to_load_saved_data = true;

    pulse_energy(1) = calc_total_energy(signal_fields{end},Nt,initial_condition.dt); % nJ
    if gain_rate_eqn.include_ASE
        energy_ASE_forward(1)  = sum(Power_ASE_forward{end}(:))/(Nt*initial_condition.dt); % W
        energy_ASE_backward(1) = sum(Power_ASE_backward{1}(:)) /(Nt*initial_condition.dt); % W
    end
    if gain_rate_eqn.verbose
        fprintf('Gain rate equation, iteration %u: pulse energy = %7.6g(nJ)\n',1,pulse_energy(1));
        if gain_rate_eqn.include_ASE
            fprintf('                                 forward  ASE power = %7.6g(mW)\n',energy_ASE_forward(1) *1e3);
            fprintf('                                 backward ASE power = %7.6g(mW)\n',energy_ASE_backward(1)*1e3);
        end
    end
    finish_iteration = false;
    for i = 2:gain_rate_eqn.max_iterations
        % -----------------------------------------------------------------
        % Propagate back and forth
        % -----------------------------------------------------------------
        % Backward propagation
                          gain_propagate('backward',...
                                         sim,gain_rate_eqn,...
                                         num_zPoints,segments,...
                                         need_to_load_saved_data,...
                                         save_points,num_zPoints_persave,...
                                         Nt,n2_prefactor,...
                                         SK_info,SRa_info,SRb_info,...
                                         omegas,D_op,...
                                         haw,hbw,...
                                         sponRS_prefactor,...
                                         At_noise,...
                                         initial_condition);

        % Forward propagation
        [T_delay_out,N] = gain_propagate('forward',...
                                         sim,gain_rate_eqn,...
                                         num_zPoints,segments,...
                                         need_to_load_saved_data,...
                                         save_points,num_zPoints_persave,...
                                         Nt,n2_prefactor,...
                                         SK_info,SRa_info,SRb_info,...
                                         omegas,D_op,...
                                         haw,hbw,...
                                         sponRS_prefactor,...
                                         At_noise,...
                                         initial_condition);

        % -----------------------------------------------------------------
        % Check convergence
        % -----------------------------------------------------------------
        pulse_energy(i) = calc_total_energy(signal_fields{end},Nt,initial_condition.dt); % nJ
        if gain_rate_eqn.include_ASE
            energy_ASE_forward(i)  = sum(Power_ASE_forward{end}(:))/(Nt*initial_condition.dt); % W
            energy_ASE_backward(i) = sum(Power_ASE_backward{1}(:)) /(Nt*initial_condition.dt); % W
        end
        if gain_rate_eqn.verbose
            fprintf('Gain rate equation, iteration %u: pulse energy = %7.6g(nJ)\n',i,pulse_energy(i));
            if gain_rate_eqn.include_ASE
                fprintf('                                 forward  ASE power (at seed output end) = %7.6g(mW)\n',energy_ASE_forward(i) *1e3);
                fprintf('                                 backward ASE power (at seed input end)  = %7.6g(mW)\n',energy_ASE_backward(i)*1e3);
            end
        end
        if pulse_energy(i-1)==0 || ... % no propagating pulse; this is to check the validity of the following pulse_energy convergence check due to the 1/pulse_energy(i-1) factor
           abs((pulse_energy(i)-pulse_energy(i-1))./pulse_energy(i-1)) < gain_rate_eqn.tol % pulse reaches a steady state
            if gain_rate_eqn.include_ASE % Both ASEs reach a steady state
                ASE_converged = ( energy_ASE_forward (i-1) ~= 0 && abs((energy_ASE_forward (i)-energy_ASE_forward (i-1))./energy_ASE_forward (i-1)) < gain_rate_eqn.tol ) && ...
                                ( energy_ASE_backward(i-1) ~= 0 && abs((energy_ASE_backward(i)-energy_ASE_backward(i-1))./energy_ASE_backward(i-1)) < gain_rate_eqn.tol );
            else % ignore ASE
                ASE_converged = true;
            end

            % Entering this if-statement means the pulse has converged.
            % If ASEs are converged as well, everything is converged and
            % the pulse propagation of this gain fiber is done.
            if ASE_converged
                finish_iteration = true;
            end
        end
        % Plot the convergence
        if gain_rate_eqn.verbose && i > 10 && ~(i==10+1 && finish_iteration) % this final condition means that the iteration finishes right when it's the 11th iteration. Don't plot anything in this situation.
            if i == 10+1 % plot it the first time; initialize the figure
                fig_gain_iterations = figure('Name','Gain iterations');
            end
            figure(fig_gain_iterations);
            if gain_rate_eqn.include_ASE
                yyaxis left;
            end
            plot(1:i,pulse_energy(1:i),'linewidth',2,'Color','b');
            set(gca,'YColor','b');
            xlabel('Iterations'); ylabel('Pulse energy (nJ)');
            if gain_rate_eqn.include_ASE
                yyaxis right;
                plot(1:i,energy_ASE_forward(1:i),'linewidth',2);
                ylabel('Forward ASE power (mW)');
            end
            set(gca,'fontsize',20);
            drawnow;
            % Close the convergence plot when it's done
            if finish_iteration
                close(fig_gain_iterations);
            end
        end
        
        % -----------------------------------------------------------------
        % Done!
        % -----------------------------------------------------------------
        if finish_iteration
            break;
        end
        % -----------------------------------------------------------------
        % Error! Iterations of the gain-fiber computation isn't finished
        % within "max_iterations" iterations.
        % -----------------------------------------------------------------
        if i == gain_rate_eqn.max_iterations
            warning('GMMNLSE_rategain:NotConvergedError',...
                    ['The iteration of forward and backward propagation of the gain fiber doesn''t converge to a steady state within %u iterations.\n',...
                    'Please run it again with a larger number of maximum iterations.'],gain_rate_eqn.max_iterations);
            break;
        end
    end
end

%% Output:
saved_zPoints = 1:num_zPoints_persave:num_zPoints;

% Change the size back to (N,num_modes)
Power_ASE_forward_out  = cellfun(@(P) permute(P,[5 3 1 2 4]), Power_ASE_forward, 'UniformOutput',false);
Power_ASE_backward_out = cellfun(@(P) permute(P,[5 3 1 2 4]), Power_ASE_backward,'UniformOutput',false);

% Transform them into arrays
signal_fields_out          = fft(cell2mat(signal_fields(:,:,saved_zPoints)));
%signal_fields_backward_out = fft(cell2mat(signal_fields_backward(:,:,saved_zPoints)));
Power_pump_forward_out     = cell2mat(Power_pump_forward (:,:,saved_zPoints));
Power_pump_backward_out    = cell2mat(Power_pump_backward(:,:,saved_zPoints));
Power_ASE_forward_out      = fftshift(cell2mat(Power_ASE_forward_out (:,:,saved_zPoints)),1);
Power_ASE_backward_out     = fftshift(cell2mat(Power_ASE_backward_out(:,:,saved_zPoints)),1);

%signal_fields_out = struct('forward', signal_fields_out,...
%                           'backward',signal_fields_backward_out);
Power_out = struct('pump',struct('forward',Power_pump_forward_out,'backward', Power_pump_backward_out),...
                   'ASE' ,struct('forward', Power_ASE_forward_out,'backward', Power_ASE_backward_out));

if gain_rate_eqn.reuse_data
    % Reverse the order and save the data for the linear-oscillator scheme for the next round
    if gain_rate_eqn.linear_oscillator
        reverse_direction = @(x) flip(x,3);
        Power_pump_forward_reuse     = reverse_direction(Power_pump_backward);
        Power_pump_backward_reuse    = reverse_direction(Power_pump_forward);
        Power_ASE_forward_reuse      = reverse_direction(Power_ASE_backward);
        Power_ASE_backward_reuse     = reverse_direction(Power_ASE_forward);
        signal_fields_backward_reuse = reverse_direction(signal_fields);
        signal_fields_reuse          = signal_fields_backward_reuse; % dummy saved "signal_fields_reuse"
        [Power_pump_forward,Power_pump_backward,...
         Power_ASE_forward,Power_ASE_backward,...
         signal_fields,signal_fields_backward] = deal(Power_pump_forward_reuse,Power_pump_backward_reuse,...
                                                      Power_ASE_forward_reuse,Power_ASE_backward_reuse,...
                                                      signal_fields_reuse,signal_fields_backward_reuse);
    end
    saved_data.signal_fields          = signal_fields;
    saved_data.signal_fields_backward = signal_fields_backward;
    saved_data.Power_pump_forward     = Power_pump_forward;
    saved_data.Power_pump_backward    = Power_pump_backward;
    saved_data.Power_ASE_forward      = Power_ASE_forward;
    saved_data.Power_ASE_backward     = Power_ASE_backward;
end

%%
function [T_delay_out,N] = gain_propagate(direction,...
                                          sim,gain_rate_eqn,...
                                          num_zPoints,segments,...
                                          need_to_load_saved_data,...
                                          save_points,num_zPoints_persave,...
                                          Nt,n2_prefactor,...
                                          SK_info,SRa_info,SRb_info,...
                                          omegas,D_op,...
                                          haw,hbw,...
                                          sponRS_prefactor,...
                                          At_noise,...
                                          initial_condition)
%GAIN_PROPAGATE Runs the corresponding propagation method based on "direction".

dt = initial_condition.dt;

T_delay_out = zeros(save_points,1);

zPoints_each_segment = segments(1);

% The 1st backward-propagation for counter/bi-pumping which happens before iteration starts.
if isequal(direction,'backward') && ~need_to_load_saved_data
    first_backward_before_iterations = true;
else
    first_backward_before_iterations = false;
end

% Pulse centering based on the moment of its intensity
if sim.pulse_centering
    % Center the pulse
    TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(initial_condition.fields).^2,[1,2])/sum(abs(initial_condition.fields).^2,[1,2]));
    % Because circshift is slow on GPU, I discard it.
    %last_result = ifft(circshift(initial_condition.fields,-tCenter));
    if ~isnan(TCenter) && TCenter ~= 0
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

% Initialization
switch direction
    case 'forward'
        [signal_fields,...
         Power_pump_forward,...
         Power_ASE_forward] = initialization('forward',...
                                             gain_rate_eqn,...
                                             1,num_segments,...
                                             Nt,num_modes,...
                                             num_zPoints,...
                                             initial_condition);
        % Put the first segment into GPU
        if sim.gpu_yes
            [signal_fields,signal_fields_backward,...
             Power_pump_forward,Power_pump_backward,...
             Power_ASE_forward,Power_ASE_backward] = mygpuArray2(1,segments(1),...
                                                                 signal_fields,signal_fields_backward,...
                                                                 Power_pump_forward,Power_pump_backward,...
                                                                 Power_ASE_forward,Power_ASE_backward);
        end
    case 'backward'
        % signal_fields is taken from the input argument
        [Power_pump_backward,...
         Power_ASE_backward] = initialization('backward',...
                                              gain_rate_eqn,...
                                              num_segments,num_segments,...
                                              Nt,num_modes,...
                                              num_zPoints,...
                                              initial_condition);
        % Put the last segment into GPU
        if sim.gpu_yes
            [signal_fields,signal_fields_backward,...
             Power_pump_forward,Power_pump_backward,...
             Power_ASE_forward,Power_ASE_backward] = mygpuArray2(num_zPoints-segments(end)+1,num_zPoints,...
                                                                 signal_fields,signal_fields_backward,...
                                                                 Power_pump_forward,Power_pump_backward,...
                                                                 Power_ASE_forward,Power_ASE_backward);
        end
end

% Initialize N to be exported, the ion density of the upper state
N = zeros([size(gain_rate_eqn.N_total),save_points,length(gain_rate_eqn.energy_levels)-1]);

if sim.progress_bar
    if ~isfield(sim,'progress_bar_name')
        sim.progress_bar_name = '';
    elseif ~ischar(sim.progress_bar_name)
        error('GMMNLSE_propagate:ProgressBarNameError',...
              '"sim.progress_bar_name" should be a string.');
    end
    h_progress_bar = waitbar(0,sprintf('%s   0.0%%',sim.progress_bar_name),...
        'Name',sprintf('Running GMMNLSE (%s): %s...',direction,sim.progress_bar_name),...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h_progress_bar,'canceling',0);
    
    % Create the cleanup object
    cleanupObj = onCleanup(@()cleanMeUp(h_progress_bar));
    
    % Use this to control the number of updated time for the progress bar below 1000 times.
    count_progress_bar = 1;
    num_progress_updates = 1000;
end

% Then start the propagation
last_N = zeros([size(gain_rate_eqn.N_total),1,1,1,1,1,length(gain_rate_eqn.energy_levels)-1]); % initial guess for solving the population during propagation
GMMNLSE_rategain_func = str2func(['stepping_',sim.step_method,'_rategain']);
for ii = 2:num_zPoints
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('GMMNLSE_propagate:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end
    
    % =====================================================================
    % GMMNLSE: Run the correct step function depending on the options chosen.
    % =====================================================================
    switch direction
        % -----------------------------------------------------------------
        % Forward propagation
        % -----------------------------------------------------------------
        case 'forward'
            Zi = ii; % z-index of the forward-propagating distance
            
            % Load the initial powers and field
            if Zi == 2 % the first/starting z-index
                last_Power_pump_forward = Power_pump_forward{1};
                last_Power_ASE_forward  = Power_ASE_forward {1};
                last_signal_fields      = signal_fields     {1}; % = initial_condition.fields
            end
            last_signal_fields_backward = signal_fields_backward{Zi-1};
            last_Power_pump_backward    = Power_pump_backward   {Zi-1};
            last_Power_ASE_backward     = Power_ASE_backward    {Zi-1};
            
            [last_signal_fields,...
             last_Power_pump_forward,last_Power_ASE_forward,...
             last_N] = GMMNLSE_rategain_func(last_signal_fields,last_signal_fields_backward,...
                                             last_Power_pump_forward,last_Power_pump_backward,...
                                             last_Power_ASE_forward,last_Power_ASE_backward,...
                                             last_N,...
                                             dt,sim,n2_prefactor,...
                                             SK_info,SRa_info,SRb_info,...
                                             omegas,D_op,...
                                             haw,hbw,...
                                             sponRS_prefactor,...
                                             At_noise,...
                                             gain_rate_eqn,...
                                             first_backward_before_iterations);
                                          
            % Apply the damped frequency window
            last_signal_fields = last_signal_fields.*sim.damped_freq_window;
            
            % Update "forward" only
            Power_pump_forward{Zi} = last_Power_pump_forward;
            Power_ASE_forward {Zi} = last_Power_ASE_forward;
            signal_fields     {Zi} = last_signal_fields;
            
            % Save N
            if Zi == num_zPoints % save the last N
                [~,~,~,...
                 ~,last_N] = solve_gain_rate_eqn('forward',...
                                                 sim,gain_rate_eqn,...
                                                 last_N,...
                                                 last_signal_fields,last_signal_fields_backward,...
                                                 last_Power_pump_forward,last_Power_pump_backward,...
                                                 last_Power_ASE_forward,last_Power_ASE_backward,...
                                                 omegas,dt,...
                                                 first_backward_before_iterations );
                         
                if sim.gpu_yes
                    N(:,:,end,:) = permute(gather(last_N),[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
                else
                    N(:,:,end,:) = permute(last_N,[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
                end
                
                N = N/gather(max(gain_rate_eqn.N_total(:)));
            end
            % Current N is computed by the next propagating step, so there is Zi-2, not Zi-1.
            % Coincidentally, this command also helps save the first N.
            if rem(Zi-2, num_zPoints_persave) == 0
                if sim.gpu_yes
                    N(:,:,int64((Zi-2)/num_zPoints_persave+1),:) = permute(gather(last_N),[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
                else
                    N(:,:,int64((Zi-2)/num_zPoints_persave+1),:) = permute(last_N,[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
                end
            end
            
            % Whenever at the end of each segment, gather the current data
            % from GPU back to the RAM and then throw those in the next
            % segment to GPU.
            if rem(Zi,zPoints_each_segment) == 0 || Zi == num_zPoints
                ready_for_next_segment_into_GPU = true;
            else
                ready_for_next_segment_into_GPU = false;
            end
        % -----------------------------------------------------------------
        % Backward propagation
        % -----------------------------------------------------------------
        case 'backward'
            Zi = num_zPoints+1 - ii; % z-index of the backward-propagating distance
            
            % Load the final powers
            if Zi == num_zPoints+1 - 2 % ii=2; the first/starting z-index
                last_Power_pump_backward = Power_pump_backward{num_zPoints};
                last_Power_ASE_backward  = Power_ASE_backward {num_zPoints};
            end
            last_Power_pump_forward     = Power_pump_forward {Zi+1};
            last_Power_ASE_forward      = Power_ASE_forward  {Zi+1};
            
            last_signal_fields          = signal_fields         {Zi+1};
            last_signal_fields_backward = signal_fields_backward{Zi+1};
            
            [last_Power_pump_backward,last_Power_ASE_backward,...
             last_N] = solve_gain_rate_eqn('backward',...
                                           sim,gain_rate_eqn,...
                                           last_N,...
                                           last_signal_fields,last_signal_fields_backward,...
                                           last_Power_pump_forward,last_Power_pump_backward,...
                                           last_Power_ASE_forward,last_Power_ASE_backward,...
                                           omegas,dt,...
                                           first_backward_before_iterations);
            
            % Update "backward" only
            Power_pump_backward{Zi} = last_Power_pump_backward;
            Power_ASE_backward {Zi} = last_Power_ASE_backward;
            
            % Whenever at the beginning of each segment, gather the current
            % data from GPU back to the RAM and then throw those in the
            % next segment to GPU.
            if rem(Zi,zPoints_each_segment) == 1
                ready_for_next_segment_into_GPU = true;
            else
                ready_for_next_segment_into_GPU = false;
            end
    end
    
    % =====================================================================
    % Some post-stepping checks, saves, and updates
    % =====================================================================
    if isequal(direction,'forward')
        % Check for any NaN elements
        if any(any(isnan(last_signal_fields))) %any(isnan(last_signal_fields),'all')
            error('SteppingCaller_rategain:NaNError',...
                  'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
        end
        
        % Center the pulse
        if sim.pulse_centering
            last_signal_fields_in_time = fft(last_signal_fields);
            TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(last_signal_fields_in_time).^2,[1,2])/sum(abs(last_signal_fields_in_time).^2,[1,2]));
            if ~isnan(TCenter) && TCenter ~= 0 % all-zero fields; for calculating ASE power only
                % Because circshift is slow on GPU, I discard it.
                %last_signal_fields = ifft(circshift(last_signal_fields_in_time,-tCenter));
                if TCenter > 0
                    last_signal_fields = ifft([last_signal_fields_in_time(1+TCenter:end,:);last_signal_fields_in_time(1:TCenter,:)]);
                elseif TCenter < 0
                    last_signal_fields = ifft([last_signal_fields_in_time(end+1+TCenter:end,:);last_signal_fields_in_time(1:end+TCenter,:)]);
                end
                if sim.gpu_yes
                    T_delay = T_delay + gather(TCenter*dt);
                else
                    T_delay = T_delay + TCenter*dt;
                end
            end
            if rem(Zi-1, num_zPoints_persave) == 0
                T_delay_out(int64((Zi-1)/num_zPoints_persave)+1) = T_delay;
            end
        end
    end
    
    % Gather the current GPU data back to RAM and those in the next segment from RAM to GPU
    if ready_for_next_segment_into_GPU
        current_segment_idx = ceil(Zi/zPoints_each_segment);
        switch direction
            case 'forward'
                next_segment_idx = current_segment_idx + 1;
            case 'backward'
                next_segment_idx = current_segment_idx - 1;
        end

        cumsum_segments = [0,cumsum(segments)];
        start_idx = cumsum_segments(current_segment_idx)+1;
        end_idx   = cumsum_segments(current_segment_idx+1);
        [signal_fields,signal_fields_backward,...
         Power_pump_forward,Power_pump_backward,...
         Power_ASE_forward,Power_ASE_backward] = mygather2(start_idx,end_idx,...
                                                           signal_fields,signal_fields_backward,...
                                                           Power_pump_forward,Power_pump_backward,...
                                                           Power_ASE_forward,Power_ASE_backward);

        if next_segment_idx > 0 && next_segment_idx <= length(segments)
            start_idx = cumsum_segments(next_segment_idx)+1;
            end_idx   = cumsum_segments(next_segment_idx+1);
            [signal_fields,signal_fields_backward,...
             Power_pump_forward,Power_pump_backward,...
             Power_ASE_forward,Power_ASE_backward] = mygpuArray2(start_idx,end_idx,...
                                                                 signal_fields,signal_fields_backward,...
                                                                 Power_pump_forward,Power_pump_backward,...
                                                                 Power_ASE_forward,Power_ASE_backward);
        end
    end
    
    % Report current status in the progress bar's message field
    if sim.progress_bar
        if num_zPoints < num_progress_updates || floor((ii-1)/((num_zPoints-1)/num_progress_updates)) == count_progress_bar
            waitbar((ii-1)/(num_zPoints-1),h_progress_bar,sprintf('%s%6.1f%%',sim.progress_bar_name,(ii-1)/(num_zPoints-1)*100));
            count_progress_bar = count_progress_bar+1;
        end
    end
end

end

end

%% initialization
function varargout = initialization(direction,...
                                    gain_rate_eqn,...
                                    segment_idx,num_segments,...
                                    Nt,num_modes,...
                                    zPoints,...
                                    initial_condition)
%INITIALIZATION initializes "signal_fields" and "Powers" based on
%"segment_idx/num_segment".
%
%   If segment_idx = 1, they need to include copump_power, initial_fields, and initial forward ASE.
%   If segment_idx = num_segment(the last segment), "Power" needs to include counterpump_power if it's nonzero.
%                                                   They also need to include initial backward ASE.
%   Otherwise, they're just structures with zero matrices.
%
%   The reason I use cell arrays instead of a matrix (N,num_modes,zPoints) for signal_fields and Power:
%       It's faster!
%       e.g. "signal_fields(:,:,zi) = signal_fields_next" is very slow.

    function output = initialize_zeros(mat_size)
        output = cell(1,1,zPoints);
        output(:) = {zeros(mat_size)};
    end

% =========================================================================
% Power
% =========================================================================
% Pump
if ismember(direction,{'forward','both'})
    Power_pump_forward = initialize_zeros(1);
end
if ismember(direction,{'backward','both'})
    Power_pump_backward = initialize_zeros(1);
end
% ASE
if ismember(direction,{'forward','both'})
    if gain_rate_eqn.ignore_ASE % Because ASE is ignored, set it a scalar zero is enough.
        Power_ASE_forward = initialize_zeros(1);
    else % include ASE
        % Make it the size of (1,1,num_modes,1,N) for "solve_gain_rate_eqn.m"
        Power_ASE_forward = initialize_zeros([1,1,num_modes,1,Nt]);
    end
end
if ismember(direction,{'backward','both'})
    if gain_rate_eqn.ignore_ASE % Because ASE is ignored, set it a scalar zero is enough.
        Power_ASE_backward = initialize_zeros(1);
    else % include ASE
        % Make it the size of (1,1,num_modes,1,N) for "solve_gain_rate_eqn.m"
        Power_ASE_backward = initialize_zeros([1,1,num_modes,1,Nt]);
    end
end

% -------------------------------------------------------------------------
% Put in the necessary information
% -------------------------------------------------------------------------
% Pump power
if ismember(direction,{'forward','both'})
    if segment_idx == 1 && ismember(gain_rate_eqn.pump_direction,{'co','bi'})
        Power_pump_forward{1} = gain_rate_eqn.copump_power;
    end
end
if ismember(direction,{'backward','both'})
    if segment_idx == num_segments && ismember(gain_rate_eqn.pump_direction,{'counter','bi'})
        Power_pump_backward{end} = gain_rate_eqn.counterpump_power;
    end
end

% ASE
if ismember(direction,{'forward','both'})
    if gain_rate_eqn.include_ASE
        % Make it the size of (1,1,num_modes,1,N) for "solve_gain_rate_eqn.m"
        Power_ASE_forward{1} = permute(initial_condition.Power.ASE.forward,[3 4 2 5 1]);
    end
end
if ismember(direction,{'backward','both'})
    if gain_rate_eqn.include_ASE
        % Make it the size of (1,1,num_modes,1,N) for "solve_gain_rate_eqn.m"
        Power_ASE_backward{end} = permute(initial_condition.Power.ASE.backward,[3 4 2 5 1]);
    end
end

% =========================================================================
% Signal fields (backward and forward)
% =========================================================================
% Forward signal field (the pulse)
if ismember(direction,{'forward','both'})
    signal_fields = initialize_zeros([Nt,num_modes]);
    if segment_idx == 1
        signal_fields{1} = initial_condition.fields;
    end
end
% Backward signal field
% (This signal_fields_backward is used only for the linear-oscillator scheme.)
if isequal(direction,'both')
    signal_fields_backward = initialize_zeros(1);
end

% =========================================================================
% Output arguments
% =========================================================================
switch direction
    case 'forward'
        varargout = {signal_fields,...
                     Power_pump_forward,...
                     Power_ASE_forward};
    case 'backward'
        varargout = {Power_pump_backward,...
                     Power_ASE_backward};
    otherwise % 'both'
        varargout = {signal_fields,signal_fields_backward,...
                     Power_pump_forward,Power_pump_backward,...
                     Power_ASE_forward,Power_ASE_backward};
end

end

%% CALC_TOTAL_ENERGY
function total_energy = calc_total_energy(Aw,Nt,dt)
%CALC_TOTAL_ENERGY
% total_energy = sum(trapz(abs(fftshift(Aw,1)).^2),2)*Nt*dt/1e3; % nJ
%
% "Aw" is the spectral intensity of the pulse under the convention of
% "ifft", so it needs "fftshift" to reorder it correctly for "trapz".
%
% In practice, because "fftshift" is slow, I use "sum" instead.
% The difference between "trapz" and "sum" is half the sum of the first and
% last elements.

center_idx = ceil(size(Aw,1)/2);
A2 = abs(Aw).^2;

total_energy = ( sum(A2(:))-sum(A2([center_idx,center_idx+1],:),[1,2])/2 )*Nt*dt/1e3; % nJ;

end

%% MYGPUARRAY
function varargout = mygpuArray(varargin)
%MYGPUARRAY It throws all the inputs to GPU
% Input arguments should be in "cell" arrays.

varargout = cell(1,nargin);
for i = 1:nargin
    varargout{i} = cellfun(@gpuArray,varargin{i},'UniformOutput',false);
end

end

%% MYGATHER
function varargout = mygather(varargin)
%MYGATHER It gathers all the inputs from GPU to the RAM
% Input arguments should be in "cell" arrays.

varargout = cell(1,nargin);
for i = 1:nargin
    varargout{i} = cellfun(@gather,varargin{i},'UniformOutput',false);
end

end

%% MYGPUARRAY2
function varargout = mygpuArray2(starti,endi,varargin)
%MYGPUARRAY2 It throws a part of the inputs to GPU, from "starti" to "endi"
% Input arguments should be in "cell" arrays.

varargout = varargin;
for i = 1:nargin-2
    x = cellfun(@gpuArray,varargin{i}(starti:endi),'UniformOutput',false);
    varargout{i}(starti:endi) = x;
end

end

%% MYGATHER2
function varargout = mygather2(starti,endi,varargin)
%MYGATHER2 It gathers a part of the inputs from GPU to the RAM, from "starti" to "endi"
% Input arguments should be in "cell" arrays.

varargout = varargin;
for i = 1:nargin-2
    x = cellfun(@gather,varargin{i}(starti:endi),'UniformOutput',false);
    varargout{i}(starti:endi) = x;
end

end

%% CLEANMEUP
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end