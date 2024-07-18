function [signal_fields_out,T_delay_out,...
          Power_out,N,...
          saved_data] = SteppingCaller_rategain_linear_oscillator(sim,gain_rate_eqn,...
                                                                  num_zPoints,save_points,num_zPoints_persave,...
                                                                  initial_condition,...
                                                                  prefactor,...
                                                                  SK_info, SRa_info,SRb_info,...
                                                                  omegas,D_op,...
                                                                  haw,hbw,...
                                                                  haw_sponRS,hbw_sponRS,sponRS_prefactor,...
                                                                  saved_data)
%STEPPINGCALLER_RATEGAIN_LINEAR_OSCILLATOR It attains the field after 
%propagation inside the gain medium solved by the rate equations.
%   
% The computation of this code is based on
%   1. Lindberg et al., "Accurate modeling of high-repetition rate ultrashort pulse amplification in optical fibers", Scientific Reports (2016)
%   2. Chen et al., "Optimization of femtosecond Yb-doped fiber amplifiers for high-quality pulse compression", Opt. Experss (2012)
%   3. Gong et al., "Numerical modeling of transverse mode competition in strongly pumped multimode fiber lasers and amplifiers", Opt. Express (2007)
%
%   Please go to "gain_info.m" file to find the information about some input arguments.
%   The info of some other input arguments are inside "GMMNLSE_propagate.m"
%
% Output:
%   signal_fields_out - the signal field; (N,num_modes,save_points)
%   Power_out - a structure with
%               pump.forward - (1,1,save_points)
%               pump.backward - (1,1,save_points)
%               ASE.forward - (N,num_modes,save_points)
%               ASE.backward - (N,num_modes,save_points)
%   N - the ion density/population of each energy level; (Nx,Nx,save_points)
%       the output of N is transformed into N/N_total, the ratio of
%       population inversion

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
    segments = [zPoints_each_segment*ones(1,num_segments-1) num_zPoints-zPoints_each_segment*(num_segments-1)];
    if segments(end) == 0
        segments = segments(1:end-1);
    end
end

%% Propagations
% =========================================================================
% Load the saved backward data
% =========================================================================
% Only the counterpumping pump power at the pulse-input end is used
% If saved backward pump comes from "GMMNLSE_rategain()", 
%       it's a cell array of counterpump power at all z points.
% If saved backward pump comes from this function "GMMNLSE_rategain_linear_oscillator()",
%       it's a scalar value of counterpump power at the pulse-input end.
% Therefore, I use a "extract_first_content()" to deal with this discrepancy.
%
% For other two saved data, backward ASE and backward signal field,
%       they're always cell arrays of data at all z points.
guess_Power_pump_backward = extract_first_content(saved_data.Power_pump_backward);
Power_ASE_backward        = saved_data.Power_ASE_backward;
signal_fields_backward    = saved_data.signal_fields_backward;

% This linear-oscillator algorithm updates only forward signal field,
% forward and backward pump powers, and forward ASE power.
% Since pump powers are saved only at output saved points, only backward
% ASE and backward signal field are needed to put into GPU for performance
% if GPU is used.
if sim.gpu_yes
    [signal_fields_backward,Power_ASE_backward] = mygpuArray2(1,segments(1),...
                                                              signal_fields_backward,Power_ASE_backward);
end

% =========================================================================
% Start the pulse propagation
% =========================================================================
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
    [signal_fields,signal_fields_backward,...
     T_delay_out,...
     Power_pump_forward,Power_pump_backward,...
     Power_ASE_forward,Power_ASE_backward,...
     N] = gain_propagate(sim,gain_rate_eqn,...
                          num_zPoints,segments,save_points,num_zPoints_persave,...
                          Nt,prefactor,SK_info,SRa_info,SRb_info,omegas,D_op,...
                          haw,hbw,...
                          haw_sponRS,hbw_sponRS,sponRS_prefactor,...
                          initial_condition,guess_Power_pump_backward,...
                          signal_fields_backward,Power_ASE_backward);

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
    while abs(Power_pump_backward{end}-gain_rate_eqn.counterpump_power)/gain_rate_eqn.counterpump_power > min(gain_rate_eqn.tol,0.01) && iterations <= gain_rate_eqn.max_iterations
        % I call this "modified" binary search because I don't use the
        % middle point as the next upper or lower boundary value; instead,
        % I use the linear interpolation to find it. It's just faster!
        % Traditional binary serarch uses
        %    binary_m = (binary_L + binary_R)/2;
        binary_m = binary_L*(binary_R_counterpump_power-gain_rate_eqn.counterpump_power)/(binary_R_counterpump_power-binary_L_counterpump_power) + binary_R*(gain_rate_eqn.counterpump_power-binary_L_counterpump_power)/(binary_R_counterpump_power-binary_L_counterpump_power);
        
        [signal_fields,signal_fields_backward,...
         T_delay_out,...
         Power_pump_forward,Power_pump_backward,...
         Power_ASE_forward,Power_ASE_backward,...
         N] = gain_propagate(sim,gain_rate_eqn,...
                              num_zPoints,segments,save_points,num_zPoints_persave,...
                              Nt,prefactor,SK_info,SRa_info,SRb_info,omegas,D_op,...
                              haw,hbw,...
                              haw_sponRS,hbw_sponRS,sponRS_prefactor,...
                              initial_condition,binary_m,...
                              signal_fields_backward,Power_ASE_backward);
        
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
saved_zPoints = 1:num_zPoints_persave:num_zPoints;

% Change the size back to (N,num_modes)
Power_ASE_forward_out   = cellfun(@(P) permute(P,[5 3 1 2 4]), Power_ASE_forward ,'UniformOutput',false);
Power_ASE_backward_out  = cellfun(@(P) permute(P,[5 3 1 2 4]), Power_ASE_backward,'UniformOutput',false);

% Transform them into arrays
signal_fields_out          = fft(cell2mat(signal_fields         (:,:,saved_zPoints)));
signal_fields_backward_out = fft(cell2mat(signal_fields_backward(:,:,saved_zPoints)));
Power_pump_forward_out     = cell2mat(Power_pump_forward);
Power_pump_backward_out    = cell2mat(Power_pump_backward);
Power_ASE_forward_out      = fftshift(cell2mat(Power_ASE_forward_out (:,:,saved_zPoints)),1);
Power_ASE_backward_out     = fftshift(cell2mat(Power_ASE_backward_out(:,:,saved_zPoints)),1);   

signal_fields_out = struct('forward', signal_fields_out,...
                           'backward',signal_fields_backward_out);
Power_out = struct('pump',struct('forward',Power_pump_forward_out,'backward',Power_pump_backward_out),...
                   'ASE', struct('forward', Power_ASE_forward_out,'backward', Power_ASE_backward_out));

% Reverse the order and save the data for the linear oscillator scheme for the next iteration
reverse_direction = @(x) flip(x,3);
Power_pump_backward = Power_pump_forward{end};
Power_ASE_backward  = reverse_direction(Power_ASE_forward);
signal_fields_backward = reverse_direction(signal_fields);
saved_data.Power_pump_backward    = Power_pump_backward;
saved_data.Power_ASE_backward     = Power_ASE_backward;
saved_data.signal_fields_backward = signal_fields_backward;

end

%%
function [signal_fields,signal_fields_backward,...
          T_delay_out,...
          Power_pump_forward,Power_pump_backward,...
          Power_ASE_forward,Power_ASE_backward,...
          N] = gain_propagate(sim,gain_rate_eqn,...
                               num_zPoints,segments,save_points,num_zPoints_persave,...
                               Nt,prefactor,SK_info,SRa_info,SRb_info,omegas,D_op,...
                               haw,hbw,...
                               haw_sponRS,hbw_sponRS,sponRS_prefactor,...
                               initial_condition,input_Power_pump_backward,...
                               signal_fields_backward,Power_ASE_backward)
%GAIN_PROPAGATE Runs the corresponding propagation method.

dt = initial_condition.dt;

T_delay_out = zeros(save_points,1);

zPoints_each_segment = segments(1);
num_modes = size(initial_condition.fields,2);

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
        T_delay = TCenter*initial_condition.dt;
    else
        T_delay = 0;
    end
else
    T_delay = 0;
end
T_delay_out(1) = T_delay;
initial_condition.fields = ifft(initial_condition.fields);

% Initialization
[signal_fields,...
 Power_pump_forward,Power_pump_backward,...
 Power_ASE_forward] = initialization(sim,gain_rate_eqn,...
                                     num_zPoints,...
                                     Nt,num_modes,...
                                     save_points,...
                                     segments(1),...
                                     initial_condition,input_Power_pump_backward);

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
    count_progress_bar = 1;
    num_progress_updates = 1000;
end

% Then start the propagation
last_N = cat(8, ones([size(gain_rate_eqn.N_total),1,1,1,1,1,1]).*gain_rate_eqn.N_total,...
               zeros([size(gain_rate_eqn.N_total),1,1,1,1,1,length(gain_rate_eqn.energy_levels)-2])); % initial guess for solving the population during propagation
switch sim.step_method
    case 'RK4IP'
        mode_str = '';
    case 'MPA'
        mode_str = 'MM';
end
GMMNLSE_rategain_func = str2func(['stepping_',sim.step_method,'_',mode_str,'rategain_linear_oscillator']);
for ii = 2:num_zPoints
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('GMMNLSE_propagate:ProgressBarBreak',...
        'The "cancel" button of the progress bar has been clicked.');
    end
    
    % =====================================================================
    % GMMNLLSE: Run the correct step function depending on the options chosen.
    % =====================================================================     
    Zi = ii;

    % Load the initial powers and field
    if Zi == 2 % the first/starting z_point
        last_Power_pump_forward  = Power_pump_forward {1};
        last_Power_pump_backward = Power_pump_backward{1};
        last_Power_ASE_forward   = Power_ASE_forward  {1};
        last_signal_fields       = signal_fields{1};
    end
    last_Power_ASE_backward     = Power_ASE_backward    {Zi-1};
    last_signal_fields_backward = signal_fields_backward{Zi-1};

    [last_signal_fields,...
     last_Power_pump_forward,last_Power_pump_backward,...
     last_Power_ASE_forward,...
     last_N] = GMMNLSE_rategain_func(last_signal_fields,last_signal_fields_backward,...
                                      last_Power_pump_forward,last_Power_pump_backward,...
                                      last_Power_ASE_forward,last_Power_ASE_backward,...
                                      last_N,...
                                      dt,sim,prefactor,...
                                      SK_info,SRa_info,SRb_info,...
                                      omegas,D_op,...
                                      haw,hbw,...
                                      haw_sponRS,hbw_sponRS,sponRS_prefactor,...
                                      gain_rate_eqn);
    % Apply the damped frequency window
    last_signal_fields = last_signal_fields.*sim.damped_freq_window;

    % Update the pump powers only at saved points
    if rem(Zi-1, num_zPoints_persave) == 0
        Power_pump_forward {int64((Zi-1)/num_zPoints_persave+1)} = last_Power_pump_forward;
        Power_pump_backward{int64((Zi-1)/num_zPoints_persave+1)} = last_Power_pump_backward;
    end
    % Update only "forward" components at all z points
    Power_ASE_forward{Zi} = last_Power_ASE_forward;
    signal_fields    {Zi} = last_signal_fields;

    % Save N
    if Zi == num_zPoints % save the last N
        [~,~,~,~,...
         ~,last_N] = solve_gain_rate_eqn_linear_oscillator(sim,gain_rate_eqn,...
                                                           last_N,...
                                                           last_signal_fields,last_signal_fields_backward,...
                                                           last_Power_pump_forward,last_Power_pump_backward,...
                                                           last_Power_ASE_forward,last_Power_ASE_backward,...
                                                           omegas,dt);

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
        if sim.gpu_yes % if using MPA, save only the last one
            N(:,:,int64((Zi-2)/num_zPoints_persave+1),:) = permute(gather(last_N),[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
        else
            N(:,:,int64((Zi-2)/num_zPoints_persave+1),:) = permute(last_N,[1,2,3,8,5,6,7,4]); % (Nx,Nx,1,num_levels-1)
        end
    end

    % Save powers and fields
    if rem(Zi,zPoints_each_segment) == 0 || Zi == num_zPoints
        ready_for_next_segment_into_GPU = true;
    else
        ready_for_next_segment_into_GPU = false;
    end

    % =====================================================================
    % Some post-stepping checks, saves, and updates
    % =====================================================================
    % Check for any NaN elements
    if any(any(isnan(last_signal_fields))) %any(isnan(last_signal_fields),'all')
        error('GMMNLSE_propagate:NaNError',...
            'NaN field encountered, aborting.\nPossible reason is that the nonlinear length is too close to the large step size.');
    end
    
    % Center the pulse
    if sim.pulse_centering
        last_signal_fields_in_time = fft(last_signal_fields);
        tCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(last_signal_fields_in_time).^2,[1,2])/sum(abs(last_signal_fields_in_time).^2,[1,2]));
        if ~isnan(tCenter) && tCenter ~= 0 % all-zero fields; for calculating ASE power only
            % Because circshift is slow on GPU, I discard it.
            %last_signal_fields = ifft(circshift(last_signal_fields_in_time,-tCenter));
            if tCenter > 0
                last_signal_fields = ifft([last_signal_fields_in_time(1+tCenter:end,:);last_signal_fields_in_time(1:tCenter,:)]);
            elseif tCenter < 0
                last_signal_fields = ifft([last_signal_fields_in_time(end+1+tCenter:end,:);last_signal_fields_in_time(1:end+tCenter,:)]);
            end
            if sim.gpu_yes
                T_delay = T_delay + gather(tCenter*dt);
            else
                T_delay = T_delay + tCenter*dt;
            end
        end
        if rem(Zi-1, num_zPoints_persave) == 0
            T_delay_out(int64((Zi-1)/num_zPoints_persave)+1) = T_delay;
        end
    end
    
    % Put current GPU data back to RAM and those in the next segment from RAM to GPU
    if ready_for_next_segment_into_GPU
        current_segment_idx = ceil(Zi/zPoints_each_segment);
        next_segment_idx = current_segment_idx + 1;
        
        cumsum_segments = [0 cumsum(segments)];
        starti = cumsum_segments(current_segment_idx)+1;
        endi   = cumsum_segments(current_segment_idx+1);
        [signal_fields,signal_fields_backward,Power_ASE_forward,Power_ASE_backward] = mygather2(starti,endi,signal_fields,signal_fields_backward,Power_ASE_forward,Power_ASE_backward);
        
        if next_segment_idx <= length(segments)
            starti = cumsum_segments(next_segment_idx)+1;
            endi   = cumsum_segments(next_segment_idx+1);
            [signal_fields_backward,Power_ASE_backward] = mygpuArray2(starti,endi,signal_fields_backward,Power_ASE_backward);
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

% Output
[Power_pump_forward,Power_pump_backward] = mygather(Power_pump_forward,Power_pump_backward);

end

%% initialization
function [signal_fields,...
          Power_pump_forward,Power_pump_backward,...
          Power_ASE_forward] = initialization(sim,gain_rate_eqn,...
                                              first_segment,...
                                              Nt,num_modes,...
                                              save_points,...
                                              segment1,...
                                              initial_condition,input_Power_pump_backward)
%INITIALIZATION initializes "signal_fields" and "Powers" based on
%"segment_idx/num_segment".
%
%   They include copump_power, counterpump_power, initial_fields, and initial forward ASE.
%
%   The reason I use cell arrays instead of a matrix (N,num_modes,zPoints) for signal_fields and Power:
%       It's faster!
%       e.g. "signal_fields(:,:,zi) = signal_fields_next" is very slow.

    function output = initialize_zeros(mat_size,zPoints)
        output = cell(1,1,zPoints);
        output(:) = {zeros(mat_size)};
    end

% =========================================================================
% Initialize with zeros
% =========================================================================
% Signal fields
signal_fields = initialize_zeros([Nt,num_modes],first_segment);
% Pump
Power_pump_forward  = initialize_zeros(1,save_points);
Power_pump_backward = initialize_zeros(1,save_points);
% ASE
if gain_rate_eqn.ignore_ASE % Because ASE is ignored, set it a scalar zero is enough.
    Power_ASE_forward = initialize_zeros(1,first_segment);
else % include ASE
    Power_ASE_forward = initialize_zeros([1,1,num_modes,1,Nt],first_segment); % Make it the size to (1,1,num_modes,1,N) for "solve_gain_rate_eqn.m"
end

% =========================================================================
% Put in the necessary information
% =========================================================================
% Signal fields
signal_fields{1} = initial_condition.fields;
% Pump power
if any(strcmp(gain_rate_eqn.pump_direction,{'co','bi'}))
    Power_pump_forward{1} = gain_rate_eqn.copump_power;
end
if any(strcmp(gain_rate_eqn.pump_direction,{'counter','bi'}))
    Power_pump_backward{1} = input_Power_pump_backward;
end
% ASE power
if gain_rate_eqn.include_ASE
    Power_ASE_forward{1} = permute(initial_condition.Power.ASE.forward,[3 4 2 5 1]);
end

% =========================================================================
% Put necessary initialized data into GPU
% =========================================================================
if sim.gpu_yes
    [Power_pump_forward,Power_pump_backward] = mygpuArray(Power_pump_forward,Power_pump_backward);
    [signal_fields,Power_ASE_forward] = mygpuArray2(1,segment1,signal_fields,Power_ASE_forward);
end

end

%% EXTRACT_CELL_CONTENT
function x = extract_first_content(x)

if iscell(x)
    x = x{1};
end

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