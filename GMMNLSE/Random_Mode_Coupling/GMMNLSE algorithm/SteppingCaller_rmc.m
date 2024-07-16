function [A_out,T_delay_out] = SteppingCaller_rmc(sim,...
                                                  G, saturation_parameter,...
                                                  num_zPoints, save_points, num_zPoints_persave,...
                                                  initial_condition,...
                                                  prefactor,...
                                                  SK_info, SRa_info, SRb_info,...
                                                  D,...
                                                  haw, hbw,...
                                                  haw_sponRS, hbw_sponRS, sponRS_prefactor)
%STEPPINGCALLER_RMC Summary of this function goes here
%   Detailed explanation goes here

Nt = size(initial_condition.fields,1);
num_modes = size(initial_condition.fields,2);
dt = initial_condition.dt;

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

last_A = ifft(initial_condition.fields);

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
    count_progress_bar = 1;
    num_progress_updates = 1000;
end

switch sim.gain_model
    case 0
        gain_str = 'nogain';
    case 1
        gain_str = 'MMGaussianGain'; % with MPA
end
GMMNLSE_func = str2func(['stepping_' sim.step_method '_', gain_str '_rmc']);

% Then start the propagation
for ii = 2:num_zPoints
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('GMMNLSE_propagate:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end

    last_A = GMMNLSE_func(last_A, dt,...
                          sim, prefactor,...
                          SK_info, SRa_info, SRb_info,...
                          D, sim.rmc.D(ii-1),...
                          haw, hbw,...
                          haw_sponRS, hbw_sponRS, sponRS_prefactor,...
                          G, saturation_parameter);
                      
    % Apply the damped frequency window
    last_A = last_A.*sim.damped_freq_window;

    % Check for any NaN elements, if desired
    if any(any(isnan(last_A))) %any(isnan(last_A),'all')
        error('SteppingCaller_rmc:NaNError',...
              'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
    end

    % Center the pulse
    if sim.pulse_centering
        last_A_in_time = fft(last_A);
        TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(last_A_in_time).^2,[1,2])/sum(abs(last_A_in_time).^2,[1,2]));
        % Because circshift is slow on GPU, I discard it.
        %last_result = ifft(circshift(last_A_in_time,-tCenter));
        if TCenter ~= 0
            if TCenter > 0
                last_A = ifft([last_A_in_time(1+TCenter:end,:);last_A_in_time(1:TCenter,:)]);
            elseif TCenter < 0
                last_A = ifft([last_A_in_time(end+1+TCenter:end,:);last_A_in_time(1:end+TCenter,:)]);
            end
            if sim.gpu_yes
                TCenter = gather(TCenter);
            end
            T_delay = T_delay + TCenter*dt;
        end
    end

    % If it's time to save, get the result from the GPU if necessary,
    % transform to the time domain, and save it
    if rem(ii-1, num_zPoints_persave) == 0
        A_out_ii = fft(last_A);

        if sim.gpu_yes
            A_out(:, :, int64((ii-1)/num_zPoints_persave)+1) = gather(A_out_ii);
        else
            A_out(:, :, int64((ii-1)/num_zPoints_persave)+1) = A_out_ii;
        end

        T_delay_out(int64((ii-1)/num_zPoints_persave)+1) = T_delay;
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

%% Helper function
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end