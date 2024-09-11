function [E_out,...
          save_z,save_dz,...
          T_delay_out] = SteppingCaller_adaptive(sim,...
                                                 save_z, save_points,...
                                                 initial_condition,...
                                                 prefactor,...
                                                 D_op, W_op,...
                                                 fr, haw, hbw, sponRS_prefactor)
%STEPPINGCALLER_ADAPTIVE It starts the pulse propagation.

[Nt,Nx,Ny,~,Np] = size(initial_condition.field);
dt = initial_condition.dt;

save_dz = zeros(save_points,1);
T_delay_out = zeros(save_points,1);

% Pulse centering based on the moment of its intensity
if sim.pulse_centering
    % Center the pulse
    TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(initial_condition.field).^2,[1,2,3,5])/sum(abs(initial_condition.field(:)).^2));
    % Because circshift is slow on GPU, I discard it.
    if TCenter ~= 0
        if TCenter > 0
            initial_condition.field = [initial_condition.field(1+TCenter:end,:,:,:,:);initial_condition.field(1:TCenter,:,:,:,:)];
        elseif TCenter < 0
            initial_condition.field = [initial_condition.field(end+1+TCenter:end,:,:,:,:);initial_condition.field(1:end+TCenter,:,:,:,:)];
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

E_out = zeros(Nt, Nx, Ny, save_points, Np);

% Start by saving the initial condition
if sim.gpu_yes
    E_out(:,:,:,1,:) = gather(initial_condition.field);
else
    E_out(:,:,:,1,:) = initial_condition.field;
end

last_E = ifft(ifft(ifft(initial_condition.field,[],1),[],2),[],3); % in k- and frequency space

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

% Then start the propagation
while z+eps(z) < save_z(end) % eps(z) here is necessary due to the numerical error
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('GMMNLSE_propagate:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end

    ever_fail = false;
    previous_E = last_E;
    previous_a5 = a5;

    success = false;
    while ~success
        if ever_fail
            last_E = previous_E;
            a5 = previous_a5;
        end

        [last_E,a5,...
         opt_dz, success] = stepping_RK4IP_nogain_adaptive(last_E,...
                                                           sim, prefactor,...
                                                           D_op, W_op,...
                                                           fr, haw, hbw, sponRS_prefactor, dt,...
                                                           a5);

        if ~success
            ever_fail = true;

            sim.dz = opt_dz;
        end
    end
    sim.last_dz = sim.dz; % previous dz
    
    % Apply the damped frequency window
    last_E = last_E.*sim.damped_window;
    
    % Check for any NaN elements
    if any(any(isnan(last_E))) %any(isnan(last_A),'all')
        error('SteppingCaller_adaptive:NaNError',...
              'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
    end

    % Center the pulse
    if sim.pulse_centering
        last_E_in_time = fft(last_E,[],1);
        TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*abs(last_E_in_time).^2,[1,2,3,5])/sum(abs(last_E_in_time(:)).^2));
        % Because circshift is slow on GPU, I discard it.
        if TCenter ~= 0
            if ~isempty(a5) % RK4IP reuses a5 from the previous step
                a5 = fft(a5);
            end
            if TCenter > 0
                if ~isempty(a5) % RK4IP reuses a5 from the previous step
                    a5 = ifft([a5(1+TCenter:end,:,:,:,:);a5(1:TCenter,:,:,:,:)]);
                end
                last_E = ifft([last_E_in_time(1+TCenter:end,:,:,:,:);last_E_in_time(1:TCenter,:,:,:,:)],[],1);
            elseif TCenter < 0
                if ~isempty(a5) % RK4IP reuses a5 from the previous step
                    a5 = ifft([a5(end+1+TCenter:end,:,:,:,:);a5(1:end+TCenter,:,:,:,:)]);
                end
                last_E = ifft([last_E_in_time(end+1+TCenter:end,:,:,:,:);last_E_in_time(1:end+TCenter,:,:,:,:)],[],1);
            end
            if sim.gpu_yes
                TCenter = gather(TCenter);
            end
            T_delay = T_delay + TCenter*dt;
        end
    end

    % Update z
    z = z + sim.dz;

    sim.dz = min([opt_dz,save_z(end)-z,sim.adaptive_dz.max_dz]);

    % If it's time to save, get the result from the GPU if necessary,
    % transform to the time domain, and save it
    if z >= save_z(save_i)-eps(z)
        E_out_ii = fft(fft(fft(last_E,[],1),[],2),[],3);
        if sim.gpu_yes
            save_dz(save_i) = gather(sim.last_dz);
            save_z(save_i) = gather(z);
            E_out(:,:,:,save_i,:) = gather(E_out_ii);
        else
            save_dz(save_i) = sim.last_dz;
            save_z(save_i) = z;
            E_out(:,:,:,save_i,:) = E_out_ii;
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
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end