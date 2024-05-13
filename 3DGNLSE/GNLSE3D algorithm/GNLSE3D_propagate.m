function foutput = GNLSE3D_propagate(fiber, initial_condition, sim)
%GMMNLSE_PROPAGATE 

%%
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end

% Load the folder
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
current_folder = current_path(1:sep_pos(end-1));

sim.cuda_dir_path = [current_folder 'cuda'];

%% Consider only the last field of the initial condition
initial_condition.field = initial_condition.field(:,:,:,end,:); % take only the last field; size: (Nt,Nx,Ny,Nz,Np)

%% Check the validity of input parameters
if sim.save_period == 0
    sim.save_period = fiber.L0;
end

num_saves_total = fiber.L0/sim.save_period;
if rem(num_saves_total,1) && rem(num_saves_total+eps(num_saves_total),1) && rem(num_saves_total-eps(num_saves_total),1)
    error('GMMNLSE_propagate:SizeIncommemsurateError',...
          'The save period is %f m and the fiber length is %f m, which are not commensurate', sim.save_period, fiber.L0)
else
    num_saves_total = round(num_saves_total);
end

%% Some Parameters
% Get the numerical parameters from the initial condition.
[Nt, Nx, Ny, ~, ~] = size(initial_condition.field);

%% Set up the GPU details
try
    sim.gpuDevice.Device = gpuDevice(sim.gpuDevice.Index); % use the specified GPU device
catch
    error('Please set the GPU you''re going to use by setting "sim.gpuDevice.Index".');
end

%% Pre-calculate the dispersion term
% The "omegas" here is actually (omega - omega0), omega: true angular frequency
%                                                 omega0: central angular frequency (=2*pi*f0)
if Nt == 1 % CW case
    dt = 0;
    omegas = 0;
    omegas_real = 2*pi*sim.f0;
else
    if sim.gpu_yes
        dt = gpuArray(initial_condition.dt);
    else
        dt = initial_condition.dt;
    end
    omegas = 2*pi*ifftshift(linspace(-floor(Nt/2), floor((Nt-1)/2), Nt))'/(Nt*dt); % in 1/ps, in the order that the fft gives
    omegas_real = omegas + 2*pi*sim.f0;
end

dx = initial_condition.dx;
if isfield(initial_condition,'dy')
    dy = initial_condition.dy;
else
    dy = dx;
end
if sim.gpu_yes
    dx = gpuArray(dx);
    dy = gpuArray(dy);
end

[D_op,W_op,...
 kc,k0,...
 fiber.n,...
 sim] = calc_D_op(sim,...
                  fiber.n,...
                  Nt,dt,...
                  Nx,dx,...
                  Ny,dy,...
                  omegas,omegas_real,...
                  initial_condition.field);

%% Pre-calculate the factor used in 3D-GNLSE
prefactor = (1i/2./kc).*(k0.^2*2.*fiber.n*fiber.n2/3);

%% Pre-compute the Raman response in frequency space
if Nt == 1 % ignore Raman scattering under CW cases
    sim.Raman_model = 0;
end

if ~isfield(fiber,'material')
    fiber.material = 'silica';
end
if sim.Raman_model == 0 % no Raman
    haw = []; hbw = [];
    fr = 0;
else
    [fiber,haw,hbw] = Raman_model(fiber,sim,Nt,dt);
    fr = fiber.fr;
end

%% Include the shot noise: one photon per mode
if Nt ~= 1 % CW case
    initial_condition.field = include_shot_noise(sim,omegas,Nt,dt,dx,dy,initial_condition.field);
end

%% Include spontaneous Raman scattering
if sim.Raman_model ~= 0 && sim.Raman_sponRS
    sponRS_prefactor = spontaneous_Raman(sim,fr,...
                                         Nt,dt,...
                                            dx,dy);
else
    sponRS_prefactor = 0; % dummy variable
end

%% Create a damped frequency window to prevent aliasing
sim.damped_window = create_damped_window(Nt,Nx,Ny);

%% Setup the exact save points
% We will always save the initial condition as well
save_points = int64(num_saves_total + 1);
save_z = double(0:save_points-1)'*sim.save_period;

%% Run the step function over each step
run_start = tic;
% -------------------------------------------------------------------------
if sim.gain_model == 2 % rate-equation-gain model
    error('It hasn''t been implemented yet.')
% -------------------------------------------------------------------------
else % No gain, Gaussian gain
    [E_out,...
     save_z,save_deltaZ,...
     T_delay_out] = SteppingCaller_adaptive(sim,...
                                            save_z,save_points,...
                                            initial_condition,...
                                            prefactor,...
                                            D_op, W_op,...
                                            fr, haw, hbw, sponRS_prefactor);
end

% -------------------------------------------------------------------------
% Just to get an accurate timing, wait before recording the time
if sim.gpu_yes
    sim.betas = gather(sim.betas);
    dx = gather(dx);
    dy = gather(dy);
    wait(sim.gpuDevice.Device);
end
fulltime = toc(run_start);

%% Save the results in a struct
foutput = struct('z', save_z,...
                 'deltaZ', save_deltaZ,...
                 'field', E_out,...
                 'dt', initial_condition.dt,...
                 'dx', dx,...
                 'dy', dy,...
                 'betas', sim.betas,...
                 'seconds', fulltime,...
                 't_delay', T_delay_out);
if sim.gain_model == 2
    foutput.Power = Power;
    if gain_rate_eqn.export_N2
        foutput.N2 = N2;
    end
end

end