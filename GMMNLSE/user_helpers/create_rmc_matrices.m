function iQ = create_rmc_matrices(fiber,sim,num_modes,save_points)
%CREATE_RMC_MATRICES Summary of this function goes here
%   Detailed explanation goes here

sim.rmc.func_blkdiag_rand = blkdiag_rand(sim.gpu_yes);

mode_profiles = load_mode_profiles(fiber,sim);

% Find the random mode coupling strength
sim.rmc.std_Q = rmc_std_Q( fiber,sim,num_modes,mode_profiles );

iQ = 1i*sim.rmc.func_blkdiag_rand.hermitian(num_modes,save_points).*sim.rmc.std_Q;

if sim.gpu_yes
    iQ = gather(iQ);
end

end

%% helper function
function mode_profiles = load_mode_profiles(fiber,sim)

load(sprintf('%smode%uwavelength%u.mat',fiber.MM_folder,sim.midx(1),int16(sim.rmc.lambda0*1e10)),'phi','x'); % load 1st mode first to get the size and dimension of the mode profile
mode_profiles.x = x; % x-position vector
mode_profiles.mode_profiles = zeros(length(x),length(x),length(sim.midx)); % initialization
mode_profiles.mode_profiles(:,:,1) = phi; % the 1st mode
for ni = 2:length(sim.midx)
    n = sim.midx(ni);
    load(sprintf('%smode%uwavelength%u.mat',fiber.MM_folder,n,int16(sim.rmc.lambda0*1e10)),'phi');
    mode_profiles.mode_profiles(:,:,ni) = phi;
end

mode_profiles.dx = abs(mode_profiles.x(2)-mode_profiles.x(1)); % unit: um
% Normalization
norms = sqrt(sum(sum( abs(mode_profiles.mode_profiles).^2 ,1),2))*mode_profiles.dx;
mode_profiles.mode_profiles = mode_profiles.mode_profiles./norms; % unit: 1/um

% Downsample the spatial profiles.
total_mode_profiles = sum(mode_profiles.mode_profiles.^2,3);
max_mode_profiles = max(mode_profiles.mode_profiles(:).^2);
idx = total_mode_profiles>max_mode_profiles/1e4;
[~,chosen_region_left_idx] = ind2sub(size(total_mode_profiles),find(idx,1));         chosen_region_left_idx  = chosen_region_left_idx  - 1;
[~,chosen_region_right_idx] = ind2sub(size(total_mode_profiles),find(idx,1,'last')); chosen_region_right_idx = chosen_region_right_idx + 1;
core_region_idx = chosen_region_left_idx:chosen_region_right_idx;
mode_profiles.x = downsample(mode_profiles.x(core_region_idx),sim.rmc.downsampling_factor);
old_large_mode_profiles = mode_profiles.mode_profiles;
mode_profiles.mode_profiles = zeros(length(mode_profiles.x),length(mode_profiles.x),length(sim.midx));
for ni = 1:length(sim.midx)
    mode_profiles.mode_profiles(:,:,ni) = downsample(downsample(old_large_mode_profiles(core_region_idx,core_region_idx,ni), sim.rmc.downsampling_factor)', sim.rmc.downsampling_factor)'; % downsample;;
end

mode_profiles.dx = mode_profiles.dx*sim.rmc.downsampling_factor;

end