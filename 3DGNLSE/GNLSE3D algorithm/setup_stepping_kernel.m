function [gpuDevice_Device,...
          cuda_MPA_psi_update] = setup_stepping_kernel(sim,Nt,Nx,Ny,M,Np)
%SETUP_STEPPING_KERNEL It sets cuda for computing sums of SR and SK terms.

% Use the specified GPU
% This needs to run at the beginning; otherwise, the already-stored values
% in GPU will be unavailable in a new GPU if the GPU device is switched.
try
    gpuDevice_Device = gpuDevice(sim.gpuDevice.Index); % use the specified GPU device
catch
    error('Please set the GPU you''re going to use by setting "sim.gpuDevice.Index".');
end

% Kernal for updating psi in MPA stepping algorithm
cuda_MPA_psi_update = setup_kernel_MPA_coeff('MPA_psi_update',sim.cuda_dir_path,Nt,Nx,Ny,M,Np);

end