function mode_time_profiles_tm = decompose_into_modes(use_gpu,normalized_mode_space_profiles_xym, full_field_txy, dx, cuda_dir_path)
%DECOMPOSE_INTO_MODES Get the modal decomposition from the full 3D field
% mode_space_profile - a (Nx, Nx, num_modes) matrix with the mode profile 
% for each mode.
%
% =========================================================================
% Input:
% use_gpu - true or false
% normalized_mode_space_profiles_xym - a (Nx, Nx, num_modes) matrix with each mode's profile in space. The units do not matter.
% full_field - a (Nt,Nx,Nx) matrix with full spatiotemporal fields in each time.
% dx - spatial grid spacing, in m
% cuda_dir_path - if you use GPU, provide the cuda (decompose_into_modes_Cartesian_basis.cu) path
% -------------------------------------------------------------------------
% Output:
% mode_time_profiles - a (Nt, num_modes) matrix with each mode's time profile.
% =========================================================================

num_modes = size(normalized_mode_space_profiles_xym, 3);
Nx = size(normalized_mode_space_profiles_xym, 1);
Nt = size(full_field_txy, 1);

if use_gpu
    % I use cuda for GPU computing because it's faster and more memory-efficient
    if ~exist(fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.ptx'), 'file')
        if ispc
            system(['nvcc -ptx "', fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.cu'), '" --output-file "', fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.ptx') '"']);
        else % unix
            % tested: Debian 9 (Stretch)
            % Cuda 8 doesn't support gcc6, beware to use gcc5 or clang-3.8.
            system(['nvcc -ccbin clang -ptx "', fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.cu'), '" --output-file "', fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.ptx') '"']);
        end
    end
    kernel = parallel.gpu.CUDAKernel(fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.ptx'), fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.cu'));
    kernel.ThreadBlockSize = [num_modes,1,1];
    kernel.GridSize =[ceil(Nt*Nx^2/kernel.ThreadBlockSize(1)), 1, 1];
    mode_time_profiles_tm = complex(zeros(Nt,num_modes,'gpuArray'));
    mode_time_profiles_tm = feval(kernel,mode_time_profiles_tm,normalized_mode_space_profiles_xym,ones(1,num_modes,'gpuArray'),permute(full_field_txy,[2,3,1]),dx,dx,num_modes,Nx,Nx,Nt);
    
    mode_time_profiles_tm = gather(mode_time_profiles_tm);
else
    normalized_mode_space_profiles_mxy = permute(normalized_mode_space_profiles_xym,[3 1 2]);
    
    % Einstein summation convention: tm=(txy)*(mxy)
    mode_time_profiles_tm = reshape(full_field_txy,[Nt,Nx^2])*reshape(normalized_mode_space_profiles_mxy,[num_modes,Nx^2]).'*dx^2;
end

end