function full_field_txy = recompose_into_space(use_gpu,normalized_mode_space_profiles_xym, mode_time_profiles_tm,cuda_dir_path)
% RECOMPOSE_INTO_SPACE Combine a set of mode time profiles with their 
% corresponding space profiles to get the full 3D spatio-temperal field.
%
% =========================================================================
% Input:
% use_gpu - true or false
% normalized_mode_space_profiles_xym - a (Nx, Nx, num_modes) matrix with each mode's profile in space. The units do not matter.
% mode_time_profiles_tm - a (Nt, num_modes) matrix with each mode's time profile.
% cuda_dir_path - if you use GPU, provide the cuda (recompose_into_space.cu) path
% -------------------------------------------------------------------------
% Output:
% full_field - a (Nt,Nx,Nx) matrix with full spatiotemporal fields in each time.
% =========================================================================

num_modes = size(normalized_mode_space_profiles_xym, 3);
Nx = size(normalized_mode_space_profiles_xym, 1);
Nt = size(mode_time_profiles_tm, 1);

if use_gpu
    % I use cuda for GPU computing because it's faster and more memory-efficient
    if ~exist(fullfile(cuda_dir_path,'recompose_into_space.ptx'), 'file')
        if ispc
            system(['nvcc -ptx "', fullfile(cuda_dir_path,'recompose_into_space.cu'), '" --output-file "', fullfile(cuda_dir_path,'recompose_into_space.ptx') '"']);
        else % unix
            % tested: Debian 9 (Stretch)
            % Cuda 8 doesn't support gcc6, beware to use gcc5 or clang-3.8.
            system(['nvcc -ccbin clang-3.8 -ptx "', fullfile(cuda_dir_path,'recompose_into_space.cu'), '" --output-file "', fullfile(cuda_dir_path,'recompose_into_space.ptx') '"']);
        end
    end
    kernel = parallel.gpu.CUDAKernel(fullfile(cuda_dir_path,'recompose_into_space.ptx'), fullfile(cuda_dir_path,'recompose_into_space.cu'));
    currentGPU = gpuDevice;
    kernel.ThreadBlockSize = [currentGPU.MaxThreadBlockSize(1),1,1];
    kernel.GridSize =[ceil(Nt*Nx^2/kernel.ThreadBlockSize(1)), 1];
    full_field_txy = complex(zeros(Nt,Nx,Nx,'gpuArray'));
    full_field_txy = feval(kernel,full_field_txy,normalized_mode_space_profiles_xym,mode_time_profiles_tm,ones(1,num_modes,'gpuArray'),num_modes,Nx,Nx,Nt);
    
    full_field_txy = gather(full_field_txy);
else
    normalized_mode_space_profiles_mxy = permute(normalized_mode_space_profiles_xym,[3 1 2]);

    % I use tensor product here, summing over index m.
    % Since MATLAB doesn't support tensor product, we need to use the trick
    % of "reshape" function with a combination of matrix multiplications.
    %
    % F_txy = sum( F_tm*P_mxy ,m)
    full_field_txy = reshape(...
        mode_time_profiles_tm*reshape(normalized_mode_space_profiles_mxy, [num_modes Nx^2]),...
                        [Nt Nx Nx]);
end

end