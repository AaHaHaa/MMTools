function full_field_txy = recompose_into_space(use_gpu, normalized_mode_space_profiles_xym, mode_time_profiles_tm, cuda_dir_path)
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
        recompile_ptx(cuda_dir_path,'recompose_into_space.cu','recompose_into_space.ptx');
    end
    
    try
        kernel = parallel.gpu.CUDAKernel(fullfile(cuda_dir_path,'recompose_into_space.ptx'), fullfile(cuda_dir_path,'recompose_into_space.cu'));
    catch
        % Compile the CUDA code again.
        % Currently found error:
        %    version mismatch due to different versions of cuda I use in Windows and Debian.
        recompile_ptx(cuda_dir_path,'recompose_into_space.cu','recompose_into_space.ptx');
        kernel = parallel.gpu.CUDAKernel(fullfile(cuda_dir_path,'recompose_into_space.ptx'), fullfile(cuda_dir_path,'recompose_into_space.cu'));
    end
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

%% helper function
function recompile_ptx(cuda_dir_path,cudaFilename,ptxFilename)

if ispc
    MATLAB_version = version('-release'); MATLAB_version = str2double(MATLAB_version(1:4));
    if MATLAB_version < 2023
        system(['nvcc -ptx "', fullfile(cuda_dir_path,cudaFilename), '" --output-file "', fullfile(cuda_dir_path,ptxFilename) '"']);
    else
        current_path = pwd;
        cd(cuda_dir_path);
        mexcuda('-ptx',cudaFilename);
        cd(current_path);
    end
else % unix
    % tested: Debian 10 (Buster)
    system(['nvcc -ccbin clang -ptx "', fullfile(cuda_dir_path,cudaFilename), '" --output-file "', fullfile(cuda_dir_path,ptxFilename) '"']);

    % For Yi-Hao's Debian PC only.
    % Matlab 2020b generates ptx of ".target sm_52" by default.
    % However, my own GPU, GeForce 660 Ti" is too old. I need to set it to sm_30 to make it work.
    % Starting from Matlab 2021a, my own GPU isn't supported.
    system(['sed -i ''s/.target sm_52/.target sm_30/g'' "' fullfile(cuda_dir_path,ptxFilename) '"']); % run the shell script to change it to sm_30
end

end