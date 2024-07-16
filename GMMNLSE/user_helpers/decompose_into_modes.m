function mode_time_profiles_tm = decompose_into_modes(use_gpu, normalized_mode_space_profiles_xym, full_field_txy, dx, cuda_dir_path)
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
        recompile_ptx(cuda_dir_path,'decompose_into_modes_Cartesian_basis.cu','decompose_into_modes_Cartesian_basis.ptx');
    end
    
    try
        kernel = parallel.gpu.CUDAKernel(fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.ptx'), fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.cu'));
    catch
        % Compile the CUDA code again.
        % Currently found error:
        %    version mismatch due to different versions of cuda I use in Windows and Debian.
        recompile_ptx(cuda_dir_path,'decompose_into_modes_Cartesian_basis.cu','decompose_into_modes_Cartesian_basis.ptx');
        kernel = parallel.gpu.CUDAKernel(fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.ptx'), fullfile(cuda_dir_path,'decompose_into_modes_Cartesian_basis.cu'));
    end
    kernel.ThreadBlockSize = [num_modes,1,1];
    kernel.GridSize =[Nt*Nx^2, 1, 1];
    mode_time_profiles_tm = complex(zeros(Nt,num_modes,'gpuArray'));
    mode_time_profiles_tm = feval(kernel,...
                                  mode_time_profiles_tm,...
                                  normalized_mode_space_profiles_xym,ones(1,num_modes,'gpuArray'),...
                                  permute(full_field_txy,[2,3,1]),...
                                  dx,dx,...
                                  num_modes,Nx,Nx,Nt);
    
    mode_time_profiles_tm = gather(mode_time_profiles_tm);
else
    normalized_mode_space_profiles_mxy = permute(normalized_mode_space_profiles_xym,[3 1 2]);
    
    % Einstein summation convention: tm=(txy)*(mxy)
    mode_time_profiles_tm = reshape(full_field_txy,[Nt,Nx^2])*reshape(normalized_mode_space_profiles_mxy,[num_modes,Nx^2]).'*dx^2;
end

end

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