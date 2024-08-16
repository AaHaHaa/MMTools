function kernel = setup_kernel_MPA_coeff(filename,cuda_dir_path,Nt,Nx,Ny,M,Np)

cudaFilename = [filename, '.cu'];
ptxFilename = [filename, '.ptx'];

if ~exist(fullfile(cuda_dir_path,ptxFilename), 'file')
    recompile_ptx(cuda_dir_path,cudaFilename,ptxFilename);
end

% Setup the kernel from the cu and ptx files
try
    kernel = parallel.gpu.CUDAKernel(fullfile(cuda_dir_path,ptxFilename), fullfile(cuda_dir_path,cudaFilename));
catch
    % Compile the CUDA code again.
    % Currently found error:
    %    version mismatch due to different versions of cuda I use in Windows and Debian.
    recompile_ptx(cuda_dir_path,cudaFilename,ptxFilename);
    kernel = parallel.gpu.CUDAKernel(fullfile(cuda_dir_path,ptxFilename), fullfile(cuda_dir_path,cudaFilename));
end

if M > kernel.MaxThreadsPerBlock
    error('setup_kernel_MPA_coeff:numThreadsError',...
          '(MPA.M+1) can''t be larger than %u.',uint32(kernel.MaxThreadsPerBlock));
end
if M > 30
    error('setup_kernel_MPA_coeff:MError',...
          '(MPA.M+1) is set to be 30 max in the cuda file.');
end

kernel.ThreadBlockSize = [M,1,1];
kernel.GridSize =[Nt*Nx*Ny*Np, 1];

end

%% Helper function
function recompile_ptx(cuda_dir_path,cudaFilename,ptxFilename)

if ispc
    system(['nvcc -ptx "', fullfile(cuda_dir_path,cudaFilename), '" --output-file "', fullfile(cuda_dir_path,ptxFilename) '"']);
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