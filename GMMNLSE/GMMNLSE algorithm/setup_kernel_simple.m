function kernel = setup_kernel_simple(filename,cuda_dir_path,total_num)

total_num = gather(total_num);

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

if total_num > kernel.MaxThreadsPerBlock
    ThreadsPerBlock = kernel.MaxThreadsPerBlock;
else
    ThreadsPerBlock = total_num;
end
kernel.ThreadBlockSize = [ThreadsPerBlock,1,1];
kernel.GridSize =[ceil(total_num/ThreadsPerBlock), 1];

end

%% Helper function
function recompile_ptx(cuda_dir_path,cudaFilename,ptxFilename)

if ispc
    MATLAB_version = version('-release'); MATLAB_version = str2double(MATLAB_version(1:4));
    if MATLAB_version < 2023
        system(['nvcc -ptx "', fullfile(cuda_dir_path,cudaFilename), '" --output-file "', fullfile(cuda_dir_path,ptxFilename) '"']);
    else
        try
            current_path = pwd;
            cd(cuda_dir_path);
            mexcuda('-ptx',cudaFilename);
            cd(current_path);
        catch ME
            %disp(ME.message);
            system(['nvcc -ptx "', fullfile(cuda_dir_path,cudaFilename), '" --output-file "', fullfile(cuda_dir_path,ptxFilename) '"']);
            cd(current_path);
        end
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