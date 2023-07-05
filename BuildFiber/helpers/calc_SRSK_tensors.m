%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script calculates the mode overlap intergrals for the Kerr
% nonlinearity, as in "Multimode Nonlinear Fibre Optics: Theory and
% Applications," P. Horak and F. Poletti
%
% If using a GPU, the CUDA toolkit needs to be installed, which requires
% visual studio as well, and the executables for the C++ compiler and the
% CUDA compiler driver need to be in the path. See the user manual for more
% details
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_SRSK_tensors(modes_used,Nx,gpu_yes,folder_name,lambda0)

    % File name parameters:
    dir_prefix = ['Fibers/', folder_name]; % folder containing the calculated modes

    %% Load the modes
    num_modes = length(modes_used); % number of modes for which the tensors should be calculated

    SR = zeros(num_modes^4, 1);

    if gpu_yes
        gd = gpuDevice();
        reset(gd); % It's always a good idea to reset the GPU before using it

        SR = gpuArray(SR);
    else % CPU
        if ispc
            [~,systemview] = memory;
            mem_planned_for_SRSK = systemview.PhysicalMemory.Available/3*2^20;
        elseif isunix
            [~,w] = unix('free -b | grep Mem'); % Display the memory info in Bytes
            stats = str2double(regexp(w, '[0-9]*', 'match'));
            %memsize = stats(1)/1e6;
            freemem = stats(end); % availabel memory
            mem_planned_for_SRSK = freemem/3;
        end
    end

    % Load each mode, and calculate the normalization constant ahead of time
    fields = zeros(Nx, Nx, num_modes); % the spatial field of each mode
    norms = zeros(num_modes, 1); % the normalization constant of each mode
    if gpu_yes
        fields = gpuArray(fields);
        norms = gpuArray(norms);
    end

    for ii = 1:num_modes
        l_str = num2str(round(lambda0*1e10));
        name = [dir_prefix '/mode',int2str(modes_used(ii)),'wavelength', l_str, '.mat'];
        load(name, 'phi');
        fields(:,:,ii) = phi;
        norms(ii) = sqrt(sum(sum(abs(phi).^2)));
        if norms(ii) == 0 % to avoid NaN when dividing norms, use 1 intentionally
            norms(ii) = 1;
        end
        disp(['Loaded mode ', int2str(modes_used(ii)), ' under ', l_str ' Ang'])
    end

    % Also load the spatial information to calculate Aeff accurately
    load(name, 'x');
    dx = (x(2)-x(1))*1e-6; % spatial step in m

    %% Calculate the overlap integrals
    % SR will hold the tensor. We only need to calculate SR

    % If using the GPU, we need to do some things differently
    if gpu_yes
        fields = permute(fields, [3 1 2]); % The order needs to go (num_modes, Nx, Nx)

        specific_filename = 'cuda/calculate_tensors_double';

        cudaFilename = [specific_filename, '.cu'];
        ptxFilename = [specific_filename, '.ptx'];

        if ~exist(ptxFilename, 'file')
            recompile_ptx(cudaFilename,ptxFilename);
        end

        % Setup the kernel from the cu and ptx files
        try
            kernel = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );
        catch
            % Compile the CUDA code again.
            % Currently found error:
            %    version mismatch due to different versions of cuda I use in Windows and Debian.
            recompile_ptx(cudaFilename,ptxFilename);
            kernel = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );
        end

        % Setup the kernel parameters
        num_threads_per_block = gd.MaxThreadBlockSize(1); % use as many threads per block as possible
        num_blocks = ceil((num_modes^4)/num_threads_per_block); % use as many blocks as needed
        kernel.ThreadBlockSize = [num_threads_per_block,1,1];
        kernel.GridSize = [num_blocks,1,1];

        % Run the CUDA code
        SR = feval(kernel, SR, fields, norms, int32(num_modes), int32(Nx));
        %wait(gd);
    else
        % If we're not using the GPU, then do all the calculations directly in MATLAB
        [midx1,midx2,midx3,midx4] = ind2sub(num_modes*ones(1,4),1:num_modes^4);

        precision = 8;

        num_segments = ceil(Nx^2*num_modes^4*precision/mem_planned_for_SRSK);
        num_each_segment = ceil(num_modes^4/num_segments);
        segments = [num_each_segment*ones(1,num_segments-1) num_modes^4-num_each_segment*(num_segments-1)];
        cs = [0 cumsum(segments)];

        for segi = 1:num_segments
            s = cs(segi)+1;
            e = cs(segi+1);
            SR(s:e) = permute(sum(sum(fields(:, :, midx1(s:e)).*fields(:, :, midx2(s:e)).*fields(:, :, midx3(s:e)).*fields(:, :, midx4(s:e)),1),2),[3 1 2])./...
                          (norms(midx1(s:e)).*norms(midx2(s:e)).*norms(midx3(s:e)).*norms(midx4(s:e)));
        end
        fprintf('Finished calculating SR for wavelength %s Ang\n',l_str);
    end
    % Give SR the correct dimensions
    if gpu_yes
        SR = gather(SR);
    end
    SR = SR/dx^2;

    %% Eliminate the zero elements
    thresholdzero = max(abs(SR))/1e5; % This is fairly arbitrary

    zero_SR = find(abs(SR)<thresholdzero);
    SR(zero_SR) = 0;
    cnt = num_modes^4 - length(zero_SR);

    fprintf('Calculated %d nonzero entries in the S_R tensor\n', cnt);

    %% Save to disk
    SR = reshape(SR, [num_modes, num_modes, num_modes, num_modes]);
    Aeff = 1./SR(1,1,1,1,:);
    Aeff(SR(1,1,1,1,:)==0) = NaN;

    save([dir_prefix '/S_tensors_' num2str(num_modes) 'modes'], 'SR', 'Aeff');

end

%% sub-function
function recompile_ptx(cudaFilename,ptxFilename)
    if ispc
        system(['nvcc -ptx ', cudaFilename, ' --output-file ', ptxFilename]);
    else % unix
        % tested: Debian 9
        % Cuda 8 doesn't support gcc6, beware to use gcc5 or clang-3.8.
        system(['nvcc -ccbin clang -ptx ', cudaFilename, ' --output-file ', ptxFilename]);
    end
end