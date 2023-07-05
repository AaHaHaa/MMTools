%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script deletes the mode profiles other than at the center
% wavelength to save the disk space.
% Please finish running calc_dispersion and calc_SRSK_tensors in case you
% need them again.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delete_modeprofiles(modes_used,folder_name,Nf,lambda0,freq_range)

    % File name parameters:
    dir_prefix = ['Fibers/', folder_name]; % folder containing the calculated modes

    %%
    num_modes = length(modes_used); % number of modes for which the tensors should be calculated

    if freq_range == 0
        Nf = 1;
    end

    % Set the range in frequency space, which is more objective
    c = 2.99792458e-4; % speed of ligth m/ps
    freq1 = c/lambda0 - freq_range/2;
    freq2 = c/lambda0 + freq_range/2;

    f = linspace(freq1,freq2,Nf)'; % THz
    lambda = c./f*1e6; % um

    for kk = 1:Nf
        if lambda(kk) ~= lambda0*1e9
            for ii = 1:num_modes
                l_str = num2str(round(lambda(kk)*10000));
                name = [dir_prefix '/mode',int2str(modes_used(ii)),'wavelength', l_str, '.mat'];
                delete(name);

                disp(['Deleted mode ', int2str(modes_used(ii)), ' under ', l_str ' Ang'])
            end
        end
    end
end