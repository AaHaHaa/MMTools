%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script builds the fiber index profile and calls the svmodes function
% to solve for the lowest "num_modes" modes over a range of frequencies.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('helpers/');

%% Set parameters (users modify only this block)
% Frequency window
Nf = 10; % number of frequency points at which the modes will be calculated; usually 20
wavelength0 = 1030e-9; % center wavelength, in m
freq_range = 100; % THz; frequency range, in m. If 0, only the center wavelength will be used. Usually 100 THz.
num_modes = 1; % number of modes to compute; you can use a large number, since this code can find the maximum supported modes itself
include_cladding_modes = false;

% Spatial profile
Nx = 400; % number of spatial grid points
spatial_window = 50; % full spatial window size, in um, usually set to 100 um

% Extra parameters:
% (1) step fiber:
% No extra parameters for the step-index fibers.
% (2) GRIN fiber:
alpha = 2; % Shape parameter

use_fiber_collection = true; % use fibers in "fiber_collections" function
if use_fiber_collection
    fiber = '1060XP';
    
    [fiber_type,core_diameter,clad_diameter,core_NA,clad_NA,fname_user_defined,alpha] = fiber_collections(fiber,wavelength0);
else
    fiber_type = 'GRIN'; % type of fiber
    core_diameter = 62.5; % core diameter of fiber, in um
    clad_diameter = 245; % cladding diameter of fiber, in um
    % Since I found out most commercial fibers show only NA, it's more convenient to use NA here.
    core_NA = 0.275;
    clad_NA = 0.22;

    fname_user_defined = 'GRIN-62.5_245_wavelength1030nm'; % the folder name
end

% Sellmeier coefficients
material = 'fused silica';
[a,b] = Sellmeier_coefficients(material);

%% Calculate the modes
core_radius = core_diameter/2; % core radius of fiber, in um
clad_radius = clad_diameter/2; % cladding radius, in um

if isempty(fname_user_defined)
    folder_name = sprintf(['Fibers/%s_wavelength%4unm_diameter' num2str(core_diameter) 'um'],fiber_type,wavelength0*1e9); % folder where the output will be stored
else
    folder_name = ['Fibers/' fname_user_defined];
end

if ~exist(folder_name,'dir')
    mkdir(folder_name);
end

if freq_range == 0
    Nf = 1;
end

% Set the range in frequency space, which is more objective
if freq_range == 0
    wavelength = wavelength0*1e6;
else
    c = 299.792458; % um/ps
    freq_min = c/wavelength0*1e-6 + freq_range/2;
    freq_max = c/wavelength0*1e-6 - freq_range/2;
    
    f = linspace(freq_min,freq_max,Nf)'; % THz
    wavelength = c./f; % um
end

% Calculate the index difference using the Sellmeier equation to generate n(lambda)
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

n_core = n_from_Sellmeier(wavelength);
n_clad = sqrt(n_core.^2-core_NA^2);
n_coat = sqrt(n_clad.^2-clad_NA^2);

% At each wavelength, calculate the modes
save_data(Nf) = struct('x',[],'phi',[],'epsilon',[],'neff',[]);
field = 'EX'; % See svmodes for details
boundary = '0000'; % See svmodes for details

%% Calculate at the max frequency first to determine the number of modes
wavelength_min = min(wavelength(wavelength>0));
n_core0 = n_from_Sellmeier(wavelength_min);
n_clad0 = sqrt(n_core0.^2-core_NA^2);
n_coat0 = sqrt(n_clad0.^2-clad_NA^2);
n0 = {n_core0,n_clad0,n_coat0};
[epsilon,...
 x, dx] = run_profile_function(fiber_type,...
                               Nx, spatial_window,...
                               core_radius, clad_radius,...
                               n0, alpha);
guess = sqrt(epsilon(Nx/2, Nx/2));

% Quickly show the index profile to make sure everything's working correctly
gg = figure;
subplot(2,1,1)
pcolor(x,x,sqrt(epsilon));
colormap(gray); colormap(flipud(colormap)); shading interp; axis square

subplot(2,1,2);
plot(x,sqrt(epsilon(:,Nx/2)));

saveas(gg,[folder_name '/fiberprofile'],'fig');
print(gg,[folder_name '/fiberprofile'],'-dpng');
close(gg)

% Actually do the calculation
t_justsolve = tic();
[phi1,neff1] = svmodes(wavelength_min,guess,num_modes,dx,dx,epsilon,boundary,field);
phi1 = phi1./sqrt(sum(sum(abs(phi1).^2,1),2))/dx;
toc(t_justsolve);

% Exclude lossy unconfined modes
threshold = core_diameter*2; % If the mode isn't confined within the core,
                             % it'll be solved as an eigenmode of the whole spatial window by the svmmodes code;
                             % therefore, its D4Sigma will be much larger than the core diameter
save_midx = 1:num_modes;
delete_midx = zeros(1,num_modes);
for midx = 1:num_modes
    mode_intensity = phi1(:,:,midx).^2;
    background_intensity = mean(mode_intensity(mode_intensity<max(mode_intensity(:))/1000));
    [D4SigmaX,D4SigmaY] = calcMFD(mode_intensity,background_intensity);
    D4Sigma = max(D4SigmaX,D4SigmaY);

    if include_cladding_modes
        [xx,yy] = meshgrid(x,x);
        core_region = sqrt(xx.^2+yy.^2)<core_radius;
        if max(max(mode_intensity(core_region)))/max(max(mode_intensity)) < 0.5
            delete_midx(midx) = midx;
        end
    else
        if D4Sigma*dx > threshold
            num_modes = midx-1;
            save_midx = 1:num_modes;
            fprintf('num_modes is reset to %u\n',num_modes);
            break
        end
    end
end
delete_midx(delete_midx==0) = [];
save_midx(delete_midx) = [];

% the fiber to draw on the mode profiles
theta = linspace(0,2*pi,1000);
core_unit_x = core_radius*cos(theta);
core_unit_y = core_radius*sin(theta);
clad_unit_x = clad_radius*cos(theta);
clad_unit_y = clad_radius*sin(theta);

% (1) Plot the mode profiles out
% (2) Figure out which modes are degenerate
degenerate_list = cell(num_modes,1);
list_i = 1; list_j = 1;
ccc = whitejet_centered(256);
for ii = 1:length(save_midx)
    midx = save_midx(ii);
    
    fname = [folder_name '/mode' num2str(midx) 'wavelength' num2str(round(wavelength0*1e10))];

    gg = figure('Position',[1 1 1200 800]);

    phi = phi1(:,:,midx);
    neff = neff1(midx);
    if abs(max(phi(:))) < abs(min(phi(:)))
        phi = -phi;
    end
    if ii>1
        isdegenerate = abs((neff-neff1(midx-1))/neff) < 1e-5;
        if isdegenerate
            list_j = list_j + 1;
            degenerate_list{list_i}(list_j) = midx;
            
            normphi = phi./sqrt(sum(sum(abs(phi).^2)))/dx; % 1/um
            normphi_p = phi1(:,:,midx-1)./sqrt(sum(sum(abs(phi1(:,:,midx-1)).^2)))/dx; % 1/um
            if sum(sum(abs(abs(normphi).^2-rot90(abs(normphi_p)).^2)))*dx^2 < 1e-2
                phi = rot90(phi1(:,:,midx-1));
            end
            neff = neff1(midx-1);
        else
            list_i = list_i + 1;
            list_j = 1;
            degenerate_list{list_i}(list_j) = midx;
        end
    else
        degenerate_list{1} = 1;
    end
    save(fname,'phi','neff','x','epsilon');

    normphi = phi./sqrt(sum(sum(abs(phi).^2)))/dx; % 1/um
    colormap_range = max(abs([max(normphi(:)),min(normphi(:))]));
    pcolor(x,x,normphi); shading interp; colormap(ccc); caxis([-colormap_range,colormap_range]); axis square;
    hold on;
    plot(core_unit_x,core_unit_y,'linewidth',2,'Color','k');
    plot(clad_unit_x,clad_unit_y,'linewidth',2,'Color','k'); hold off;
    title(['n_{eff} = ' num2str(neff)]);

    % Save the file with identifying information
    print(gg,fname,'-dpng');
    close(gg);
end
degenerate_list = degenerate_list(~cellfun(@isempty,degenerate_list));

%% Calculate the other frequencies
parfor kk = 1:Nf
    wavelength_kk = wavelength(kk); % wavelength
    
    % Build the index profile. The funcation can be arbitrary, and can take any extra parameters
    n_kk = {n_core(kk),n_clad(kk),n_coat(kk)};
    [epsilon,...
     x, dx] = run_profile_function(fiber_type,...
                                   Nx, spatial_window,...
                                   core_radius, clad_radius,...
                                   n_kk, alpha);
    guess = sqrt(epsilon(Nx/2, Nx/2));

    % Actually do the calculation
    t_justsolve = tic();
    [phi1,neff1] = svmodes(wavelength_kk,guess,num_modes,dx,dx,epsilon,boundary,field);
    toc(t_justsolve);
    
    for midx = save_midx
        if midx>1
            isdegenerate = abs((neff1(midx)-neff1(midx-1))/neff1(midx)) < 1e-5;
            if isdegenerate
                normphi = phi1(:,:,midx)./sqrt(sum(sum(abs(phi1(:,:,midx)).^2)))/dx; % 1/um
                normphi_p = phi1(:,:,midx-1)./sqrt(sum(sum(abs(phi1(:,:,midx-1)).^2)))/dx; % 1/um
                if sum(sum(abs(abs(normphi).^2-rot90(abs(normphi_p)).^2)))*dx^2 < 1e-2
                    phi1(:,:,midx) = rot90(phi1(:,:,midx-1));
                end
                neff1(midx) = neff1(midx-1);
            end
        end
    end
    
    for mode_group = 1:length(degenerate_list)
        neff1(degenerate_list{mode_group}) = mean(neff1(degenerate_list{mode_group}));
    end
    
    save_data(kk).x = x;
    save_data(kk).phi = phi1(:,:,save_midx)./sqrt(sum(sum(abs(phi1(:,:,save_midx)).^2,1),2))/dx; % 1/um
    save_data(kk).epsilon = epsilon;
    save_data(kk).neff = neff1(save_midx);
end

%% Save data
for kk = 1:Nf
    for ii = 1:length(save_midx)
        midx = save_midx(ii);
        
        wavelength_kk = wavelength(kk);
        fname = [folder_name '/mode' num2str(midx) 'wavelength' num2str(round(wavelength_kk*10000))];
        ki_save_data = save_data(kk);
        ki_save_data.phi = ki_save_data.phi(:,:,ii); % pick only the ii-th mode
        ki_save_data.neff = ki_save_data.neff(ii); % pick only the ii-th mode
        save(fname,'-struct','ki_save_data');
    end
end

% Save fiber information
fid = fopen(sprintf('%s/fiber_info.txt',folder_name),'w');
fprintf(fid,'Material: %s\nFiber type: %s\ncore NA: %6.5f\nclad NA: %6.5f\ncore diameter: %6.2fum\nclad diameter: %6.2fum\nboundary: %s\nfield: %s\n\nNx: %d\nspatial_window: %dum\n\n',material,fiber_type,core_NA,clad_NA,core_diameter,clad_diameter,boundary,field,Nx,spatial_window);
if isequal(fiber_type,'GRIN')
    fprintf(fid,'alpha: %6.5f\n',alpha);
end
fprintf(fid,'save_midx: ');
fprintf(fid,'%u ',save_midx);
fprintf(fid,'\n');
fclose(fid);

%% helper function
function [epsilon, x, dx] = run_profile_function(fiber_type, Nx, spatial_window, core_radius, clad_radius, n, alpha)

profile_function = str2func(sprintf('build_%s',fiber_type)); % function that builds the fiber
switch fiber_type
    case 'step'
        [epsilon, x, dx] = profile_function(Nx, spatial_window, core_radius, clad_radius, n);
    case 'GRIN'
        [epsilon, x, dx] = profile_function(Nx, spatial_window, core_radius, clad_radius, n, alpha);
end

end