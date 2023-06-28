%clearvars; close all;

addpath('../../../MMTools/GMMNLSE/user_helpers/')

MM_folder = '../Fibers/GRIN-YDF-30_400_wavelength1030nm/';
lambda0 = 1030e-9;
core_diameter = 30;
num_modes = 10;
mode_profiles = zeros(400, 400, num_modes);
for ni = 1:num_modes
    n = ni;
    load(sprintf('%smode%uwavelength%u.mat',MM_folder,n,round(lambda0*1e10)),'phi','x');
    mode_profiles(:,:,ni) = phi;
end
mode_profiles = mode_profiles./sqrt(sum(sum(abs(mode_profiles).^2,1),2));

% the fiber core to draw on the mode profiles
radius = core_diameter/2;
theta = linspace(0,2*pi,1000);
unit_x = radius*cos(theta);
unit_y = radius*sin(theta);

%fraction = sqrt([0.1,ones(1,num_modes-1)*0.9/(num_modes-1)]);%.*exp(1i*2*pi*rand(1,num_modes));
E = output_field.fields(:,:,50);
fraction2 = sqrt(sum(abs(E).^2,1)).*exp(1i*angle(E(Nt/2,:)));
full_field_txy = recompose_into_space(false, mode_profiles, fraction2, '');

figs = figure;
pcolor(x,x,squeeze(sum(abs(full_field_txy).^2,1)));
ccc = whitejet_lower(256);
shading interp; colormap(ccc);
set(figs,'Color',[1,1,1]);
hold on; plot(unit_x,unit_y,'linewidth',2,'Color','k'); hold off;