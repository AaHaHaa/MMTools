% I copied this from my paper's code.
% During writing my multimode gain paper, I used the notation N2 for
% upper-state population, which is now changed into N1 in the new code.
% The population is stored as prop_output.population in the new code,
% rather than prop_output.N2 in the old code.

close all; clearvars;

load(sprintf('MM_YDFA20_10000ps.mat'));

%% Profile
MM_folder = '../../Fibers/GRIN_168_400_wavelength1030nm/';
mode_profiles = zeros(400, 400, length(sim.midx));
for ni = 1:length(sim.midx)
    n = sim.midx(ni);
    load(sprintf('%smode%uwavelength%u.mat',MM_folder,n,round(sim.lambda0*1e10)),'phi','x');
    mode_profiles(:,:,ni) = phi;
end
mode_profiles = mode_profiles./sqrt(sum(sum(abs(mode_profiles).^2,1),2));

% the fiber core to draw on the mode profiles
radius = gain_rate_eqn.core_diameter/2;
theta = linspace(0,2*pi,1000);
unit_x = radius*cos(theta);
unit_y = radius*sin(theta);

save_point = size(output_field.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    % Electric field
    full_field_txy = recompose_into_space(false, mode_profiles, output_field.fields(:,:,i), '');

    [xx,yy] = meshgrid(x,x);
    plot_x = [x-x(end),x-x(1)+mean(diff(x))];
    [plot_xx,plot_yy] = meshgrid(plot_x,plot_x);
    I = squeeze(sum(abs(full_field_txy).^2,1));
    plot_I = interp2(xx,yy,I,plot_xx,plot_yy,'spline');
    
    figs = figure;
    subplot(1,2,1);
    pcolor(plot_x,plot_x,plot_I);
    ccc = whitejet_lower(256);
    shading interp; colormap(ccc);
    max_I = max(max( squeeze(sum(abs(full_field_txy).^2,1)) ));
    caxis([0,max_I]);
    pbaspect([1 1 1]);
    xlim([-45,45]);
    ylim([-45,45]);
    title('I');
    set(gca,'fontsize',20);
    %set(gca,'Position',[0,0,1,1]);
    set(gca,'XTick',[],'YTick',[],'Color','None','XColor','None','YColor','None');

    hold on; plot(unit_x,unit_y,'linewidth',5,'Color','k'); hold off;
    
    % Inversion profile
    subplot(1,2,2);
    N2_x = gain_rate_eqn.mode_profile_dx*(-size(output_field.N2,1)/2:size(output_field.N2,1)/2-1);
    [N2_xx,N2_yy] = meshgrid(N2_x,N2_x);
    plot_x = [N2_x-N2_x(end),N2_x-N2_x(1)+mean(diff(N2_x))];
    [plot_xx,plot_yy] = meshgrid(plot_x,plot_x);
    plot_N2 = interp2(N2_xx,N2_yy,output_field.N2(:,:,i),plot_xx,plot_yy,'spline',0);
    pcolor(plot_xx,plot_yy,plot_N2);
    shading interp; colormap(jet);
    max_N2 = max(max(output_field.N2(:,:,i)));
    caxis([0,max_N2]);
    pbaspect([1 1 1]);
    xlim([-45,45]);
    ylim([-45,45]);
    title('N_1');
    set(gca,'fontsize',20);
    %set(gca,'Position',[0,0,1,1]);
    set(gca,'XTick',[],'YTick',[],'Color','None','XColor','None','YColor','None');
    
    hold on; plot(unit_x,unit_y,'linewidth',5,'Color','k'); hold off;

    set(figs,'Color',[1,1,1]);
    
    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('Field_N1');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);