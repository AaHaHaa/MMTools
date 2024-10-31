close all; clearvars;

num_modes = 21;
save_num = 101;
file_num = 20;

all_each_energy_ratio = zeros(save_num,num_modes,file_num);
all_total_peak_power = zeros(save_num,1,file_num);
for i = 1:file_num
    load(sprintf('MM_YDFA%u_10000ps.mat',i));
    
    each_energy = permute(trapz(abs(output_field.fields).^2)*dt/1e3,[3,2,1]); % nJ
    each_energy_ratio = each_energy./sum(each_energy,2);
    
    total_peak_power = squeeze(max(sum(abs(output_field.fields).^2,2),[],1)/1e6); % MW
    
    all_each_energy_ratio(:,:,i) = each_energy_ratio;
    all_total_peak_power(:,:,i) = total_peak_power;
end
mean_each_energy_ratio = mean(all_each_energy_ratio,3);
mean_total_peak_power = mean(all_total_peak_power,3);

z = output_field.z;

%% Plot results
figure;
plot(z,mean_each_energy_ratio(:,1),'linewidth',2,'Color','k');
hold on;
ccc = distinguishable_colors(num_modes);
for i = 2:num_modes
    plot(z,mean_each_energy_ratio(:,i),'linewidth',2,'Color',ccc(i-1,:));
end
plot(z,squeeze(all_each_energy_ratio(:,1,:)),'linewidth',2,'Color',[0.753,0.753,0.753]);
for i = 1:file_num
    plot(z,squeeze(all_each_energy_ratio(:,2:end,i)),'linewidth',2,'Color',[0.467,0.502,0.565]);
end
ylim([0,1]);
ylabel('Mode occupation');
set(gca,'fontsize',20,'Color','None');
set(gca,'XTick',[],'YTick',0:0.5:1);
line_order = get(gca,'Children'); set(gca,'Children',flipud(line_order));
fig_pos = get(gcf,'Position'); set(gcf,'Units','Pixels','Position',[round(fig_pos(1)*0.5),fig_pos(2)*0.5,fig_pos(3)*2,fig_pos(4)*0.8]);
ax_pos = get(gca,'Position'); set(gca,'Position',[0.08,ax_pos(2),ax_pos(3)*1.1,ax_pos(4)]); set(gca,'Units','Pixels');
print(gcf,'Energy ratio (10000ps).pdf','-dpdf','-bestfit');

figure;
plot(z,mean_total_peak_power,'linewidth',2,'Color','k');
hold on;
ccc = distinguishable_colors(num_modes-1);
for i = 1:file_num
    plot(z,all_total_peak_power(:,:,i),'linewidth',2,'Color',[0.753,0.753,0.753]);
end
xlabel('Propagation distance (m)');
ylabel('Peak power (MW)');
set(gca,'fontsize',20,'Color','None');
ylim([0,0.25]);
set(gca,'XTick',0:2:10,'YTick',0:0.1:0.2);
line_order = get(gca,'Children'); set(gca,'Children',flipud(line_order));
fig_pos = get(gcf,'Position'); set(gcf,'Units','Pixels','Position',[round(fig_pos(1)*0.5),fig_pos(2)*0.5,fig_pos(3)*2,fig_pos(4)]);
ax_pos = get(gca,'Position'); set(gca,'Position',[0.08,ax_pos(2),ax_pos(3)*1.1,ax_pos(4)]); set(gca,'Units','Pixels');
print(gcf,'Peak power (10000ps).pdf','-dpdf','-bestfit');

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

probe_z = [1,9,51,101];
for i = 1:length(probe_z)
    % Electric field
    full_field_txy = recompose_into_space(false, mode_profiles, output_field.fields(:,:,probe_z(i)), '');

    [xx,yy] = meshgrid(x,x);
    plot_x = [x-x(end),x-x(1)+mean(diff(x))];
    [plot_xx,plot_yy] = meshgrid(plot_x,plot_x);
    I = squeeze(sum(abs(full_field_txy).^2,1));
    plot_I = interp2(xx,yy,I,plot_xx,plot_yy,'spline');
    
    figure;
    pcolor(plot_x,plot_x,plot_I);
    ccc = whitejet_lower(256);
    shading interp; colormap(ccc);
    max_I = max(max( squeeze(sum(abs(full_field_txy).^2,1)) ));
    caxis([0,max_I]);
    pbaspect([1 1 1]);
    xlim([-45,45]);
    ylim([-45,45]);
    set(gca,'Position',[0,0,1,1]);
    fig_pos = get(gcf,'Position'); set(gcf,'Position',[fig_pos(1:2),fig_pos(4),fig_pos(4)]);
    set(gca,'XTick',[],'YTick',[],'Color','None','XColor','None','YColor','None');

    hold on; plot(unit_x,unit_y,'linewidth',5,'Color','k'); hold off;

    print(gcf,sprintf('E%u (10000ps).png',probe_z(i)),'-dpng');
    
    % Inversion profile
    figure;
    N2_x = gain_rate_eqn.mode_profile_dx*(-size(output_field.N2,1)/2:size(output_field.N2,1)/2-1);
    [N2_xx,N2_yy] = meshgrid(N2_x,N2_x);
    plot_x = [N2_x-N2_x(end),N2_x-N2_x(1)+mean(diff(N2_x))];
    [plot_xx,plot_yy] = meshgrid(plot_x,plot_x);
    plot_N2 = interp2(N2_xx,N2_yy,output_field.N2(:,:,probe_z(i)),plot_xx,plot_yy,'spline',0);
    pcolor(plot_xx,plot_yy,plot_N2);
    shading interp; colormap(jet);
    max_N2 = max(max(output_field.N2(:,:,probe_z(i))));
    caxis([0,max_N2]);
    pbaspect([1 1 1]);
    xlim([-45,45]);
    ylim([-45,45]);
    set(gca,'Position',[0,0,1,1]);
    fig_pos = get(gcf,'Position'); set(gcf,'Position',[fig_pos(1:2),fig_pos(4),fig_pos(4)]);
    set(gca,'XTick',[],'YTick',[],'Color','None','XColor','None','YColor','None');
    
    hold on; plot(unit_x,unit_y,'linewidth',5,'Color','k'); hold off;
    
    print(gcf,sprintf('N2%u (10000ps).png',probe_z(i)),'-dpng');
end

%% Calculate the absorption
all_absorption = zeros(101,34);
transition = 0;
for probe_z = 1:101
    N2_x = gain_rate_eqn.mode_profile_dx*(-size(output_field.N2,1)/2:size(output_field.N2,1)/2-1);
    N2_ratio = squeeze(output_field.N2(round(size(gain_rate_eqn.N_total,1)/2),:,probe_z));
    N_total = max(gain_rate_eqn.N_total(:));
    N2 = permute(N2_ratio*N_total,[1,3,2]);
    absorption = fftshift(permute((gain_rate_eqn.cross_sections.emission+gain_rate_eqn.cross_sections.absorption).*N2-gain_rate_eqn.cross_sections.absorption*N_total,[5,3,1,2,4]),1)*1e6;
    
    all_absorption(probe_z,:) = absorption(1,:);
    
    if absorption(1,17) < 0 && transition == 0
        transition = probe_z;
    end
end

figure;
pcolor(output_field.z,N2_x,-10*log10(exp(1))*all_absorption.');
shading interp; colormap(jet); cb = colorbar;
ylim([-40,40]);
xlim([output_field.z(1),output_field.z(end)]);
caxis([0,max(-10*log10(exp(1))*all_absorption(:))]);
hold on; plot(output_field.z(transition*ones(length(N2_x),1)),N2_x,'linewidth',4,'Color','w','LineStyle','--'); hold off;
set(gca,'fontsize',30);
set(gca,'YTick',[-40,0,40],'XTick',[]);
ylabel('x (µm)');
%xlabel('Propagation distance (m)');
fig_pos = get(gcf,'Position'); set(gcf,'Units','Pixels','Position',[round(fig_pos(1)*0.5),fig_pos(2)*0.5,fig_pos(3)*2.4,fig_pos(4)*0.7]);
ax_pos = get(gca,'Position'); set(gca,'Position',[0.09,ax_pos(2),ax_pos(3)*1.08,ax_pos(4)]); set(gca,'Units','Pixels');
cb.Label.String = 'Loss (1/m in dB)';
print(gcf,'loss_effect (10000ps).png','-dpng');