% This code plots the result from "nonlinear_gain_coupling".

addpath('../../');
close all; clearvars;

load('nonlinear_gain_coupling2.mat');

%% Results
energy_linear = permute(trapz(abs(prop_output_linear.fields).^2,1)*dt/1e3,[3 2 1]);
energy_circular = permute(trapz(abs(prop_output_circular.fields).^2,1)*dt/1e3,[3 2 1]);
energy_elliptical = cell(9,1);
for i = 1:9
    energy_elliptical{i} = permute(trapz(abs(prop_output_elliptical(i).fields).^2,1)*dt/1e3,[3 2 1]);
end

[phi_linear,theta_linear] = calc_ellipticity( prop_output_linear.fields(:,:,end),0);
[phi_circular,theta_circular] = calc_ellipticity( prop_output_circular.fields(:,:,end),1);
phi_elliptical = cell(9,1); theta_elliptical = cell(9,1);
for i = 1:9
    [phi_elliptical{i},theta_elliptical{i}] = calc_ellipticity( prop_output_elliptical(i).fields(:,:,end),sim_e(i).ellipticity);
end

% Transform into the circular basis
P = @(r)[1  1i*r;...
         r   -1i]/sqrt(1+r^2); % the transformation matrix from basis (x,y) into one with ellipticity "r"
P_circ = @(r) P(1)/P(r); % finally transform into the circular basis

prop_output_linear_circ = prop_output_linear;
prop_output_elliptical_circ = prop_output_elliptical;
prop_output_linear_circ.fields = gather(permute(pagefun(@mtimes, P(1),gpuArray(permute(prop_output_linear.fields,[2 1 3]))),[2 1 3]));
for i = 1:9
    prop_output_elliptical_circ(i).fields = gather(permute(pagefun(@mtimes, P_circ(sim_e(i).ellipticity),gpuArray(permute(prop_output_elliptical(i).fields,[2 1 3]))),[2 1 3]));
end

% Under circular basis
energy_linear_circ = permute(trapz(abs(prop_output_linear_circ.fields).^2,1)*dt/1e3,[3 2 1]);
ratio_linear_circ = energy_linear_circ(:,1)./energy_linear_circ(:,2);
energy_elliptical_circ = cell(9,1);
ratio_elliptical_circ = cell(9,1);
for i = 1:9
    energy_elliptical_circ{i} = permute(trapz(abs(prop_output_elliptical_circ(i).fields).^2,1)*dt/1e3,[3 2 1]);
    ratio_elliptical_circ{i} = energy_elliptical_circ{i}(:,1)./energy_elliptical_circ{i}(:,2);
end
ratio_circular = energy_circular(:,1)./energy_circular(:,2);

% total_field
total_field_linear = sum(abs(prop_output_linear.fields(:,:,end)).^2,2);
total_field_circular = sum(abs(prop_output_circular.fields(:,:,end)).^2,2);
total_field_elliptical = cell(9,1);
for i = 1:9
    total_field_elliptical{i} = sum(abs(prop_output_elliptical(i).fields(:,:,end)).^2,2);
end

%% Plot Results

c = 299792458; % m/s

% linear polarization
fig = figure('Name','linear');
fp = get(fig,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(fig,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
subplot(4,3,1);
h = plot(prop_output_linear.z,energy_linear(:,1));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('x polarization');
subplot(4,3,4);
h = plot(prop_output_linear.z,energy_linear(:,2));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('y polarization');
subplot(4,3,2);
h = plot(prop_output_linear.z,energy_linear_circ(:,1));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_+ polarization');
subplot(4,3,5);
h = plot(prop_output_linear.z,energy_linear_circ(:,2));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_- polarization');
subplot(4,3,[3 6]);
h = plot(prop_output_linear.z,ratio_linear_circ);
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Ratio');
title('Energy ratio (\sigma_+/\sigma_-)');
subplot(4,3,7:9);
yyaxis left
plot(t,abs(prop_output_linear_circ.fields(:,:,end)).^2);
hold on; plot(t,total_field_linear,'r'); hold off;
ylabel('Intensity (W)');
legend('\sigma_+','\sigma_-','total field');
yyaxis right
h = plot(t,[phi_linear 270*ones(Nt,1) 90*ones(Nt,1)]);
set(h(1),'linewidth',2);
xlim([-10 0]);
xlabel('time (ps)'); ylabel('the phase difference (deg)');
title('Final fields');
subplot(4,3,10:12);
spectrum_linear = abs(fftshift(ifft(prop_output_linear_circ.fields(:,:,end)),1)).^2;
h = plot(lambda(lambda>0),spectrum_linear(lambda>0,:)*(Nt*dt)^2*c/1e6./lambda(lambda>0).^2);
set(h,'linewidth',2);
xlim([950,1200]);
ylabel('Intensity (nJ/nm)');
xlabel('Wavelength (nm)');
title('Spectrum');

% circular polarization
fig = figure('Name','circular');
fp = get(fig,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(fig,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
subplot(4,2,1);
h = plot(prop_output_circular.z,energy_circular(:,1));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_+ polarization');
subplot(4,2,3);
h = plot(prop_output_circular.z,energy_circular(:,2));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_- polarization');
subplot(4,2,[2 4]);
h = plot(prop_output_linear.z,ratio_circular);
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Ratio');
title('Energy ratio (\sigma_+/\sigma_-)');
subplot(4,2,5:6);
yyaxis left
plot(t,abs(prop_output_circular.fields(:,:,end)).^2);
hold on; plot(t,total_field_circular,'r'); hold off;
ylabel('Intensity (W)');
legend('\sigma_+','\sigma_-','total field');
yyaxis right
h = plot(t,[phi_circular 270*ones(Nt,1) 90*ones(Nt,1)]);
set(h(1),'linewidth',2);
xlim([-10 0]);
xlabel('time (ps)'); ylabel('the phase difference (deg)');
title('Final fields');
subplot(4,2,7:8);
spectrum_circular = abs(fftshift(ifft(prop_output_circular.fields(:,:,end)),1)).^2;
h = plot(lambda(lambda>0),spectrum_circular(lambda>0,:)*(Nt*dt)^2*c/1e6./lambda(lambda>0).^2);
set(h,'linewidth',2);
xlim([950,1200]);
ylabel('Intensity (nJ/nm)');
xlabel('Wavelength (nm)');
title('Spectrum');

% elliptical polarization
for i = 1:9
    fig = figure('Name',sprintf('ellipticity: %2.1f',sim_e(i).ellipticity));
    fp = get(fig,'position');
    screen_size = get(0,'ScreenSize');
    original_top = screen_size(4)-fp(2)-fp(4);
    set(fig,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
    subplot(4,3,1);
    h = plot(prop_output_linear.z,energy_elliptical{i}(:,1));
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
    title('the 1st polarization');
    subplot(4,3,4);
    h = plot(prop_output_linear.z,energy_elliptical{i}(:,2));
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
    title('the 2nd polarization');
    subplot(4,3,2);
    h = plot(prop_output_linear.z,energy_elliptical_circ{i}(:,1));
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
    title('\sigma_+ polarization');
    subplot(4,3,5);
    h = plot(prop_output_linear.z,energy_elliptical_circ{i}(:,2));
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
    title('\sigma_- polarization');
    subplot(4,3,[3 6]);
    h = plot(prop_output_linear.z,ratio_elliptical_circ{i});
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Ratio');
    title('Energy ratio (\sigma_+/\sigma_-)');
    subplot(4,3,7:9);
    yyaxis left
    plot(t,abs(prop_output_elliptical_circ(i).fields(:,:,end)).^2);
    hold on; plot(t,total_field_elliptical{i},'r'); hold off;
    ylabel('Intensity (W)');
    legend('\sigma_+','\sigma_-','total field');
    yyaxis right
    h = plot(t,[phi_elliptical{i} 270*ones(Nt,1) 90*ones(Nt,1)]);
    set(h(1),'linewidth',2);
    xlim([-10 0]);
    xlabel('time (ps)'); ylabel('the phase difference (deg)');
    title('Final fields');
    subplot(4,3,10:12);
    spectrum_elliptical_circ = abs(fftshift(ifft(prop_output_elliptical_circ(i).fields(:,:,end)),1)).^2;
    h = plot(lambda(lambda>0),spectrum_elliptical_circ(lambda>0,:)*(Nt*dt)^2*c/1e6./lambda(lambda>0).^2);
    set(h,'linewidth',2);
    xlim([950,1200]);
    ylabel('Intensity (nJ/nm)');
    xlabel('Wavelength (nm)');
    title('Spectrum');
end