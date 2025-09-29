close all; clearvars;

load('SPM.mat');

%% Plot
% Time
figure;
h = plot(t,abs(prop_output.fields(:,:,end)).^2);
xlim([-3,3]);
xlabel('t');
ylabel('Power');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
figure;
h = plot(f-sim.f0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2);
xlim([-20,20]);
xlabel('\nu-\nu_0');
ylabel('PSD');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-3,3]);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x,y] = meshgrid(f-sim.f0,prop_output.z(2:end));
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-20,20]);
shading interp; colormap(jet);
xlabel('\nu-\nu_0');
ylabel('z/LD');
title('Spectrum during propagation');
set(gca,'fontsize',14);

% Movie (spectrogram)
pause(1);

% High-resolution spectrogram
log_yes = false; % Use "log_yes = true" to see, under log scale, how spectral interference generates temporal fringes.
save_point = size(prop_output.fields,3);
Frame1(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    [~,~,~,figs,ax] = calc_spectrogram(t,f,prop_output.fields(:,1,i),true,[-3,3],[940,1120],400,400,true,true,log_yes);
    set(figs,'Color',[1,1,1]);

    Frame1(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame1,20);

exportVideo = VideoWriter('SPM_spectrogram');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame1);
close(exportVideo);

% Movie (spectrum)
pause(1);

spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2;
spectrum = spectrum/max(spectrum(:));
Frame2(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    
    figs = figure;
    h = plot(f-sim.f0,spectrum(:,:,i),'linewidth',2,'Color','b');
    xlim([-20,20]);
    xlabel('\nu-\nu_0');
    ylabel('PSD');
    title('Spectrum');
    set(gca,'fontsize',20);
    set(figs,'Color',[1,1,1]);

    Frame2(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame2,20);

exportVideo = VideoWriter('SPM_spectrum');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame2);
close(exportVideo);