%clearvars; close all;

% High-resolution spectrogram
log_yes = false; % Use "log_yes = true" to see, under log scale, how spectral interference generates temporal fringes.
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    % Spectrum
    figs = figure;
    plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,i)),1)).^2*c./lambda.^2,'linewidth',2);
    xlabel('Wavelength (nm)');
    ylabel('PSD (a.u.)');
    xlim([1500,2000]);
    ylim([0,45]);
    title('Spectrum');
    set(gca,'fontsize',20);
    l = legend('x','y'); set(l,'location','northwest');
    set(figs,'Color',[1,1,1]);

    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('vector_Raman');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);