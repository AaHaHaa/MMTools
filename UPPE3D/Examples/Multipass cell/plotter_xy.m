function fig = plotter_xy(fig,...
                          A,...
                          z,MFD,...
                          Nt,Nx,dx,lambda)

% Make them column vectors
z = z(:);
MFD = MFD(:);

x = (-floor(Nx/2):floor((Nx-1)/2))*dx*1e6; % spatial vector

if isempty(fig)
    fig = figure;
else
    figure(fig); % reuse the previous figure handle; plot in the same figure
end
subplot(2,2,1)
input_A = A(Nt/2,:,:,1);
pcolor(x,x,abs(squeeze(input_A)).^2);
shading interp; colormap(jet);
xlabel('x (\mum)');
xlim([-600,600]);
ylim([-600,600]);
title('Beam profile input');

subplot(2,2,2)
output_A = A(Nt/2,:,:,end);
pcolor(x,x,abs(squeeze(output_A)).^2);
shading interp; colormap(jet);
xlabel('x (\mum)');
xlim([-600,600]);
ylim([-600,600]);
title('Beam profile output');

subplot(2,2,3)
plot(z*1e2,MFD*1e3,'Color','k','linewidth',2);
xlabel('z (cm)');
title('Beam size evolution')

subplot(2,2,4)
spectrum = abs(fftshift(ifft(A(:,:,:,end)),1)).^2./lambda.^2; % "./lambda" is to make it into the wavelength domain (not with the correct unit but the correct relative strength; we'll normalize it later)
% remove the weak spectral signal from spatial integration
multiplication_ratio = 3;
max_spectrum = max(spectrum(:));
spectrum = spectrum./max_spectrum; % make spectrum from 0-1
spectrum = spectrum.^multiplication_ratio*max_spectrum;
center_spectrum = spectrum(:,floor(Nx/2)+1,floor(Nx/2)+1);
avg_spectrum = sum(spectrum,[2,3]);
avg_spectrum = avg_spectrum/max(avg_spectrum); % normalized
center_spectrum = center_spectrum/max(center_spectrum); % normalized
plot(lambda,avg_spectrum,'Color','b','linewidth',2);
hold on;
plot(lambda,center_spectrum,'Color','r','linewidth',2);
hold off;
xlabel('Wavelength (nm)');
xlim([950,1100]);
legend('Avg spectrum','Center spectrum');
title('Spectral domain')

drawnow; % plot it during running the code

end