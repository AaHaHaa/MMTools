function fig = plotter_xy(fig,...
                          A,...
                          z,MFD,...
                          Nt,Nx,dx,lambda)

% Make them column vectors
z = z(:);
MFD = MFD(:);

x = (-Nx/2:Nx/2-1)*dx*1e6; % spatial vector

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
output_center_spectrum = abs(fftshift(ifft(A(:,Nx/2,Nx/2,end)),1)).^2./lambda.^2; % "./lambda" is to make it into the wavelength domain (not with the correct unit but the correct relative strength; we'll normalize it later)
output_center_spectrum = output_center_spectrum/max(output_center_spectrum); % normalized
plot(lambda,output_center_spectrum,'Color','b','linewidth',2);
xlabel('Wavelength (nm)');
xlim([950,1100]);
title('Spectral domain')

drawnow; % plot it during running the code

end