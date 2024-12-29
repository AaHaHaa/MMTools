function Frame = animator(Frame,...
                          A,...
                          z,MFD,start_idx,...
                          Nt,dt,Nx,dx,lambda)

% Make them column vectors
z = z(:);
MFD = MFD(:);

x = (-Nx/2:Nx/2-1)*dx*1e6; % spatial vector

c = 299792.458; % nm/ps
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

for j = 1:size(A,4)-1
    if exist('fig','var')
        figure(fig);
    else
        fig = figure;
    end

    subplot(2,2,1);
    A_spatial = A(Nt/2,:,:,j+1);
    pcolor(x,x,abs(squeeze(A_spatial)).^2);
    shading interp; colormap(jet);
    xlabel('x (\mum)');
    xlim([-160,160]);
    ylim([-160,160]);

    subplot(2,2,2);
    plot(z(1:start_idx+j)*1e2,MFD(1:start_idx+j)*1e3,'Color','k','linewidth',2);
    xlabel('z (cm)');
    ylabel('MFD (\mum)');
    ylim([80,160]);

    subplot(2,2,[3,4]);
    spectrum = abs(fftshift(ifft(A(:,Nx/2,Nx/2,j+1),[],1),1)).^2*factor_correct_unit.*factor; % in wavelength domain
    plot(lambda,spectrum,'Color','b','linewidth',2);
    xlabel('Wavelength (nm)');
    ylabel('PSD (nJ/nm/m^2)');
    xlim([950,1100]);
    ylim([0,max(spectrum)]);

    set(fig,'Color',[1,1,1]);

    Frame(j) = getframe(fig);
end
close(fig);

end