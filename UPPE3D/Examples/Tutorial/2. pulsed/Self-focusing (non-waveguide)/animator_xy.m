function Frame = animator_xy(A,...
                             z,MFD,...
                             Nt,dt,Nx,dx,lambda,...
                             L0)

% Make them column vectors
z = z(:);
MFD = MFD(:);

x = (-floor(Nx/2):floor((Nx-1)/2))*dx*1e6; % um; spatial vector

c = 299792.458; % nm/ps
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

Frame(size(A,4)) = struct('cdata',[],'colormap',[]);
for j = 1:size(A,4)
    if exist('fig','var')
        figure(fig);
    else
        fig = figure;
    end

    subplot(2,2,1);
    A_spatial = A(Nt/2,:,:,j);
    pcolor(x,x,abs(squeeze(A_spatial)).^2);
    shading interp; colormap(jet);
    xlabel('x (\mum)');
    xlim([-1,1]*1e2);
    ylim([-1,1]*1e2);

    subplot(2,2,2);
    plot(z(1:j)*1e2,MFD(1:j),'Color','k','linewidth',2);
    xlabel('z (cm)');
    ylabel('MFD (\mum)');
    xlim([0,L0*1e2]); % cm
    ylim([40,105]); % um

    subplot(2,2,[3,4]);
    spectrum = abs(fftshift(ifft(A(:,:,:,j),[],1),1)).^2*factor_correct_unit.*factor; % nJ/nm/m^2
    avg_spectrum = sum(spectrum,[2,3])*dx^2; % nJ/nm
    plot(lambda,avg_spectrum,'Color','b','linewidth',2);
    hold on;
    plot(lambda,spectrum(:,floor(Nx/2)+1,floor(Nx/2)+1)/max(spectrum(:,floor(Nx/2)+1,floor(Nx/2)+1))*max(avg_spectrum),'Color','r','linewidth',2);
    hold off;
    xlabel('Wavelength (nm)');
    ylabel('PSD (nJ/nm)');
    xlim([1020,1040]);
    ylim([0,max(avg_spectrum)]);
    legend('Avg spectrum','Center spectrum (norm.)');

    set(fig,'Color',[1,1,1]);

    Frame(j) = getframe(fig);
end
close(fig);

end