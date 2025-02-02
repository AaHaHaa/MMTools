function Frame = animator_r(Frame,...
                            A,...
                            z,MFD,start_idx,...
                            Nt,dt,r,lambda,...
                            plate_z)

% Make them column vectors
z = z(:);
MFD = MFD(:);

r = r*1e6; % um

c = 299792.458; % nm/ps
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

for j = 1:size(A,3)-1
    if exist('fig','var')
        figure(fig);
    else
        fig = figure;
    end

    subplot(2,2,1);
    A_spatial = A(Nt/2,:,j+1);
    plot(r,abs(A_spatial.').^2/max(abs(A_spatial).^2),'linewidth',2,'Color','b');
    xlabel('x (\mum)');
    xlim([0,220]);

    subplot(2,2,2);
    plot_plate_MFD = [120;185];
    for i = 1:length(plate_z)-1
        plot(plate_z(i)*1e2*[1;1],plot_plate_MFD,'Color','r','linewidth',1);
        hold on;
    end
    plot(z(1:start_idx+j)*1e2,MFD(1:start_idx+j)*1e3,'Color','k','linewidth',2);
    hold off;
    xlabel('z (cm)');
    ylabel('MFD (\mum)');
    xlim([0,plate_z(end)*1e2]); % cm
    ylim(plot_plate_MFD);

    subplot(2,2,[3,4]);
    spectrum = abs(fftshift(ifft(A(:,:,j+1),[],1),1)).^2*factor_correct_unit.*factor; % nJ/nm/m^2
    % remove the weak spectral signal from spatial integration
    multiplication_ratio = 3;
    max_spectrum = max(spectrum(:));
    spectrum = spectrum./max_spectrum; % make spectrum from 0-1
    spectrum = spectrum.^multiplication_ratio*max_spectrum;
    avg_spectrum = 2*pi*trapz(r,spectrum.*r,2); % nJ/nm
    plot(lambda,avg_spectrum,'Color','b','linewidth',2);
    hold on;
    plot(lambda,spectrum(:,1)/max(spectrum(:,1))*max(avg_spectrum),'Color','r','linewidth',2);
    hold off;
    xlabel('Wavelength (nm)');
    ylabel('PSD (nJ/nm)');
    xlim([950,1100]);
    ylim([0,max(avg_spectrum)]);
    legend('Avg spectrum','Center spectrum (norm.)');

    set(fig,'Color',[1,1,1]);

    Frame(j) = getframe(fig);
end
close(fig);

end