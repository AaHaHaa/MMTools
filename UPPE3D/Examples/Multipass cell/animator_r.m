function Frame = animator_r(Frame,...
                            A,l0,...
                            z,MFD,start_idx,...
                            Nt,dt,r,lambda,...
                            plate_z)

% Make them column vectors
z = z(:);
MFD = MFD(:);

r = r*1e6; % um

% The field maximum should be at the spatial center, which requires special
% calculation for the FHATHA, Hankel transform with high accuracy.
A0 = Hankel_f_at_0(A,l0);

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
    plot(r,abs(A_spatial.').^2/max(abs(A0(Nt/2,:,j+1).').^2),'linewidth',2,'Color','b');
    xlabel('x (\mum)');
    xlim([0,600]);

    subplot(2,2,2);
    plot_plate_MFD = [70;1000];
    plot(plate_z(1)*1e2*[1;1],plot_plate_MFD,'Color','m','linewidth',0.3);
    plot(plate_z(2)*1e2*[1;1],plot_plate_MFD,'Color','b','linewidth',1.5);
    hold on;
    for i = 3:length(plate_z)
        switch mod(i-2,6)
            case {1,2,4,5}
                plot(plate_z(i)*1e2*[1;1],plot_plate_MFD,'Color','m','linewidth',0.5);
            case {0,3}
                plot(plate_z(i)*1e2*[1;1],plot_plate_MFD,'Color','b','linewidth',1.5);
        end
        hold on;
    end
    plot(z(1:start_idx+j)*1e2,MFD(1:start_idx+j)*1e3,'Color','k','linewidth',2);
    hold off;
    xlabel('z (cm)');
    ylabel('MFD (\mum)');
    xlim([0,plate_z(end)*1e2]); % cm
    ylim(plot_plate_MFD);

    subplot(2,2,[3,4]);
    spectrum = sum(abs(fftshift(ifft(A(:,1,j+1),[],1),1)).^2,[2,3])*factor_correct_unit.*factor; % in wavelength domain
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