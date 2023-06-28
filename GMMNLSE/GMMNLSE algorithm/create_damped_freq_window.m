function damped_freq_window = create_damped_freq_window(Nt)
% This function creates a sharp damped frequency window to remove anything
% around the edges, especially the high-frequency edge of the window.
% I use a super-Gaussian function here.

f = fftshift((1:Nt)',1);
fc = floor(Nt/2)+1;
ffwhm = Nt*0.9;
f0 = ffwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665
gexpo = 2*20; % 10 is to make it a sharp window

% A much sharper damped window is used to for the low-frequency side;
% otherwise, it'll easily remove the long-wavelength signal we want to see.
damped_freq_window = exp(-(f-fc).^gexpo/(2*f0^gexpo)).^20; % 20 is to make it a sharp window
damped_freq_window(damped_freq_window>0.99) = 1;
%damped_freq_window_low = exp(-(f-fc).^(gexpo*2)/(2*f0^(gexpo*2))).^4; % (gexpo*2) and 4 are to make it a much sharper window
%damped_freq_window(fc:end) = damped_freq_window_low(fc:end);
%damped_freq_window(fc:end) = 1;

end