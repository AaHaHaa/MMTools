function damped_window = create_damped_window_xy(Nt,Nx,Ny)
% This function creates sharp damped spectral and spatial windows to remove
% anything around the edges, especially the high-frequency edge of the 
% window.
% I use a super-Gaussian function here.

gexpo = 2*20; % 20 is to make it a sharp window

x = 1:Nx;
xc = floor(Nx/2)+1;
xfwhm = Nx*0.82;
x0 = xfwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665

y = permute(1:Ny,[1,3,2]);
yc = floor(Ny/2)+1;
yfwhm = Ny*0.82;
y0 = yfwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665

f = ifftshift((1:Nt)',1);
fc = floor(Nt/2)+1;
ffwhm = Nt*0.82;
f0 = ffwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665

% A much sharper damped window is used only for the low-frequency side;
% otherwise, it'll easily remove the long-wavelength signal we want to see.
damped_x_window = exp(-(x-xc).^gexpo/(2*x0^gexpo)).^20; % 20 is to make it a sharp window
damped_x_window(damped_x_window>0.99) = 1;

damped_y_window = exp(-(y-yc).^gexpo/(2*y0^gexpo)).^20; % 20 is to make it a sharp window
damped_y_window(damped_y_window>0.99) = 1;

damped_freq_window = exp(-(f-fc).^gexpo/(2*f0^gexpo)).^20; % 20 is to make it a sharp window
damped_freq_window(damped_freq_window>0.99) = 1;
damped_freq_window(fc:end) = 1;

% Multiply them together
damped_window = damped_freq_window.*damped_x_window.*damped_y_window;

end