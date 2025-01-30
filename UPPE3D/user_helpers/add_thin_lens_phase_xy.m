function Ef_out = add_thin_lens_phase_xy(Ef,dx,dy,wavelength,focal_length)
%ADD_LENS_PHASE_XY It adds a thin-lens phase for a 3D electric field
%   Ef: electric field in the frequency domain(=ifft(Et)); size: (Nt,Nx,Ny,Nz)
%   dx: (m)
%   dy: (m)
%   wavelength: wavelength in the same order as Ef (m)
%   focal_length: focal length of a lens

[~,Nx,Ny,~] = size(Ef);

x = (-floor(Nx/2):floor((Nx-1)/2))*dx;
y = permute((-floor(Ny/2):floor((Ny-1)/2))*dy,[1,3,2]);

lens_phase = -pi./wavelength/focal_length.*(x.^2+y.^2);

Ef_out = Ef.*exp(1i*lens_phase);

end