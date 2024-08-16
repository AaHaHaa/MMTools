function Ef_out = add_thin_lens_phase(Ef,dx,dy,wavelength,focal_length)
%ADD_LENS_PHASE It adds a thin-lens phase for a 3D electric field
%   Ef: electric field in the frequency domain(=ifft(Et)); size: (Nt,Nx,Ny,Nz)
%   dx: (m)
%   dy: (m)
%   wavelength: wavelength in the same order as Ef (m)
%   focal_length: focal length of a lens

[~,Nx,Ny,~] = size(Ef);

x = (-Nx/2:Nx/2-1)*dx;
y = permute((-Ny/2:Ny/2-1)*dy,[1,3,2]);

lens_phase = -pi./wavelength/focal_length.*(x.^2+y.^2);

Ef_out = Ef.*exp(1i*lens_phase);

end