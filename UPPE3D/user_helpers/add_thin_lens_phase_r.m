function Ef_out = add_thin_lens_phase_r(Ef,r,wavelength,focal_length)
%ADD_LENS_PHASE_R It adds a thin-lens phase for a 3D electric field that is
% radially symmetric
%   Ef: electric field in the frequency domain(=ifft(Et)); size: (Nt,Nx,Ny,Nz)
%   r: (1,Nr); radial sampling positions (m)
%   wavelength: wavelength in the same order as Ef (m)
%   focal_length: focal length of a lens

lens_phase = -pi./wavelength/focal_length.*r.^2;

Ef_out = Ef.*exp(1i*lens_phase);

end