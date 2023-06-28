function mode_freq_profiles_wm = decompose_into_modes_cylindrical_symmetry(normalized_mode_space_profiles_wmxyp, full_field_wxyp, r, dr, sim)
%DECOMPOSE_INTO_MODES Get the modal decomposition from the full 3D field
% mode_space_profile - a (Nx, Nx, num_modes) matrix with the mode profile 
% for each mode.
%
% =========================================================================
% Input:
% normalized_mode_space_profiles_xym - a (Nr, 1, num_modes) matrix with each mode's profile in space. The units do not matter.
% dx - spatial grid spacing, in m
% full_field - a (Nt,Nr,1) matrix with full spatiotemporal fields in each time.
% -------------------------------------------------------------------------
% Output:
% mode_time_profiles - a (Nt, num_modes) matrix with each mode's time profile.
% =========================================================================

normalized_mode_space_profiles_xymwp = permute(normalized_mode_space_profiles_wmxyp,[3,4,2,1,5]);
full_field_xywp = permute(full_field_wxyp,[2,3,1,4]);

num_modes = size(normalized_mode_space_profiles_wmxyp, 2);
Nr = size(normalized_mode_space_profiles_wmxyp, 3);
Ngrid = Nr;
Nw = size(full_field_xywp, 3);
Nw_profiles = size(normalized_mode_space_profiles_wmxyp,1);
if sim.scalar
    Np = 1;
else
    Np = 2;
end

% Einstein summation convention: wm=(wxy)*(mxy)
mode_freq_profiles_wm =  2*pi*permute(sum(pagefun(@mtimes,reshape(full_field_xywp,[1,Ngrid,Nw,Np]),reshape(conj(normalized_mode_space_profiles_xymwp),[Ngrid,num_modes,Nw_profiles,Np]).*r*dr),4),[3,2,1]).*sim.norm_spatial_modes;

end