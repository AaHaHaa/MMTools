__global__ void decompose_into_modes_Cartesian_basis(double2* field_wm,
                                                     const double* mode_space_profiles_xym, const double* norm_spatial_modes,
                                                     const double2* full_field_wxy,
                                                     const double dx, const double dy,
                                                     const unsigned int sM, const unsigned int sX, const unsigned int sY, const unsigned int sW) {
    unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (thread_idx >= (sW*sM)) return;

    const unsigned int sXY = sX*sY;
    const unsigned int mi = thread_idx % sM;
    const unsigned int wi = thread_idx / sM;

    for (unsigned int XYi = 0; XYi<sXY; XYi++) {
        field_wm[wi+sW*mi].x = field_wm[wi+sW*mi].x + mode_space_profiles_xym[XYi+sXY*mi]*full_field_wxy[wi+sW*XYi].x;
        field_wm[wi+sW*mi].y = field_wm[wi+sW*mi].y + mode_space_profiles_xym[XYi+sXY*mi]*full_field_wxy[wi+sW*XYi].y;
    }
    field_wm[wi+sW*mi].x = field_wm[wi+sW*mi].x*dx*dy*norm_spatial_modes[mi];
    field_wm[wi+sW*mi].y = field_wm[wi+sW*mi].y*dx*dy*norm_spatial_modes[mi];
}