__global__ void recompose_into_space(double2* full_field_xyw,
                                     const double* normalized_mode_space_profiles_xym, const double2* field_wm, const double* norm_spatial_modes,
                                     const unsigned int sM, const unsigned int sX, const unsigned int sY, const unsigned int sW) {
    unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (thread_idx >= (sW*sX*sY)) return;

    const unsigned int wi = thread_idx % sW;
    const unsigned int XYi = thread_idx / sW;
    const unsigned int sXY = sX*sY;

    for (unsigned int mi = 0; mi<sM; mi++) {
        full_field_xyw[thread_idx].x = full_field_xyw[thread_idx].x + normalized_mode_space_profiles_xym[XYi+sXY*mi]*field_wm[wi+sW*mi].x/norm_spatial_modes[mi];
        full_field_xyw[thread_idx].y = full_field_xyw[thread_idx].y + normalized_mode_space_profiles_xym[XYi+sXY*mi]*field_wm[wi+sW*mi].y/norm_spatial_modes[mi];
    }

}