__global__ void decompose_into_modes_Cartesian_basis(double2* field_wm,
                                                     const double* normalized_mode_space_profiles_xym, const double* norm_spatial_modes, const double2* full_field_xyw,
                                                     const double dx, const double dy,
                                                     const unsigned int sM, const unsigned int sX, const unsigned int sY, const unsigned int sW) {
    unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (threadIdx.x >= sM) return;
    if (blockIdx.x >= (sX*sY*sW)) return;
    if (thread_idx >= (sX*sY*sW*sM)) return;

    const unsigned int sXY = sX*sY;
    const unsigned int XYi = blockIdx.x % sXY;
    const unsigned int wi = blockIdx.x/sXY;

    if (XYi == 0) {
        for (unsigned int i = 0; i<sXY; i++) {
            field_wm[wi+sW*threadIdx.x].x = field_wm[wi+sW*threadIdx.x].x + normalized_mode_space_profiles_xym[i+sXY*threadIdx.x]*full_field_xyw[i+sXY*wi].x;
            field_wm[wi+sW*threadIdx.x].y = field_wm[wi+sW*threadIdx.x].y + normalized_mode_space_profiles_xym[i+sXY*threadIdx.x]*full_field_xyw[i+sXY*wi].y;
        }
        field_wm[wi+sW*threadIdx.x].x = field_wm[wi+sW*threadIdx.x].x*dx*dy*norm_spatial_modes[threadIdx.x];
        field_wm[wi+sW*threadIdx.x].y = field_wm[wi+sW*threadIdx.x].y*dx*dy*norm_spatial_modes[threadIdx.x];
    }
}