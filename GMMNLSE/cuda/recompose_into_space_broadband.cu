__global__ void recompose_into_space_broadband(double2* full_field_xypw,
                                               const double* normalized_mode_space_profiles_xypwm, const double2* field_wm, const double* norm_spatial_modes,
                                               const unsigned int sSPw, const unsigned int sM, const unsigned int sX, const unsigned int sY, const unsigned int sW, const unsigned int sP) {
    unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (threadIdx.x >= (sX*sY*sP)) return;
    if (blockIdx.x >= sW) return;
    if (thread_idx >= (sW*sX*sY*sP)) return;

    const unsigned int N = (sSPw*sX*sY*sP);

    const unsigned int max_num_modes = 20;
    __shared__ double2 field_wi[max_num_modes];
    __shared__ double norm_wi[max_num_modes];
    if (threadIdx.x < sM) { // Required: (sX*sY*sP) < sM
        field_wi[threadIdx.x].x = field_wm[blockIdx.x+sW*threadIdx.x].x;
        field_wi[threadIdx.x].y = field_wm[blockIdx.x+sW*threadIdx.x].y;

        norm_wi[threadIdx.x] = norm_spatial_modes[blockIdx.x+sW*threadIdx.x];
    }
    __syncthreads();


    if (sSPw != 1) {
        for (unsigned int mi = 0; mi<sM; mi++) {
            full_field_xypw[thread_idx].x = full_field_xypw[thread_idx].x + normalized_mode_space_profiles_xypwm[thread_idx+N*mi]*field_wi[mi].x/norm_wi[mi];
            full_field_xypw[thread_idx].y = full_field_xypw[thread_idx].y + normalized_mode_space_profiles_xypwm[thread_idx+N*mi]*field_wi[mi].y/norm_wi[mi];
        }
    } else {
        for (unsigned int mi = 0; mi<sM; mi++) {
            full_field_xypw[thread_idx].x = full_field_xypw[thread_idx].x + normalized_mode_space_profiles_xypwm[threadIdx.x+N*mi]*field_wi[mi].x/norm_wi[mi];
            full_field_xypw[thread_idx].y = full_field_xypw[thread_idx].y + normalized_mode_space_profiles_xypwm[threadIdx.x+N*mi]*field_wi[mi].y/norm_wi[mi];
        }
    }
}