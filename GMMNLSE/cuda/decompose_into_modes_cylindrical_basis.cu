__global__ void decompose_into_modes_cylindrical_basis(double2* field_wm,
                                                       const double* mode_space_profiles_xym, const double* norm_spatial_modes, const double2* full_field_xyw,
                                                       const double* r, const double dr, const double dtheta,
                                                       const unsigned int sM, const unsigned int sX, const unsigned int sY, const unsigned int sW) {
    unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x; // blockDim.x = sX*sY

    if (threadIdx.x >= (sX*sY)) return; // XYi = threadIdx.x
    if (blockIdx.x >= sW) return; // wi = blockIdx.x
    if (thread_idx >= (sX*sY*sW)) return;

    const unsigned int ri = threadIdx.x%sX;

    const unsigned int max_num_r = 100;
    const unsigned int max_sXsYsM = 1024*20;
    __shared__ double2 sum_field_xym[max_sXsYsM]; //[sX*sY*sM]
    __shared__ double this_r[max_num_r];
    if (threadIdx.x < sX)
        this_r[threadIdx.x] = r[threadIdx.x];
    __syncthreads();

    // blockDim.x = sX*sY
    for (unsigned int mi = 0; mi<sM; mi++) {
        sum_field_xym[threadIdx.x+blockDim.x*mi].x = mode_space_profiles_xym[threadIdx.x+blockDim.x*mi]*full_field_xyw[thread_idx].x*this_r[ri];
        sum_field_xym[threadIdx.x+blockDim.x*mi].y = mode_space_profiles_xym[threadIdx.x+blockDim.x*mi]*full_field_xyw[thread_idx].y*this_r[ri];
    }
    __syncthreads();

    // sM <= blockDim.x = sX*sY
    // Each thread in each block deals with the summation of each mode with its corresponding frequency now,
    // that is, threadIdx.x = mi here.
    if (threadIdx.x < sM) {
        for (unsigned int XYi = 0; XYi<blockDim.x; XYi++) {
            field_wm[blockIdx.x+sW*threadIdx.x].x = field_wm[blockIdx.x+sW*threadIdx.x].x + sum_field_xym[XYi+blockDim.x*threadIdx.x].x;
            field_wm[blockIdx.x+sW*threadIdx.x].y = field_wm[blockIdx.x+sW*threadIdx.x].y + sum_field_xym[XYi+blockDim.x*threadIdx.x].y;
        }
        field_wm[blockIdx.x+sW*threadIdx.x].x = field_wm[blockIdx.x+sW*threadIdx.x].x*dr*dtheta*norm_spatial_modes[threadIdx.x];
        field_wm[blockIdx.x+sW*threadIdx.x].y = field_wm[blockIdx.x+sW*threadIdx.x].y*dr*dtheta*norm_spatial_modes[threadIdx.x];
    }
}