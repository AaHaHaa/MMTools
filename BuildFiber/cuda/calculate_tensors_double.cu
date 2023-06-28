
__global__ void calculate_tensors(double* SR, const double* fields, const double* norms, const int num_modes, const int Nx) {
    unsigned int full_thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    // Calculate the index
    unsigned int nmp4 = num_modes*num_modes*num_modes*num_modes;
    unsigned int Nxnm = Nx*num_modes;

    if (full_thread_idx >= nmp4) {
        return;
    }

    // Turn linear index into components
    unsigned int midx1 = full_thread_idx % num_modes;
    unsigned int midx2 = (full_thread_idx/num_modes) % num_modes;
    unsigned int midx3 = (full_thread_idx/num_modes/num_modes) % num_modes;
    unsigned int midx4 = (full_thread_idx/num_modes/num_modes/num_modes);

    // Compute the sum
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Nx; j++) {
            SR[full_thread_idx] += fields[midx1+i*num_modes+j*Nxnm]*fields[midx2+i*num_modes+j*Nxnm]*fields[midx3+i*num_modes+j*Nxnm]*fields[midx4+i*num_modes+j*Nxnm];
        }
    }

    // Normalize
    SR[full_thread_idx] /= norms[midx1]*norms[midx2]*norms[midx3]*norms[midx4];
}
