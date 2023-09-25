#define MAX_NUM_MODES 32 // the maximum number of modes for this cuda = sqrt(MaxThreadsPerBlock)
                         //                                           = sqrt(1024) for our Titan XP GPU

__global__ void GMMNLSE_nonlinear_sum_MMGaussianGain_with_polarization(double2* Kerr, double2* Ra, double2* Rb, double2* transfer_matrix,
                                                                       const double2* A_t,
                                                                       const double2* Bmn,
                                                                       const double* SK,  const unsigned char* SK_nonzero_midx1234s,  const unsigned int* SK_beginning_nonzero,  const unsigned int* SK_ending_nonzero,
                                                                       const double* SRa, const unsigned char* SRa_nonzero_midx1234s, const unsigned int* SRa_beginning_nonzero, const unsigned int* SRa_ending_nonzero,
                                                                       const double* SRb, const unsigned char* SRb_nonzero_midx1234s, const unsigned int* SRb_beginning_nonzero, const unsigned int* SRb_ending_nonzero,
                                                                       const bool include_Raman, const bool include_anisoRaman,
                                                                       const unsigned int N, const unsigned int M,
                                                                       const unsigned int NUM_MODES,
                                                                       const unsigned int NUM_OPERATIONS) {
    const unsigned int midx1 = threadIdx.x / NUM_MODES;
    const unsigned int midx2 = threadIdx.x - midx1*NUM_MODES;

    const unsigned int NMIdx = blockIdx.x / NUM_OPERATIONS;
    const unsigned int OPERATIONIdx = blockIdx.x - NMIdx*NUM_OPERATIONS;

    const unsigned int Midx = NMIdx / N;
    const unsigned int Nidx = NMIdx - Midx*N;

    const unsigned int NM = N*M;
    const unsigned int NMMODES = NM*NUM_MODES;
    const unsigned int NUM_SPATIAL_MODES = NUM_MODES / 2;
    const unsigned int SMODES2 = NUM_SPATIAL_MODES*NUM_SPATIAL_MODES;

    // Preload A_t to improve the performance (avoiding accessing the global memory too many times)
    __shared__ double2 this_A[MAX_NUM_MODES];
    if (midx1 == 0) this_A[midx2] = A_t[Nidx+Midx*N+midx2*NM];
    __syncthreads();

    const unsigned int this_SK_beginning_nonzero = SK_beginning_nonzero[midx1+midx2*NUM_MODES];
    const unsigned int this_SK_ending_nonzero = SK_ending_nonzero[midx1+midx2*NUM_MODES];
    const unsigned int this_SRa_beginning_nonzero = SRa_beginning_nonzero[midx1+midx2*NUM_MODES];
    const unsigned int this_SRa_ending_nonzero = SRa_ending_nonzero[midx1+midx2*NUM_MODES];
    const unsigned int this_SRb_beginning_nonzero = SRb_beginning_nonzero[midx1+midx2*NUM_MODES];
    const unsigned int this_SRb_ending_nonzero = SRb_ending_nonzero[midx1+midx2*NUM_MODES];

    unsigned int midx3, midx4;
    double a, b, c, d, e, f, pcdef, ncdef;
    switch (OPERATIONIdx) {
        case 0: // compute the Kerr term
            if (this_SK_beginning_nonzero > 0) {
                a = this_A[midx2].x;
                b = this_A[midx2].y;

                double2 this_Kerr;
                this_Kerr.x = 0; this_Kerr.y = 0; // initialized
                for (int i = this_SK_beginning_nonzero-1; i < this_SK_ending_nonzero-1; i++) {
                    midx3 = SK_nonzero_midx1234s[2+i*4]-1;
                    midx4 = SK_nonzero_midx1234s[3+i*4]-1;
            
                    c = this_A[midx3].x;
                    d = this_A[midx3].y;
                    e = this_A[midx4].x;
                    f = this_A[midx4].y;
            
                    pcdef = SK[i]*(c*e+d*f);
                    if (midx3 == midx4 || (int(midx3 & 1) != int(midx4 & 1)) ) {
                        if (midx3 == midx4) { // c=e, d=f --> ncdef=0
                            this_Kerr.x += a*pcdef;
                            this_Kerr.y += b*pcdef;
                        } else {
                            ncdef = SK[i]*(c*f-d*e);
                            this_Kerr.x += a*pcdef+b*ncdef;
                            this_Kerr.y += b*pcdef-a*ncdef;
                        }
                    } else {
                        this_Kerr.x += a*pcdef*2;
                        this_Kerr.y += b*pcdef*2;
                    }
                }
                Kerr[Nidx+Midx*N+midx1*NM+midx2*NMMODES] = this_Kerr;
            }
            break;

        case 1: // compute the SRa tensors, isotropic Raman response
            if (include_Raman && this_SRa_beginning_nonzero > 0) {
                double2 this_Ra;
                this_Ra.x = 0; this_Ra.y = 0; // initialized
                for (int i = this_SRa_beginning_nonzero-1; i < this_SRa_ending_nonzero-1; i++) {
                    midx3 = SRa_nonzero_midx1234s[2+i*4]-1;
                    midx4 = SRa_nonzero_midx1234s[3+i*4]-1;
        
                    c = this_A[midx3].x;
                    d = this_A[midx3].y;
                    e = this_A[midx4].x;
                    f = this_A[midx4].y;

                    if (midx3 == midx4 || (int(midx3 & 1) != int(midx4 & 1)) ) {
                        if (midx3 == midx4) { // c=e, d=f
                            this_Ra.x += SRa[i]*(c*e+d*f);
                        } else {
                            this_Ra.x += SRa[i]*(c*e+d*f);
                            this_Ra.y += SRa[i]*(d*e-c*f);
                        }
                    } else {
                        this_Ra.x += SRa[i]*(c*e+d*f)*2;
                    }
                }
                Ra[Nidx+Midx*N+midx1*NM+midx2*NMMODES] = this_Ra;
            }
            break;

        case 2: // compute the transfer matrix
            if (this_SRa_beginning_nonzero > 0 && Nidx == 0) { // "Nidx==0" is to limit this to run only once through different Nidx
                unsigned int spatial_midx1, spatial_midx2, spatial_midx3, spatial_midx4;
                unsigned int polar_idx_T, polar_idx_B;
                unsigned int idx2_T, idx3_T, idx4_T;
                unsigned int idx2_B, idx3_B, idx4_B;

                // the index related to polarizations
                polar_idx_T = midx1 % 2;
                // recovery from total indices back to spatial mode indices
                spatial_midx1 = midx1 / 2;
                spatial_midx2 = midx2 / 2;
                // the dimension of the computation/assignment below
                idx2_T = spatial_midx1*M;
                idx3_T = spatial_midx2*M*NUM_SPATIAL_MODES;
                idx4_T = polar_idx_T*M*SMODES2;

                double2 this_transfer_matrix;
                this_transfer_matrix.x = 0; this_transfer_matrix.y = 0; // initialized
                for (int i = this_SRa_beginning_nonzero-1; i < this_SRa_ending_nonzero-1; i++) {
                    midx3 = SRa_nonzero_midx1234s[2+i*4]-1;
                    midx4 = SRa_nonzero_midx1234s[3+i*4]-1;

                    if (midx3 == midx4 || (int(midx3 & 1) != int(midx4 & 1)) ) {
                        // the index related to polarizations
                        polar_idx_B = midx3 % 2;
                        // recovery from total indices back to spatial mode indices
                        spatial_midx3 = midx3 / 2;
                        spatial_midx4 = midx4 / 2;
                        // the dimension of the computation/assignment below
                        idx2_B = spatial_midx3*M;
                        idx3_B = spatial_midx4*M*NUM_SPATIAL_MODES;
                        idx4_B = polar_idx_B*M*SMODES2;

                        this_transfer_matrix.x += SRa[i]*Bmn[Midx+idx2_B+idx3_B+idx4_B].x;
                        this_transfer_matrix.y += SRa[i]*Bmn[Midx+idx2_B+idx3_B+idx4_B].y;
                    } else {
                        // the index related to polarizations
                        polar_idx_B_1 = midx3 & 1;                      polar_idx_B_2 = midx4 & 1;
                        // recovery from total indices back to spatial mode indices
                        spatial_midx3_1 = midx3 >> 1;                   spatial_midx3_2 = midx4 >> 1;
                        spatial_midx4_1 = midx4 >> 1;                   spatial_midx4_2 = midx3 >> 1;
                        // the dimension of the computation/assignment below
                        idx2_B_1 = spatial_midx3_1*M;                   idx2_B_2 = spatial_midx3_2*M;
                        idx3_B_1 = spatial_midx4_1*M*NUM_SPATIAL_MODES; idx3_B_2 = spatial_midx4_2*M*NUM_SPATIAL_MODES;
                        idx4_B_1 = polar_idx_B_1*M*SMODES2;             idx4_B_2 = polar_idx_B_2*M*SMODES2;

                        this_transfer_matrix.x += SRa[i]*(Bmn[Midx+idx2_B_1+idx3_B_1+idx4_B_1].x+Bmn[Midx+idx2_B_2+idx3_B_2+idx4_B_2].x);
                        this_transfer_matrix.y += SRa[i]*(Bmn[Midx+idx2_B_1+idx3_B_1+idx4_B_1].y+Bmn[Midx+idx2_B_2+idx3_B_2+idx4_B_2].y);
                    }
                }
                transfer_matrix[Midx+idx2_T+idx3_T+idx4_T] = this_transfer_matrix;
            }
            break;

        case 3: // compute the SRb tensors, anisotropic Raman response
            if (include_anisoRaman && this_SRb_beginning_nonzero > 0) {
                double2 this_Rb;
                this_Rb.x = 0; this_Rb.y = 0; // initialized
                for (int i = this_SRb_beginning_nonzero-1; i < this_SRb_ending_nonzero-1; i++) {
                    midx3 = SRb_nonzero_midx1234s[2+i*4]-1;
                    midx4 = SRb_nonzero_midx1234s[3+i*4]-1;
        
                    c = this_A[midx3].x;
                    d = this_A[midx3].y;
                    e = this_A[midx4].x;
                    f = this_A[midx4].y;
        
                    if (midx3 == midx4 || (int(midx3 & 1) != int(midx4 & 1)) ) {
                        if (midx3 == midx4) { // c=e, d=f
                            this_Rb.x += SRb[i]*(c*e+d*f);
                        } else {
                            this_Rb.x += SRb[i]*(c*e+d*f);
                            this_Rb.y += SRb[i]*(d*e-c*f);
                        }
                    } else {
                        this_Rb.x += SRb[i]*(c*e+d*f)*2;
                    }
                }
                Rb[Nidx+Midx*N+midx1*NM+midx2*NMMODES] = this_Rb;
            }
            break;
    }
}
