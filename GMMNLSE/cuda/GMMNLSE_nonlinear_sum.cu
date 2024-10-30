#define MAX_NUM_MODES 32 // the maximum number of modes for this cuda = sqrt(MaxThreadsPerBlock)
                         //                                           = sqrt(1024) for our Titan XP GPU

__global__ void GMMNLSE_nonlinear_sum(double2* Kerr, double2* Ra,
                                      const double2* At,
                                      const double* SK, const double* SRa,
                                      const unsigned char* nonzero_midx1234s,
                                      const unsigned int* beginning_nonzero, const unsigned int* ending_nonzero,
                                      const bool include_Raman,
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

    // Preload At to improve the performance (avoiding accessing the global memory too many times)
    __shared__ double2 this_At[MAX_NUM_MODES];
    if (midx1 == 0) this_At[midx2] = At[Nidx+Midx*N+midx2*NM];
    __syncthreads();

    const unsigned int this_beginning_nonzero = beginning_nonzero[midx2+midx1*NUM_MODES];
    const unsigned int this_ending_nonzero = ending_nonzero[midx2+midx1*NUM_MODES];

    unsigned int midx3, midx4;
    double c, d, e, f;
    switch (OPERATIONIdx) {
        case 0: // compute the Kerr term
            if (this_beginning_nonzero > 0) {
                double a, b, pcdef;
                a = this_At[midx2].x;
                b = this_At[midx2].y;

                double2 this_Kerr;
                this_Kerr.x = 0; this_Kerr.y = 0; // initialized
                for (int i = this_beginning_nonzero-1; i < this_ending_nonzero-1; i++) {
                    midx3 = nonzero_midx1234s[2+i*4]-1;
                    midx4 = nonzero_midx1234s[3+i*4]-1;
                    
                    c = this_At[midx3].x;
                    d = this_At[midx3].y;
                    e = this_At[midx4].x;
                    f = this_At[midx4].y;
            
                    pcdef = SK[i]*(c*e+d*f);
                    if (midx3 == midx4) { // d*e-c*f= 0
                        this_Kerr.x += a*pcdef;
                        this_Kerr.y += b*pcdef;
                    } else { // (d*e-c*f) + (c <--> e, d <--> f) = 0
                        this_Kerr.x += a*pcdef*2;
                        this_Kerr.y += b*pcdef*2;
                    }
                }
                Kerr[Nidx+Midx*N+midx1*NM+midx2*NMMODES] = this_Kerr;
            }
            break;

        case 1: // compute the Raman term
            if (include_Raman && this_beginning_nonzero > 0) {
                double2 this_Ra;
                this_Ra.x = 0; this_Ra.y = 0; // initialized
                for (int i = this_beginning_nonzero-1; i < this_ending_nonzero-1; i++) {
                    midx3 = nonzero_midx1234s[2+i*4]-1;
                    midx4 = nonzero_midx1234s[3+i*4]-1;
                    
                    c = this_At[midx3].x;
                    d = this_At[midx3].y;
                    e = this_At[midx4].x;
                    f = this_At[midx4].y;

                    if (midx3 == midx4) { // d*e-c*f= 0
                        this_Ra.x += SRa[i]*(c*e+d*f);
                    } else { // (d*e-c*f) + (c <--> e, d <--> f) = 0
                        this_Ra.x += SRa[i]*(c*e+d*f)*2;
                    }
                }
                Ra[Nidx+Midx*N+midx1*NM+midx2*NMMODES] = this_Ra;
            }
            break;
    }
}