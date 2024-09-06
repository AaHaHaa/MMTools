#define MAX_NUM_MODES 32 // the maximum number of modes for this cuda = sqrt(MaxThreadsPerBlock)
                         //                                           = sqrt(1024) for our Titan XP GPU

__global__ void GMMNLSE_sponRS_sum(double2* Ra,
                                   const double2* At, const double2* At_sponRS,
                                   const double* SRa,
                                   const unsigned char* nonzero_midx1234s,
                                   const unsigned int* beginning_nonzero, const unsigned int* ending_nonzero,
                                   const unsigned int N, const unsigned int M,
                                   const unsigned int NUM_MODES) {
    const unsigned int midx1 = threadIdx.x / NUM_MODES;
    const unsigned int midx2 = threadIdx.x - midx1*NUM_MODES;

    const unsigned int Midx = blockIdx.x / N;
    const unsigned int Nidx = blockIdx.x - Midx*N;

    const unsigned int NM = N*M;
    const unsigned int NMMODES = NM*NUM_MODES;

    // Preload At to improve the performance (avoiding accessing the global memory too many times)
    __shared__ double2 this_At[MAX_NUM_MODES], this_At_sponRS[MAX_NUM_MODES];
    if (midx1 == 0) {
        this_At[midx2] = At[Nidx+Midx*N+midx2*NM];
        this_At_sponRS[midx2] = At_sponRS[Nidx+Midx*N+midx2*NM];
    }
    __syncthreads();

    const unsigned int this_beginning_nonzero = beginning_nonzero[midx2+midx1*NUM_MODES];
    const unsigned int this_ending_nonzero = ending_nonzero[midx2+midx1*NUM_MODES];

    unsigned int midx3, midx4;
    double c, d, e, f; // this_At
    double p, q, r, s; // this_At_sponRS
    // compute the spontaneous Raman term
    if (this_beginning_nonzero > 0) {
        double2 this_Ra;
        this_Ra.x = 0; this_Ra.y = 0; // initialized
        for (int i = this_beginning_nonzero-1; i < this_ending_nonzero-1; i++) {
            midx3 = nonzero_midx1234s[2+i*4]-1;
            midx4 = nonzero_midx1234s[3+i*4]-1;

            c = this_At[midx3].x;
            d = this_At[midx3].y;
            e = this_At[midx4].x;
            f = this_At[midx4].y;

            p = this_At_sponRS[midx3].x;
            q = this_At_sponRS[midx3].y;
            r = this_At_sponRS[midx4].x;
            s = this_At_sponRS[midx4].y;

            if (midx3 == midx4) {
                this_Ra.x += SRa[i]*( (p*r+q*s)   + (c*r+d*s)*2 );
            } else {
                this_Ra.x += SRa[i]*( (p*r+q*s)*2 + (c*r+d*s)*2+(e*p+f*q)*2 );
            }
        }
        Ra[Nidx+Midx*N+midx1*NM+midx2*NMMODES] = this_Ra;
    }
}