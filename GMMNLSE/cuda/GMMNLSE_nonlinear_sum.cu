#define MAX_NUM_MODES 32 // the maximum number of modes for this cuda = sqrt(MaxThreadsPerBlock)
                         //                                           = sqrt(1024) for our Titan XP GPU

__global__ void GMMNLSE_nonlinear_sum(double2* Kerr, double2* Ra, double2* Ra_sponRS,
                                      const double2* At, const double2* At_noise,
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
    __shared__ double2 this_At[MAX_NUM_MODES], this_At_noise[MAX_NUM_MODES];
    switch (OPERATIONIdx) {
        case 0: // For Kerr interactions, noise photon is included directly for accurately computing noise-seeded processes
            if (midx1 == 0) {
                this_At[midx2].x = At[Nidx+Midx*N+midx2*NM].x + At_noise[Nidx+Midx*N+midx2*NM].x;
                this_At[midx2].y = At[Nidx+Midx*N+midx2*NM].y + At_noise[Nidx+Midx*N+midx2*NM].y;
            }
            break;
        case 1:
            if (midx1 == 0) this_At[midx2] = At[Nidx+Midx*N+midx2*NM];
            break;
        case 2:
            if (midx1 == 0) {
                this_At[midx2] = At[Nidx+Midx*N+midx2*NM];
                this_At_noise[midx2] = At_noise[Nidx+Midx*N+midx2*NM];
            }
            break;
    }
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

        case 2: // compute the spontaneous Raman term
            if (include_Raman && this_beginning_nonzero > 0) {
                double p, q, r, s; // this_At_noise
                double2 this_Ra_sponRS;
                this_Ra_sponRS.x = 0; this_Ra_sponRS.y = 0; // initialized
                for (int i = this_beginning_nonzero-1; i < this_ending_nonzero-1; i++) {
                    midx3 = nonzero_midx1234s[2+i*4]-1;
                    midx4 = nonzero_midx1234s[3+i*4]-1;
        
                    c = this_At[midx3].x;
                    d = this_At[midx3].y;
                    e = this_At[midx4].x;
                    f = this_At[midx4].y;
        
                    p = this_At_noise[midx3].x;
                    q = this_At_noise[midx3].y;
                    r = this_At_noise[midx4].x;
                    s = this_At_noise[midx4].y;
        
                    if (midx3 == midx4) {
                        this_Ra_sponRS.x += SRa[i]*( (p*r+q*s)   + (c*r+d*s)*2 );
                    } else {
                        this_Ra_sponRS.x += SRa[i]*( (p*r+q*s)*2 + (c*r+d*s)*2+(e*p+f*q)*2 );
                    }
                }
                Ra_sponRS[Nidx+Midx*N+midx1*NM+midx2*NMMODES] = this_Ra_sponRS;
            }
            break;
    }
}