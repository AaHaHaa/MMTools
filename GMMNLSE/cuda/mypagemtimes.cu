/* ------------------------------------------------------------------------
* ----------------------------- mypagemtimes ------------------------------
* -------------------------------------------------------------------------
* This function does simultaneous many matrix multiplications, mtimes, for A_out = D*A_in.
* D is a linear operator that is applied to the modes A_in for each N and M pages.
* In MATLAB language, this can be done with its pagewised matrix multiplication:
*    A_out = permute(pagemtimes(D,permute(A_in,[2,4,3,1])),[4,1,3,2]);
*
* It's parallelized w.r.t. each entry of A_in, so the number of parallelization is N*NUM_MODES*M, the size of A_in.
*
* Input related arguments:
*    A_in:  (N, NUM_MODES, M); input fields
*    D: (NUM_MODES, NUM_MODES, M); linear operator
*    N: a scalar integer; the number of time/frequency points
*    M: a scalar integer; the number of parallelization of the MPA stepping algorithm
*    NUM_MODES: a scalar integer; the number of modes
* Output related arguments:
*    A_out: (N, NUM_MODES, M); output fields
* -----------------------------------------------------------------------*/
#define MAX_NUM_MODES 32 // the maximum number of modes for this cuda = MaxThreadsPerBlock (max)
                         //                                           = 1024 for our Titan XP GPU
                         // I pick 32 here because the bottleneck lies in GMMNLSE_nonlinear_sum.cu where only sqrt(1024) modes are allowed.

__global__ void mypagemtimes(double2* A_out,
                             const double2* A_in,
                             const double2* D,
                             const unsigned int N, const unsigned int M,
                             const unsigned int NUM_MODES) {
    if (blockIdx.x >= N*M) return;
    
    const unsigned int midx = threadIdx.x;
    const unsigned int Midx = blockIdx.x / N;
    const unsigned int Nidx = blockIdx.x - Midx*N;
    
    const unsigned int MODE2 = NUM_MODES*NUM_MODES;
    const unsigned int NMODE = N*NUM_MODES;

    // Preload A_in and D to improve the performance (avoiding accessing the global memory too many times)
    __shared__ double2 this_A[MAX_NUM_MODES];
    this_A[midx] = A_in[Nidx+midx*N+Midx*NMODE];

    __shared__ double2 this_D[MAX_NUM_MODES*MAX_NUM_MODES];
    for (int i = 0; i < NUM_MODES; i++) {
        this_D[i+midx*NUM_MODES] = D[i+midx*NUM_MODES+Midx*MODE2];
    }
    __syncthreads();

    double2 A_out_NMi;
    A_out_NMi.x = 0; A_out_NMi.y = 0; // initialized
    // Compute D*A_in
    for (int i = 0; i < NUM_MODES; i++) {
        A_out_NMi.x += (this_D[i*NUM_MODES+midx].x*this_A[i].x - this_D[i*NUM_MODES+midx].y*this_A[i].y);
        A_out_NMi.y += (this_D[i*NUM_MODES+midx].x*this_A[i].y + this_D[i*NUM_MODES+midx].y*this_A[i].x);
    }
    A_out[Nidx+midx*N+Midx*NMODE] = A_out_NMi;
}