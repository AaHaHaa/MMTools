/* ------------------------------------------------------------------------
* ---------------------------- MPA_psi_update -----------------------------
* -------------------------------------------------------------------------
* In the MPA stepping algorithm for pulse propagation, "psi" needs to be calculated.
*
* "psi" is calculated with a multistep Adams-Moulton method for each parallelization plane.
* A M-step Adams-Moulton method (with the summation of M's derivatives) has a (M+1)th-order accuracy.
* With it, the latter psi will not have a bigger error.
* The typical trapezoidal quadrature is based on Euler method which has only a first-order accuracy.
* The error builds up for latter parallization planes.
* The highest-order accuracy is 2nd order O(h^2) from
*   y(n+1)=y(n)+1/2*h*( f(t(n+1),y(n+1)) + f(t(n),y(n)) ),
* so the error estimate of adaptive-step-size control uses cubic root: error = O(h^3).
*
* Below is the MATLAB code for this update:
*     for Midx = 1:sim.MPA.M
*         psi(:,:,Midx+1) = psi(:,:,1) + sum(nonlinear.*sim.MPA.coeff(Midx,:,1:Midx+1),3);
*     end
*
* However, because of the "for" loop, it's slow. I use this cuda file to finish the computation once.
*
* It's parallelized w.r.t. each entry of psi(:,:,2:end), so the number of parallelization is N*NUM_MODES*sim.MPA.M.
*
* Input related arguments:
*    psi:  (N, NUM_MODES, M+1)
*    nonlinear: (NUM_MODES, NUM_MODES, M+1); linear operator
*    MPA_coeff: coefficients for doing the summation in the Adam-Moulton method
*    N: a scalar integer; the number of time/frequency points
*    M: a scalar integer; the number of parallelization of the MPA stepping algorithm
*    NUM_MODES: a scalar integer; the number of modes
* -----------------------------------------------------------------------*/
#define MAX_NUM_MODES 32 // the maximum number of modes for this cuda = sqrt(MaxThreadsPerBlock)
                         //                                           = sqrt(1024) for our Titan XP GPU
#define MAX_M 30 // the maximum number of parallelization planes

__global__ void MPA_psi_update(double2* psi,
                               const double2* nonlinear,
                               const double* MPA_coeff,
                               const unsigned int N, const unsigned int M, const unsigned int NUM_MODES) {
    // Calculate the index
    const unsigned int Midx = threadIdx.x / NUM_MODES;
    const unsigned int midx = threadIdx.x - Midx*NUM_MODES;    
    const unsigned int Nidx = blockIdx.x;
    
    const unsigned int NMODE = N*NUM_MODES;

    // Preload to improve the performance (avoiding accessing the global memory too many times)
    __shared__ double2 this_nonlinear[MAX_M*MAX_NUM_MODES];
    this_nonlinear[Midx+midx*M] = nonlinear[Nidx+midx*N+Midx*NMODE];
    __shared__ double this_MPA_coeff[(MAX_M-1)*MAX_M];
    if (midx == 0) {
        for (int i=0; i < M-1; i++) {
            this_MPA_coeff[i+Midx*(M-1)] = MPA_coeff[i+Midx*(M-1)];
        }
    }
    __syncthreads();

    // Calculate the summation with Adam-Moulton method
    double2 nonlinear_sum;
    nonlinear_sum.x = 0; nonlinear_sum.y = 0; // initialized
    if (Midx != 0) {
        for (int i = 0; i <= Midx; i++) {
            nonlinear_sum.x += this_MPA_coeff[(Midx-1)+i*(M-1)]*this_nonlinear[i+midx*M].x;
            nonlinear_sum.y += this_MPA_coeff[(Midx-1)+i*(M-1)]*this_nonlinear[i+midx*M].y;
        }
        psi[Nidx+midx*N+Midx*NMODE].x = psi[Nidx+midx*N].x + nonlinear_sum.x;
        psi[Nidx+midx*N+Midx*NMODE].y = psi[Nidx+midx*N].y + nonlinear_sum.y;
    }
}