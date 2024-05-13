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
*         psi(:,:,:,Midx+1,:) = psi(:,:,:,1,:) + sum(nonlinear_kw(:,:,:,1:Midx+1,:).*sim.MPA.coeff(Midx,:,:,1:Midx+1),4);
*     end
*
* However, because of the "for" loop, it's slow. I use this cuda file to finish the computation once.
*
* It's parallelized w.r.t. each entry of psi(:,:,:,2:end,:), so the number of parallelization is Nf*NNx*Ny*sim.MPA.M*Np.
*
* Input related arguments:
*    psi:  (Nf, Nx, Ny, M+1, Np)
*    nonlinear: (Nf, Nx, Ny, M+1, Np); nonlinear term in 3D-GNLSE
*    MPA_coeff: coefficients for doing the summation in the Adam-Moulton method
*    Nf: a scalar integer; the number of time/frequency points
*    Nx: a scalar integer; the number of spatial points in x
*    Ny: a scalar integer; the number of spatial points in y
*    M: a scalar integer; the number of parallelization of the MPA stepping algorithm
*    Np: a scalar integer; the number of polarization modes
* -----------------------------------------------------------------------*/
#define MAX_M 30 // the maximum number of parallelization planes

__global__ void MPA_psi_update(double2* psi,
                               const double2* nonlinear,
                               const double* MPA_coeff,
                               const unsigned int Nf, const unsigned int Nx, const unsigned int Ny, const unsigned int M, const unsigned int Np) {
    // Calculate the index
    const unsigned int Midx = threadIdx.x;

    const unsigned int Npidx = blockIdx.x / (Nf*Nx*Ny);
    const unsigned int Nyidx = (blockIdx.x-Npidx*Nf*Nx*Ny) / (Nf*Nx);
    const unsigned int Nxidx = (blockIdx.x-Npidx*Nf*Nx*Ny-Nyidx*Nf*Nx) / Nf;
    const unsigned int Nfidx = blockIdx.x-Npidx*Nf*Nx*Ny-Nyidx*Nf*Nx-Nxidx*Nf;
    
    const unsigned int idx = Nfidx+Nxidx*Nf+Nyidx*Nf*Nx+Midx*Nf*Nx*Ny+Npidx*Nf*Nx*Ny*M;
    const unsigned int idx0 = Nfidx+Nxidx*Nf+Nyidx*Nf*Nx+Npidx*Nf*Nx*Ny*M;

    // Preload to improve the performance (avoiding accessing the global memory too many times)
    __shared__ double2 this_nonlinear[MAX_M];
    this_nonlinear[Midx] = nonlinear[idx];
    __shared__ double this_MPA_coeff[(MAX_M-1)*MAX_M];
    for (int i=0; i < M-1; i++) {
        this_MPA_coeff[i+Midx*(M-1)] = MPA_coeff[i+Midx*(M-1)];
    }
    __syncthreads();

    // Calculate the summation with Adam-Moulton method
    double2 nonlinear_sum;
    nonlinear_sum.x = 0; nonlinear_sum.y = 0; // initialized
    if (Midx != 0) {
        for (int i = 0; i <= Midx; i++) {
            nonlinear_sum.x += this_MPA_coeff[(Midx-1)+i*(M-1)]*this_nonlinear[i].x;
            nonlinear_sum.y += this_MPA_coeff[(Midx-1)+i*(M-1)]*this_nonlinear[i].y;
        }
        psi[idx].x = psi[idx0].x + nonlinear_sum.x;
        psi[idx].y = psi[idx0].y + nonlinear_sum.y;
    }
}