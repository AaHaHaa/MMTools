// the following two lines are to use the M_PI constant=3.14...
#define _USE_MATH_DEFINES
#include <math.h>

/* ------------------------------------------------------------------------
* ------------------------------ interp1_D2 -------------------------------
* -------------------------------------------------------------------------
* This function does linear interpolation to complex numbers.
* It interpolates with the absolute values and phases of complex numbers, 
* rather than real and imaginary parts.
* I implemented this function because I realized that interpolation with 
* real and imaginary parts creates artifitial oscillations in absolute values.
* To solve this directly, interpolation should be done with absolute values 
* and phases intead. Because it uses phases, phase unwrapping is required; 
* a tolerance of phase unwrapping needs to be given, which is typically 0.5.
*
* It's parallelized w.r.t. each query point, so the number of parallelization is total_size.
*
* Input related arguments:
*    x0: (Nx0, 1) or (Nx0, M); sample points
*    A0: (Nx0, 1) or (Nx0, M); sample values
*    multiple_x0: true or false;
*                 whether each column of A0 corresponds to different x0 or not.
*                 If it's true, both x0 and A0 need to have the size (Nx0, M).
*    Nx0: a scalar integer; the length of x0
*    tol: a scalar between 0 and 1; the tolerance of the phase for phase unwrapping
*                                   (typically 0.5 which means 2*pi*0.5 = pi).
* Output related arguments:
*    x: (Nx, 1) or (Nx, M); query points
*    A: (Nx, 1) or (Nx, M); query values
*    multiple_x: true or false;
*                whether each column of A corresponds to different x or not.
*                If it's true, both x and A need to have the size (Nx, M).
*                If both multiple_x0 and multiple_x are true, each column has an independent computation of interpolations.
*    Nx: a scalar integer; the length of x
*    total_size: a scalar integer; the total size of A = Nx or Nx*M
* -----------------------------------------------------------------------*/
__global__ void interp1_D2(const double* x0, const double2* A0, const bool multiple_x0, const unsigned int Nx0,
                           const double* x,  double2* A, const bool multiple_x, const unsigned int Nx,
                           const unsigned int total_size,
                           const double tol) {
    const unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (thread_idx >= total_size) return;

    const unsigned int higherDim_idx = thread_idx / Nx;
    const unsigned int Ni = thread_idx - higherDim_idx*Nx; // thread_idx % Nx

    unsigned int base_idx_x0;
    if (multiple_x0) {
        base_idx_x0 = Nx0*higherDim_idx;
    } else {
        base_idx_x0 = 0;
    }
    const unsigned int base_idx_A0 = Nx0*higherDim_idx;

    unsigned int base_idx_x;
    if (multiple_x) {
        base_idx_x = Nx*higherDim_idx;
    } else {
        base_idx_x = 0;
    }

   /* -----------------------------------------------------------------------
    * 1) Find which interval targeted x should fit in with binary search
    * ----------------------------------------------------------------------- */
    unsigned int upper_idx, lower_idx;

    unsigned int idx;
    double target_x0;
    const double target_x = x[base_idx_x+Ni];
    unsigned int min_idx = 0;
    unsigned int max_idx = Nx0-1;
    double min_x = x0[base_idx_x0+min_idx];
    double max_x = x0[base_idx_x0+max_idx];
    
    // target_x is out of range
    if (target_x > max_x) {
        upper_idx = max_idx;
        lower_idx = max_idx-1;
    } else if (target_x < min_x) {
        upper_idx = min_idx+1;
        lower_idx = min_idx;
    } else {
        // target_x is within range
        if (target_x == min_x) {
            A[thread_idx] = A0[base_idx_A0+min_idx];
            return;
        }
        if (target_x == max_x) {
            A[thread_idx] = A0[base_idx_A0+max_idx];
            return;
        }
        while (true) {
            idx = floor(((double)max_idx+(double)min_idx)/2);
            target_x0 = x0[base_idx_x0+idx];
            if (target_x0 > target_x) {
                if (target_x > x0[base_idx_x0+idx-1]) {
                    upper_idx = idx;
                    lower_idx = idx-1;
                    break;
                } else if (target_x == x0[base_idx_x0+idx-1]) {
                    A[thread_idx] = A0[base_idx_A0+idx-1];
                    return;
                }
                max_idx = idx;
                max_x = target_x0;
            } else if (target_x0 < target_x) {
                if (target_x < x0[base_idx_x0+idx+1]) {
                    upper_idx = idx+1;
                    lower_idx = idx;
                    break;
                } else if (target_x == x0[base_idx_x0+idx+1]) {
                    A[thread_idx] = A0[base_idx_A0+idx+1];
                    return;
                }
                min_idx = idx;
                min_x = target_x0;
            } else { // target_x0 == target_x
                A[thread_idx] = A0[base_idx_A0+idx];
                return;
            }
        }
    }

   /* -----------------------------------------------------------------------
    * 2) start the interpolation
    * ----------------------------------------------------------------------- */
    // coefficients for interpolation
    const double c1 = (target_x-x0[base_idx_x0+lower_idx])/(x0[base_idx_x0+upper_idx]-x0[base_idx_x0+lower_idx]);
    const double c2 = (x0[base_idx_x0+upper_idx]-target_x)/(x0[base_idx_x0+upper_idx]-x0[base_idx_x0+lower_idx]);

    const double r0_upper = sqrt(pow(A0[base_idx_A0+upper_idx].x,2) + pow(A0[base_idx_A0+upper_idx].y,2));
    const double r0_lower = sqrt(pow(A0[base_idx_A0+lower_idx].x,2) + pow(A0[base_idx_A0+lower_idx].y,2));
    double angle0_upper = atan2(A0[base_idx_A0+upper_idx].y,A0[base_idx_A0+upper_idx].x);
    double angle0_lower = atan2(A0[base_idx_A0+lower_idx].y,A0[base_idx_A0+lower_idx].x);
    if ((angle0_upper - angle0_lower) > 2*M_PI*tol) {
        angle0_upper -= 2*M_PI;
    } else if ((angle0_lower - angle0_upper) > 2*M_PI*tol) {
        angle0_upper += 2*M_PI;
    }
    
    const double target_r = r0_upper*c1 + r0_lower*c2;
    double target_angle;
    if (A0[base_idx_A0+upper_idx].y == 0 && A0[base_idx_A0+lower_idx].y == 0) { // real value
        target_angle = 0;
    } else {
        target_angle = angle0_upper*c1 + angle0_lower*c2;
    }
    A[thread_idx].x = target_r*cos(target_angle);
    A[thread_idx].y = target_r*sin(target_angle);
}