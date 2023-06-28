#include "mySpline_helper.cu"

/* ------------------------------------------------------------------------
* ------------------------------- mySpline --------------------------------
* -------------------------------------------------------------------------
* This function calculates the spline interpolation of the num_order-th order.
* 
* For each interval of sampling points, the interpolation is found by fitting
* a polynomial of order 2*(num_order+1)-1 to the sampling values and their derivatives.
* Its interpolation aims to maintain the C^num_order continuity of the sampling points.
*
* It's parallelized w.r.t. each query point, so the number of parallelization is Nx.
*
* E.g. 
*    Linear interpolation is "num_order = 0" of mySpline(); not implemented here.
*    Cubic spline interpolation is "num_order = 1" of mySpline().
*
* Input related arguments:
*    x0: (Nx0, 1); sample points
*    A0: (Nx0, num_order+1); sample values
*    Nx0: a scalar integer; the length of x0
*    num_order: a scalar integer from 1 to 6 (up to 6th-order spline is implemented)
* Output related arguments:
*    x: (Nx, 1); query points
*    A: (Nx, 1); query values
*    Nx: a scalar integer; the length of x
* -----------------------------------------------------------------------*/
__global__ void mySpline(const double* x0, const double2* A0, const unsigned int Nx0,
                         const double* x,  double2* A, const unsigned int Nx,
                         const unsigned int num_order) {
    const unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (thread_idx >= Nx) return;

    const unsigned int Ni = thread_idx;

   /* -----------------------------------------------------------------------
    * 1) Find which interval targeted x should fit in with binary search
    * ----------------------------------------------------------------------- */
    unsigned int upper_idx, lower_idx;

    unsigned int idx;
    double target_x0;
    const double target_x = x[Ni];
    unsigned int min_idx = 0;
    unsigned int max_idx = Nx0-1;
    double min_x = x0[min_idx];
    double max_x = x0[max_idx];
    
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
            A[thread_idx] = A0[min_idx];
            return;
        }
        if (target_x == max_x) {
            A[thread_idx] = A0[max_idx];
            return;
        }
        while (true) {
            idx = floor(((double)max_idx+(double)min_idx)/2);
            target_x0 = x0[idx];
            if (target_x0 > target_x) {
                if (target_x > x0[idx-1]) {
                    upper_idx = idx;
                    lower_idx = idx-1;
                    break;
                } else if (target_x == x0[idx-1]) {
                    A[thread_idx] = A0[idx-1];
                    return;
                }
                max_idx = idx;
                max_x = target_x0;
            } else if (target_x0 < target_x) {
                if (target_x < x0[idx+1]) {
                    upper_idx = idx+1;
                    lower_idx = idx;
                    break;
                } else if (target_x == x0[idx+1]) {
                    A[thread_idx] = A0[idx+1];
                    return;
                }
                min_idx = idx;
                min_x = target_x0;
            } else { // target_x0 == target_x
                A[thread_idx] = A0[idx];
                return;
            }
        }
    }

   /* -----------------------------------------------------------------------
    * 2) start the interpolation
    * ----------------------------------------------------------------------- */
    switch (num_order) {
        case 1:
            A[thread_idx] = call_order1(target_x-x0[lower_idx],x0[upper_idx]-x0[lower_idx], \
                                        A0[lower_idx],A0[upper_idx], \
                                        A0[lower_idx+Nx0],A0[upper_idx+Nx0]);
            break;
        case 2:
            A[thread_idx] = call_order2(target_x-x0[lower_idx],x0[upper_idx]-x0[lower_idx], \
                                        A0[lower_idx],A0[upper_idx], \
                                        A0[lower_idx+Nx0],A0[upper_idx+Nx0], \
                                        A0[lower_idx+Nx0*2],A0[upper_idx+Nx0*2]);
            break;
        case 3:
            A[thread_idx] = call_order3(target_x-x0[lower_idx],x0[upper_idx]-x0[lower_idx], \
                                        A0[lower_idx],A0[upper_idx], \
                                        A0[lower_idx+Nx0],A0[upper_idx+Nx0], \
                                        A0[lower_idx+Nx0*2],A0[upper_idx+Nx0*2], \
                                        A0[lower_idx+Nx0*3],A0[upper_idx+Nx0*3]);
            break;
        case 4:
            A[thread_idx] = call_order4(target_x-x0[lower_idx],x0[upper_idx]-x0[lower_idx], \
                                        A0[lower_idx],A0[upper_idx], \
                                        A0[lower_idx+Nx0],A0[upper_idx+Nx0], \
                                        A0[lower_idx+Nx0*2],A0[upper_idx+Nx0*2], \
                                        A0[lower_idx+Nx0*3],A0[upper_idx+Nx0*3], \
                                        A0[lower_idx+Nx0*4],A0[upper_idx+Nx0*4]);
            break;
        case 5:
            A[thread_idx] = call_order5(target_x-x0[lower_idx],x0[upper_idx]-x0[lower_idx], \
                                        A0[lower_idx],A0[upper_idx], \
                                        A0[lower_idx+Nx0],A0[upper_idx+Nx0], \
                                        A0[lower_idx+Nx0*2],A0[upper_idx+Nx0*2], \
                                        A0[lower_idx+Nx0*3],A0[upper_idx+Nx0*3], \
                                        A0[lower_idx+Nx0*4],A0[upper_idx+Nx0*4], \
                                        A0[lower_idx+Nx0*5],A0[upper_idx+Nx0*5]);
            break;
        case 6:
            A[thread_idx] = call_order6(target_x-x0[lower_idx],x0[upper_idx]-x0[lower_idx], \
                                        A0[lower_idx],A0[upper_idx], \
                                        A0[lower_idx+Nx0],A0[upper_idx+Nx0], \
                                        A0[lower_idx+Nx0*2],A0[upper_idx+Nx0*2], \
                                        A0[lower_idx+Nx0*3],A0[upper_idx+Nx0*3], \
                                        A0[lower_idx+Nx0*4],A0[upper_idx+Nx0*4], \
                                        A0[lower_idx+Nx0*5],A0[upper_idx+Nx0*5], \
                                        A0[lower_idx+Nx0*6],A0[upper_idx+Nx0*6]);
            break;
    }
}