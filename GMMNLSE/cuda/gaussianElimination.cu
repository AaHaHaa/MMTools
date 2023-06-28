__device__ void gaussianElimination(double*, const unsigned int, double*);
__device__ void swap_row(double*, const unsigned int, const unsigned int, const unsigned int);
__device__ void forwardElim(double*, const unsigned int);
__device__ void backSub(const double*, const unsigned int, double*);


// function to get matrix content
__device__ void gaussianElimination(double* mat, const unsigned int N, double* x) { // mat[N][N+1]
    
    /* reduction into r.e.f. */
    forwardElim(mat,N);

    /* get solution to system and print it using backward substitution */
    backSub(mat,N,x);
}

// function for elementary operation of swapping two rows
__device__ void swap_row(double* mat, const unsigned int N, const unsigned int i, const unsigned int j) {

    for (int k=0; k<=N; k++) {
        double temp = mat[i*(N+1)+k];
        mat[i*(N+1)+k] = mat[j*(N+1)+k];
        mat[j*(N+1)+k] = temp;
    }
}

// function to reduce matrix to r.e.f.
__device__ void forwardElim(double* mat, const unsigned int N) {
    for (int k=0; k<N; k++) {
        // Initialize maximum value and index for pivot
        int i_max = k;
        int v_max = abs(mat[i_max*(N+1)+k]);
  
        /* find greater amplitude for pivot if any */
        for (int i = k+1; i < N; i++) {
            if (abs(mat[i*(N+1)+k]) > v_max) {
                v_max = abs(mat[i*(N+1)+k]);
                i_max = i;
            }
        }
  
        /* Swap the greatest value row with current row */
        if (i_max != k)
            swap_row(mat, N, k, i_max);

        for (int i=k+1; i<N; i++) {
            /* factor f to set current row kth element to 0,
             * and subsequently remaining kth column to 0 */
            double f = mat[i*(N+1)+k]/mat[k*(N+1)+k];

            /* subtract fth multiple of corresponding kth row element*/
            for (int j=k+1; j<=N; j++)
                mat[i*(N+1)+j] -= mat[k*(N+1)+j]*f;

            /* filling lower triangular matrix with zeros*/
            mat[i*(N+1)+k] = 0;
        }
    }
}

// function to calculate the values of the unknowns
__device__ void backSub(const double* mat, const unsigned int N, double* x) {
    //double x[N];  // An array to store solution

    /* Start calculating from last equation up to the first */
    for (int i = N-1; i >= 0; i--) {
        /* start with the RHS of the equation */
        x[i] = mat[i*(N+1)+N];

        /* Initialize j to i+1 since matrix is upper triangular*/
        for (int j=i+1; j<N; j++) {
            /* subtract all the lhs values except the coefficient of the 
             * variable whose value is being calculated */
            x[i] -= mat[i*(N+1)+j]*x[j];
        }
  
        /* divide the RHS by the coefficient of the unknown being calculated */
        x[i] = x[i]/mat[i*(N+1)+i];
    }
} 