/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"
#include <cblas.h>

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	double *C = (double*) malloc(N * N * sizeof(double));
	double *temp = (double*) malloc(N * N * sizeof(double));

	/* 1. Calculate B^t x B^t and store the result in temp */
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, N, N, N, 1.0, B, N, B, N, 0.0, temp, N);
    
    /* 2. Calculate A × B and store the result in B,
		  considering that A is triangular matrix */
	cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, 1.0, A, N, B, N);

    /* 3. Calculate B(new B, which is A x B) x A^t and store the result in C */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, B, N, A, N, 0.0, C, N);

	/* 4. Add to C(A x B x A^t) the temp(B^t x B^t) variable */
    cblas_daxpy(N * N, 1.0, temp, 1, C, 1);

    free(temp);
    return C;
}
