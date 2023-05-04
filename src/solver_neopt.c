/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double *C = (double*) malloc(N * N * sizeof(double));
	double *temp = (double*) malloc(N * N * sizeof(double));

	/* Calculate A × B and store the result in `temp`. The matrices are not
	   traversed entirely, as A is upper triangular. */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0.0;
			for (int k = i; k < N; k++) {
				sum += A[i * N + k] * B[k * N + j];
			}
			temp[i * N + j] = sum;
		}
	}

	/* Calculate temp x A^t and store the result in C. The matrices are not
	   traversed entirely, as A^t is lower triangular. */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0.0;
			for (int k = j; k < N; k++) {
				sum += temp[i * N + k] * A[j * N + k];
			}
			C[i * N + j] = sum;
		}
	}

	/* Calculate B^t x B^t and store the result in temp. */
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double sum = 0.0;
			for (int k = 0; k < N; k++) {
				sum += B[k * N + i] * B[j * N + k];
			}
			temp[i * N + j] = sum;
		}
	}

	/* Calculate the final matrix by adding C and temp matrices */
	for (int i = 0; i < N * N; i++) {
		C[i] = C[i] + temp[i];
	}

	free(temp);
	return C;
}
