/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double *C = (double*) malloc(N * N * sizeof(double));
	double *temp = (double*) malloc(N * N * sizeof(double));
	double *B_trans = (double*) malloc(N * N * sizeof(double));

	/* Compute B transposed in order to use it in further calculation
	   and optimize the caching process */
	double *temp_B_trans = B_trans;
	for (int i = 0; i < N; ++i) {
		double *temp2_B_trans = temp_B_trans;
		for (int j = 0; j < N; ++j) {
			*temp2_B_trans = B[j * N + i];
			++temp2_B_trans;
		}

		temp_B_trans += N;
	}

	/* Calculate A × B and store the result in `temp`. The matrices are not
	   traversed entirely, as A is upper triangular. */
	double *temp_A = A;
	double *tmp_temp = temp;
	for (int i = 0; i < N; ++i) {
		double *temp2_A = temp_A;
		double *tmp2_temp = tmp_temp;
		for (int j = 0; j < N; ++j) {
			register double *temp_B_trans = &(B_trans[j * N]) + i;
			register double *temp_A2 = temp2_A + i;

			register double sum = 0.0;
			for (register int k = i; k < N; ++k) {
				sum += *temp_A2 * *temp_B_trans;
				++temp_A2;
				++temp_B_trans;
			}

			*tmp2_temp = sum;
			++tmp2_temp;
		}

		temp_A += N;
		tmp_temp += N;
	}

	/* Calculate temp x A^t and store the result in C. The matrices are not
	   traversed entirely, as A^t is lower triangular. */
	tmp_temp = temp;
	double *tmp_C = C;
	for (int i = 0; i < N; ++i) {
		double *tmp2_C = tmp_C;
		double *temp_A = A;
		for (int j = 0; j < N; ++j) {
			register double *temp2_A = temp_A + j;
			register double *tmp_temp2 = tmp_temp + j;

			register double sum = 0.0;
			for (register int k = j; k < N; ++k) {
				sum += *tmp_temp2 * *temp2_A;
				++tmp_temp2;
				++temp2_A;
			}
	
			*tmp2_C = sum;
			++tmp2_C;

			temp_A += N;
		}

		tmp_temp += N;
		tmp_C += N;
	}

	/* Calculate B^t x B^t and store the result in temp. */
	tmp_temp = temp;
	temp_B_trans = B_trans;
	for (int i = 0; i < N; ++i) {
		double *tmp2_temp = tmp_temp;
		double *tmp_B = B;
		for (int j = 0; j < N; ++j) {
			register double sum = 0.0;

			register double *tmp2_B = tmp_B;
			register double *temp2_B_trans = temp_B_trans;
			for (register int k = 0; k < N; ++k) {
				sum += *temp2_B_trans * *tmp2_B;
				++temp2_B_trans;
				++tmp2_B;
			}

			*tmp2_temp = sum;
			++tmp2_temp;

			tmp_B += N;
		}

		temp_B_trans += N;
		tmp_temp += N;
	}

	/* Calculate the final matrix by adding C and temp matrices */
	for (int i = 0; i < N * N; ++i) {
		C[i] = C[i] + temp[i];
	}

	free(B_trans);
	free(temp);
	
	return C;
}
