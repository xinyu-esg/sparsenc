#include "galois.h"

// perform forward substitution on a matrix to transform it to a upper triangular structure
//static long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT A[][ncolA], GF_ELEMENT B[][ncolB])
long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT *A[], GF_ELEMENT *B[])
{
	//printf("entering foward_substitute()...\n");
	long operations = 0;
	int i, j, k, m, n, p;
	int pivot;
	GF_ELEMENT quotient;

	// transform A into upper triangular structure by row operation
	int boundary = nrow >= ncolA ? ncolA : nrow;

	for (i=0; i<boundary; i++) {
		int has_a_dimension = 1;			// to indicate if this column is all-zeros

		if (A[i][i] == 0) {
			has_a_dimension = 0;
			if (i == nrow-1)
				break;
			else {
				for (pivot=i+1; pivot<nrow; pivot++) {
					if (A[pivot][i] != 0) {
						has_a_dimension = 1;
						break;
					}
				}
			}

			// if this column is an all-zeros column, just skip this column
			if (has_a_dimension == 0)
				continue;
			else {
				// do the swap
				GF_ELEMENT tmp2;
				for (m=0; m<ncolA; m++) {
					tmp2 = A[i][m];
					A[i][m] = A[pivot][m];
					A[pivot][m] = tmp2;
				}

				// swap B accordingly
				for (m=0; m<ncolB; m++) {
					tmp2 = B[i][m];
					B[i][m] = B[pivot][m];
					B[pivot][m] = tmp2;
				}
				for (j=i+1; j<nrow; j++) {
					if (A[j][i] == 0)
						continue;			// the matrix we are dealing here is high likely being sparse, so we use this to avoid unnecessary operations
					quotient = galois_divide(A[j][i], A[i][i], GF_ORDER);
					operations += 1;
					// eliminate the items under row i at col i
					galois_multiply_add_region(A[j], A[i], quotient, ncolA, GF_ORDER);
					operations += (ncolA-i);
					// simultaneously do the same thing on right matrix B
					galois_multiply_add_region(B[j], B[i], quotient, ncolB, GF_ORDER);
					operations += ncolB;
				}
			}
		} else {
			for (j=i+1; j<nrow; j++) {
				if (A[j][i] == 0)
					continue;			// the matrix we are dealing here is high likely being sparse, so we use this to avoid unnecessary operations
				quotient = galois_divide(A[j][i], A[i][i], GF_ORDER);
				operations += 1;
				// eliminate the items under row i at col i
				galois_multiply_add_region(A[j], A[i], quotient, ncolA, GF_ORDER);
				operations += (ncolA-i);
				// simultaneously do the same thing on right matrix B
				galois_multiply_add_region(B[j], B[i], quotient, ncolB, GF_ORDER);
				operations += ncolB;
			}
		}
	}
	return operations;
}
				
// perform back-substitution on full-rank upper trianguler matrix A
long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT *A[], GF_ELEMENT *B[])
{
	//printf("entering back_substitute()...\n");
	long operations = 0;

	// 1, transform the upper triangular matrix A into diagonal.
	int i, j, k, l;
	for (i=ncolA-1; i>=0; i--) {
		// eliminate all items above A[i][i]
		for (j=0; j<i; j++) {
			if (A[j][i] == 0)
				continue;				// the matrix we are dealing here is high likely being sparse, so we use this to avoid unnecessary operations
			GF_ELEMENT quotient = galois_divide(A[j][i], A[i][i], GF_ORDER);
			operations += 1;

			A[j][i] = galois_sub(A[j][i], galois_multiply(A[i][i], quotient, GF_ORDER));
			operations += 1;

			// doing accordingly to B
			galois_multiply_add_region(B[j], B[i], quotient, ncolB, GF_ORDER);
			operations += ncolB;
		}
	}

	// 2, transform matrix A into identity matrix
	for (l=0; l<ncolA; l++) {
		if (A[l][l] != 0) {
			for (k=0; k<ncolB; k++) {
				B[l][k] = galois_divide(B[l][k], A[l][l], GF_ORDER);
				operations += 1;
			}
			A[l][l] = 1;
		}
	}

	return operations;
}

/*
 * Reduce square left-hand side matrix, A, to reduced echelon form (REF), and its right-hand side accordingly:
 *   A x = B
 */
long long matrix_to_REF(int ncolA, int ncolB, GF_ELEMENT *A[], GF_ELEMENT *B[])
{
	long long operations = 0;
	int i, j, k, l;
	GF_ELEMENT quotient;

	// Partially diagonalize the LDM
	int consecutive = 1;
	for (k=ncolA-1; k>=0; k--) {
		if (A[k][k] == 0) {
			consecutive = 0;
			continue;
		}

		// eliminate elements above the nonzero diagonal elements
		for (l=0; l<k; l++) {
			if (A[l][k] == 0)
				continue;
				
			quotient = galois_divide(A[l][k], A[k][k], GF_ORDER);
			operations += 1;
			A[l][k] = 0;
			// 注意后面有的列可能并非为零列(也就是对角元素非零)，这些也要eliminate
			for (int m=k+1; m<ncolA; m++) {
				if (A[m][m] == 0) {
					A[l][m] = galois_add(A[l][m], galois_multiply(A[k][m], quotient, GF_ORDER));
					operations += 1;
				}
			}
			galois_multiply_add_region(B[l], B[k], quotient, ncolB, GF_ORDER);
			operations += ncolB;
		}
	}
	if (consecutive == 1) 
		printf("Class %d is self-decodable.\n", i);
	return operations;
}

