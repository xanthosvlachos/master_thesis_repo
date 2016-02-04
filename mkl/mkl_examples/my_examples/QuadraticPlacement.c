#include <stdio.h>
// #include "mkl.h"
#include "mkl_lapacke.h"

#define NR_ELEMENTS 3
#define NR_B_COEFF 2
#define N_ARRAY_ELEMENTS(array) (sizeof(array)/sizeof((array)[0]))	// ONLY for statically-declared arrays (But we will make sure that that's the case here)
#define NR_ARRAY_ELEMENTS(array) ((&array)[1] - array)				// ONLY for statically-declared arrays (But we will make sure that that's the case here)

/* Auxiliary routines prototypes */
extern void print_matrix(char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda);
extern void print_int_vector(char* desc, MKL_INT n, MKL_INT* a);

/* -------------------------------------------------------------------
 * The LAPACK routines are prefixed by a letter defining the data type:
 *
 * Letter	Description
 * s		Indicates single precision: REAL, or REAL*4 in Fortran, float in C.
 * d		Indicates double precision: DOUBLE PRECISION, or REAL*8 in Fortran, double in C.
 * c		Indicates single precision complex: COMPLEX, or COMPLEX*8 in Fortran, single precision complex data type[1] in C.
 * z		Indicates double precision complex: DOUBLE COMPLEX, or COMPLEX*16 in Fortran, double precision complex data type[1] in C.
 * i		Indicates integer: INTEGER in Fortran, int in C.
 * ?		Indicates a place holder for a single letter in several routine names. For example, ?getrf refers to the routines sgetrf, dgetrf, cgetrf, and zgetrf.
 * ------------------------------------------------------------------
 */

// int main(int argc, char *argv[]) {
int main(void) {
	int i;
	MKL_INT ipiv[NR_ELEMENTS]; // Integer Pivot Indices // Whatever that may mean
	MKL_INT info;
	// Variables right here
	// int CoreWidth = 100;
	// int CoreHeight = 100;

	int FixPoints[]		= {100, 175, 200, 225}; //, 145, 250, 160, 150};
	int NetList[]		= {1, 2, 2, 3}; //, 1, 3, 1, 4, 3, 4, 1, 5, 2, 5};
	int FixPNetList[]	= {1, 1, 3, 2}; //, 2, 3, 3, 4, 4, 1, 4, 2, 4, 3, 4, 4, 5, 3};

	// Count Element Wires
	int ElementConnections[NR_ELEMENTS] = {[0 ... (NR_ELEMENTS - 1)] = 0};
	for (i = 0; i < NR_ARRAY_ELEMENTS(FixPNetList); i += 2) {
		ElementConnections[FixPNetList[i] - 1]++;
		printf(" %d, ",ElementConnections[i]);
	}
	printf("\n\n\n");
	for (i = 0; i < NR_ARRAY_ELEMENTS(NetList); i++) {
		ElementConnections[NetList[i] - 1]++;
		printf(" %d, ",ElementConnections[i]);
	}
	printf("\n\n\n");

	double matrixA[NR_ELEMENTS * NR_ELEMENTS] = {[0 ... (NR_ELEMENTS * NR_ELEMENTS - 1)] = 0}; // also {0}; is valid
	for (i = 0; i < NR_ARRAY_ELEMENTS(NetList); i += 2) {
		matrixA[(NetList[i] - 1)* NR_ELEMENTS + NetList[i + 1] - 1] = -1;
		matrixA[(NetList[i + 1] - 1)* NR_ELEMENTS + NetList[i] - 1] = -1;
	}

	for (i = 0; i < NR_ELEMENTS; i++) {
		matrixA[i * NR_ELEMENTS + i] = ElementConnections[i];
	}

	print_matrix("Matrix A", NR_ELEMENTS, NR_ELEMENTS, matrixA, NR_ELEMENTS);

	double vectorB[NR_ELEMENTS * 2] = {[0 ... (NR_ELEMENTS * 2 - 1)] = 0}; // also {0}; is valid
	for (i = 0; i < NR_ARRAY_ELEMENTS(FixPNetList); i += 2) {
		vectorB[(FixPNetList[i] - 1) * NR_B_COEFF + 0] += FixPoints[FixPNetList[i] - 1];
		vectorB[(FixPNetList[i] - 1) * NR_B_COEFF + 1] += FixPoints[FixPNetList[i]];
	}

	//print_matrix("Vector B", NR_ELEMENTS, NR_B_COEFF, vectorB, NR_B_COEFF);
	printf("\nLAPACKE_dgesv(row-major, high-level): Quadratic Placement Results\n");
	/* Paradime call
	 * info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,	// You can call it LAPACK_COLUMN_MAJOR, but C-logic ain't
	 *						n,					// *extrapolated* Number of Rows
	 *						nrhs,				// CBLAS level 3 function: Matrix x Matrix <-- Number of concated vectors
	 *						a,					// Pointer to A matrix
	 *						lda,				// *extrapolated* Number of Columns -- possibly it is
	 *						ipiv,				// *--Something--*
	 *						b,					// Pointer to B vector - matrix
	 *						ldb					// *extrapolated* Number of Columns
	 *						);
	 */
	info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, NR_ELEMENTS, NR_B_COEFF, matrixA, NR_ELEMENTS, ipiv, vectorB, NR_B_COEFF);

	/* Check for the exact singularity */
	if (info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
		printf("the solution could not be computed.\n");
		return (1);
	}

	/* Print solution */
	//print_matrix("Solution", NR_ELEMENTS, NR_B_COEFF, vectorB, NR_B_COEFF);

	/* Print details of LU factorization */
	//print_matrix("Details of LU factorization", NR_ELEMENTS, NR_ELEMENTS, matrixA, NR_ELEMENTS);

	/* Print pivot indices */
	//print_int_vector("Pivot indices", NR_ELEMENTS, ipiv);
	return (0);
}

/* Auxiliary routine: printing a matrix */
void print_matrix(char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda) {
	MKL_INT i, j;
	printf("\n %s\n", desc);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) printf(" %6.2f", a[i * lda + j]);
		printf("\n");
	}
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector(char* desc, MKL_INT n, MKL_INT* a) {
	MKL_INT j;
	printf("\n %s\n", desc);
	for (j = 0; j < n; j++) printf(" %6i", a[j]);
	printf("\n");
}
