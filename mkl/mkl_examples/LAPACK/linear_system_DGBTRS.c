
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, double* a, int lda );
extern void print_int_vector( char* desc, int n, int* a );

/* Parameters */
#define N 3
#define NRHS 1
#define LDA N
#define LDB N
#define ZERO 0
#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102


/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, z = ZERO;
	/* Local arrays */
	MKL_INT ipiv[N];
	double a[LDA*N] = {
	   1.00, 1.00, 1.00,
	   0.00, 2.00, 5.00,
	   2.00, 5.00, -1.00
	};
	double b[LDB*NRHS] = {
	    6.00,
	    -4.00,
	    27.00
	};
	/* Executable statements */
	printf( " DGBTRS Example Program Results\n" );
	//Computes the LU factorization of a general m-by-n band matrix.

	
	
	info = LAPACKE_dgbtrf(LAPACK_ROW_MAJOR, n, n, z, z, a, lda, ipiv);
	if( info > 0 ) {
		exit( 1 );
	}
	print_matrix( "LU factored matrix", n, n, a, lda );
	/* Solve the equations A*X = B */
	info = LAPACKE_dgbtrs( LAPACK_ROW_MAJOR, 'N', n, z, z, nrhs, a, lda, ipiv, b, ldb);
	/* Check for the exact singularity */
	if( info > 0 ) {
		exit( 1 );
	}
	
	
	
	/* Print solution */
	print_matrix( "Solution", n, nrhs, b, ldb );
	/* Print pivot indices */
	print_int_vector( "Pivot indices", n, ipiv );
	exit( 0 );
} 

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, double* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, int n, int* a ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
}
