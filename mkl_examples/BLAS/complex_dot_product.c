#include "mkl.h"
#include <stdio.h>
#include <stdlib.h>

#define N 5

int main() {
	int n = N, inca = 1, incb = 1, i;
	MKL_Complex16 a[N], b[N], c;
	for(i = 0; i < n; i++) {
		a[i].real = (double)i;
		a[i].imag = (double)i * 2.0;
		b[i].real = (double)(n - 1);
		b[i].imag = (double)i * 2.0;
	}
	
	for(i = 0; i < n; i++)
		printf("A[i].real = %6.2f, A[i].imag = %6.2f, B[i].real = %6.2f, B[i].imag = %6.2f\n", a[i].real, a[i].imag, b[i].real, b[i].imag);

	zdotc(&c, &n, a, &inca, b, &incb);
	printf("The complex dot product is: (%6.2f, %6.2f) \n", c.real, c.imag);

	return 0;
}
