#include <stdio.h>
#include "mkl.h"
#include <math.h>
#include "functions.h"

int main(void) {
	int i;
	int size = 1000;
	int s_e_point = 50;
	int fl_par = 50;
	
	double *x;

	 //Allocate memmory
    x = (double *) malloc(sizeof (double) * size);
    if (x == NULL) exit(99);
	
	x[0] = 0;
	for( i = 1; i < size; i++) {
		x[i] = x[i-1] + 0.1; 
	}

	sigmoid(x, size, s_e_point, fl_par);

	return (0);
}


