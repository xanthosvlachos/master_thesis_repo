#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"
#include "functions.h"

double* sigmoid(double *x, int size_x, int s_e_point, int fl_par)
{
    double *ones;
    double *expon;
    double *return_value;

    int i = 0;
    MKL_INT len = size_x;
    
    fl_par *= -1;
    
    //Set accuracy mode to enchanced, for better performance
    vmlSetMode(VML_EP);

    //Allocate memmory
    ones = (double *) mkl_malloc(sizeof (double) * len, 64);
    if (ones == NULL) exit(1);
    expon = (double *) mkl_malloc(sizeof (double) * len, 64);
    if (expon == NULL) exit(2);
    return_value = (double *) mkl_malloc(sizeof (double) * len, 64);
    if (return_value == NULL) exit(3);
    
    //Initialize ones array to a CONSTANT VALUE
    for(i=0; i<size_x;i++)
        ones[i] = 1;

    //find better way to add integer (s_e_point) and multiply by integer (fl_param) to each array element
    for(i = 0; i < size_x; i++) {
        x[i] += s_e_point;
        x[i] *= fl_par;
    }

    //Exponential calculation using MKL Vector Mathematical Functions
    vdExp( len, x, expon);

    //find better way to add (+1) to each array element
    for(i = 0; i < size_x; i++) {
        x[i] += 1;
    }

    vdDiv( len, ones, x, return_value);

    printf("\n\n");
    for( i = 0; i < len; i++) {
        printf(" R_V[%d] = %7.9f \n",i,return_value[i]); 
    }
    printf("\n\n");
    printf("LEN = %d\n\n",len);


    mkl_free (ones);
    mkl_free (expon);
    mkl_free (return_value);

    return return_value;
}