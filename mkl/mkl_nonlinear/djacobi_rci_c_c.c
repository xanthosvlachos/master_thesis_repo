<!DOCTYPE html
  PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<!-- saved from url=(0014)about:internet -->
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
</head>
<body>
<pre>&nbsp;
/*******************************************************************************
* Copyright 2004-2015 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
*   Content : DJACOBI RCI example
*
******************************************************************************** 
*/
/* 
  The program computes the Jacobi matrix of the function on the basis of RCI
  using the central difference.
*/

#include &lt;stdlib.h&gt;
#include &lt;stdio.h&gt;
#include "mkl_service.h"

#include "mkl_rci.h"
#include "mkl_types.h"

int main ()
{
  /* user�s objective function */
  extern void extended_powell (MKL_INT *, MKL_INT *, double *, double *);
  /* n - number of function variables
     m - dimension of function value */
  MKL_INT n = 4, m = 4;
  /* jacobi matrix */
  double *a;
  /* solution vector. contains values x for f(x)
     temporary arrays f1 &amp; f2 that contain f1 = f(x+eps) | f2 = f(x-eps) */
  double *x, *f1, *f2;
  /* precisions for jacobi_matrix calculation */
  double eps = 0.000001;
  /* jacobi-matrix solver handle */
  _JACOBIMATRIX_HANDLE_t handle;
  /* controls of rci cycle */
  MKL_INT successful, rci_request, i;

  if ((a = (double *) mkl_malloc(sizeof (double) * n * m, 64)) == NULL)
    {
      printf ("\n#fail: jacobi matrix allocation failed\n");
    }
  if ((x = (double *) mkl_malloc(sizeof (double) * n, 64)) == NULL)
    {
      printf ("\n#fail: vector X allocation failed\n");
      mkl_free (a);
      return 1;
    }
  if ((f1 = (double *) mkl_malloc(sizeof (double) * n, 64)) == NULL)
    {
      printf ("\n#fail: vector F1 allocation failed\n");
      mkl_free (x);
      mkl_free (a);
      return 1;
    }
  if ((f2 = (double *) mkl_malloc(sizeof (double) * n, 64)) == NULL)
    {
      printf ("\n#fail: vector F2 allocation failed\n");
      mkl_free (f1);
      mkl_free (x);
      mkl_free (a);
      return 1;
    }

  /* set the x values */
  for (i = 0; i &lt; n; i++)
    x[i] = 10.0;
  /* initialize solver (allocate memory, set initial values) */
  if (djacobi_init (&amp;handle, &amp;n, &amp;m, x, a, &amp;eps) != TR_SUCCESS)
    {
      /* if function does not complete successfully then print error message */
      printf ("\n#fail: error in djacobi_init\n");
      fflush (0);
      MKL_Free_Buffers ();
      return 1;
    }
  /* set initial rci cycle variables */
  rci_request = 0;
  successful = 0;
  /* rci cycle */
  while (successful == 0)
    {
      /* call solver */
      if (djacobi_solve (&amp;handle, f1, f2, &amp;rci_request) != TR_SUCCESS)
        {
      /* if function does not complete successfully then print error message */
          printf ("\n#fail: error in djacobi_solve\n");
          fflush (0);
          MKL_Free_Buffers ();
          return 1;
        }
      if (rci_request == 1)
        {
          /* calculate function value f1 = f(x+eps) */
          extended_powell (&amp;m, &amp;n, x, f1);
        }
      else if (rci_request == 2)
        {
          /* calculate function value f1 = f(x-eps) */
          extended_powell (&amp;m, &amp;n, x, f2);
        }
      else if (rci_request == 0)
        /* exit rci cycle */
        successful = 1;
    } /* rci cycle */
  /* free handle memory */
  if (djacobi_delete (&amp;handle) != TR_SUCCESS)
    {
      /* if function does not complete successfully then print error message */
      printf ("\n#fail: error in djacobi_delete\n");
      fflush (0);
      MKL_Free_Buffers ();
      return 1;
    }
  MKL_Free_Buffers ();
  mkl_free (f2);
  mkl_free (f1);
  mkl_free (x);
  mkl_free (a);
  printf ("#pass\n");
  fflush (0);
  return 0;
}

/* 
    routine for extended Powell function calculation
    m in : dimension of function value
    n in : number of function variables
    x in : vector for function calculation
    f out: function value f(x)
*/
void extended_powell (MKL_INT * m, MKL_INT * n, double *x, double *f)
{
  MKL_INT i;
  for (i = 0; i &lt; (*n) / 4; i++)
    {
      f[4 * i] = x[4 * i] + 10.0 * x[4 * i + 1];
      f[4 * i + 1] = 2.2360679774 * (x[4 * i + 2] - x[4 * i + 3]);
      f[4 * i + 2] = (x[4 * i + 1] - 2.0 * x[4 * i + 2]) * 
                     (x[4 * i + 1] - 2.0 * x[4 * i + 2]);
      f[4 * i + 3] = 3.1622776601 * (x[4 * i] - x[4 * i + 3]) * 
                                    (x[4 * i] - x[4 * i + 3]);
    }
  return;
}
&nbsp;</pre> 
</body>

</html>