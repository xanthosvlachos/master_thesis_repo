#include <stdio.h>
#include <time.h>
#include "mkl_vml.h"
#include "mkl.h"
#include <limits.h>
#include <math.h>

#define LEN 3


void f_wirelength(double* X, double* Y) {
  const MKL_INT len = LEN;
  //const MKL_INT one = 1;
  int i;
  double exp_sum, log_e_sum;
  double wirelength;
  
  double gamma_p[LEN];
  double gamma_n[LEN];

  double Xp[len], Xn[len], Yp[len], Yn[len];
  double e_Xn[len], e_Xp[len], e_Yn[len], e_Yp[len];
  double S_e_Xn, S_e_Xp, S_e_Yn, S_e_Yp;
  

  //Initialize gamma_n and gamma_p to a CONSTANT VALUE
  for(i=0; i<len;i++) {
    gamma_p[i] = 5;
    gamma_n[i] = -5;
  }
  
  //Set accuracy mode to enchanced, for better performance
  vmlSetMode(VML_EP);
  
  //A) Divide each element of my vectors by gamma
  vdDiv( len, X, gamma_p, Xp );
  vdDiv( len, X, gamma_n, Xn );
  vdDiv( len, Y, gamma_p, Yp );
  vdDiv( len, Y, gamma_n, Yn );  
  
  printf("A) Divide each element of my vectors X and Y by gamma and get Xn, Xp, Yn, Yp\n\n");
  printf("    X        Xn        Xp        Y         Yn        Yp    \n");
  printf("===============================================================================\n");
  for(i=0;i<len;i++) {
    printf("% 7.6f % 7.6f % 7.6f % 7.6f % 7.6f % 7.6f\n",X[i],Xn[i],Xp[i],Y[i],Yn[i],Yp[i]);
  }
  printf("\n");
  
  
  //B)Calculate exponential for the elements of each of the above vectors
  vdExp( len, Xp, e_Xp);
  vdExp( len, Xn, e_Xn);
  vdExp( len, Yp, e_Yp);
  vdExp( len, Yn, e_Yn);
  
  printf("B) Calculate exponential for the elements of each of the above vectors\n\n");
  printf("    Xn       e_Xn       Xp       e_Xp        Yn       e_Yn      Yp       e_Yp    \n");
  printf("===============================================================================\n");
  for(i=0;i<len;i++) {
    printf("% 7.6f % 7.6f % 7.6f % 7.6f % 7.6f % 7.6f % 7.6f % 7.6f\n",Xn[i],e_Xn[i],Xp[i],e_Xp[i],Yn[i],e_Yn[i],Yp[i],e_Yp[i]);
  }
  printf("\n");
  
  //C) Sum reduction of the four above vectors
  for(i=0;i<len;i++) {
    S_e_Xp += e_Xp[i];
    S_e_Xn += e_Xn[i];
    S_e_Yp += e_Yp[i];
    S_e_Yn += e_Yn[i];
  }

  printf("C) Sum reduction of the four above vectors\n\n");
  printf("  S_e_Xn    S_e_Xp    S_e_Yn    S_e_Yp    \n");
  printf("===============================================================================\n");
  printf("% 7.6f % 7.6f % 7.6f % 7.6f\n",S_e_Xn,S_e_Xp,S_e_Yn,S_e_Yp);
  printf("\n");
  
  //D) Sum of the exponentials above
  exp_sum = S_e_Xp + S_e_Xn + S_e_Yp + S_e_Yn;
  printf("Exponential of the above sum: %7.6f\n",exp_sum);
  printf("\n");
  
  //E) Log(exp_sum)
  log_e_sum = log(exp_sum);
  printf("Log of the above exp_sum: %7.6f\n",log_e_sum);
  printf("\n");

  //F) Calculate Wirelength
  wirelength = gamma_p[0]*log_e_sum;
  printf("WIRELENGTH: %7.6f\n",wirelength);
  printf("\n");
  
  return;
}

int main()
{
  int i;
  double X[LEN],Y[LEN];
  
  srand(time(NULL));
  
  /*for(i=0;i<LEN;i++) {
    X[i] = i+1;//(double)rand()/INT_MAX;
    Y[i] = i+2;//(double)rand()/INT_MAX;
  }*/
  X[0] = 125;
  Y[0] = 187.5;
  X[1] = 150;
  Y[1] = 200;
  X[2] = 175;
  Y[2] = 212.5;
  
  for(i=0;i<LEN;i++) {
    printf("%7.6f %7.6f\n" , X[i], Y[i]);
  }
  
  f_wirelength(X, Y);
 
  return 0;
}