#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>

//    {3.55,  2.90, -5.44,  9.99,  7.99,  4.44}

//#define M 6
//#define N 4
#define M 5000
#define N 5000
#define LDA M
#define LDU M
#define LDVT N


extern int pinv(int m,int n,double** a,double** inv);

int main(){
  /*
  double A[M*N]= {
            8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
            9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
            9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
            5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
            3.16,  7.98,  3.01,  5.80,  4.27, -5.31 
  };
  */
  /*
  double A[M][N]= {
    {8.79,  6.11, -9.15,  9.57, -3.49,  9.84},
    {9.93,  6.91, -7.93,  1.64,  4.02,  0.15},
    {9.83,  5.04,  4.86,  8.83,  9.80, -8.99},
    {5.45, -0.27,  4.85,  0.74, 10.00, -6.02},
    {3.16,  7.98,  3.01,  5.80,  4.27, -5.31}
  };
  */
  /*
   double A[N][M]= {
    {7.52, -0.76,  5.13, -4.75,  1.33, -2.40},
    {-1.10,  0.62,  6.62,  8.52,  4.91, -6.77},
    {-7.95,  9.34, -5.66,  5.75, -5.49,  2.34},
    {1.08, -7.10,  0.87,  5.30, -3.52,  3.95}
  };
  */
  /*
    double A[N][M]= {
    { 1, 0, 0, 0, 2},
    { 0, 0, 3, 0, 0},
    { 0, 0, 0, 0, 0},
    { 0, 4, 0, 0, 0}
  };
  */
  double **a;
  double **inv;
  int m,n,jr,jc;
  int stat;

  a=malloc(M*sizeof(double*));
  for( jr=0; jr<M; jr++)a[jr]=(double*)malloc(N*sizeof(double));

  //for( jr=0; jr<M; jr++)for( jc=0; jc<N; jc++)a[jr][jc]=A[jc][jr];
  srand(time(NULL));
  for( jr=0; jr<M; jr++)for( jc=0; jc<N; jc++)a[jr][jc]=(double)rand()/(double)RAND_MAX;

  inv=malloc(N*sizeof(double*));
  for( jr=0; jr<N; jr++)inv[jr]=(double*)malloc(M*sizeof(double));

  stat=pinv((int)M,(int)N,a,inv);

  /*  for( jr=0; jr<N; jr++ ){
    for( jc=0; jc<M; jc++)printf("%9.4f",inv[jr][jc]);
    printf("\n");
  }

    printf("\n");
    printf("\n");

  int jr1,jc1;
  double **mprod;
  mprod=malloc(N*sizeof(double*));
  for( jr=0; jr<N; jr++)mprod[jr]=(double*)malloc(N*sizeof(double));

  for( jr=0; jr<N; jr++)for(jc=0; jc<N; jc++){
      mprod[jr][jc]=0;
      for(jr1=0; jr1<M; jr1++)mprod[jr][jc]+=inv[jr][jr1]*a[jr1][jc];
    }

  for( jr=0; jr<N; jr++ ){
    for( jc=0; jc<N; jc++)printf("%9.4f",mprod[jr][jc]);
    printf("\n");
  }
  */
}
