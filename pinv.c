#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>

#define SV_THRESH .05
#define DIFF_THRESH 0.
#define MIN(q,p) (((q)<(p))?(q):(p))
#define MAX(q,p) (((q)>(p))?(q):(p))

/*

ordering of input matrix is [r1c1, r2c1, r3c1,...rnc1, r1c2, r2c2...

output matrix "inv[]" must be declared to be [n][m]

*/

/* DGESDD prototype */

/*extern void dgesdd( char* jobz, int* m, int* n, double* a,
                int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
                double* work, int* lwork, int* iwork, int* info );
*/

int pinv(int m,int n,double* a[],double* inv[]){
  int info, lwork;
  double wkopt;
  int *iwork;
  double* work;
  int lda=m;
  int ldu=m;
  int ldvt=n;
  FILE *svfl;

  double t1,t2;

  double *aa;
  int jr,jc;
  int sixfour=64;

  /*  mkl_set_dynamic(0); */
  /* mkl_set_num_threads(6); */

  aa=(double*)mkl_calloc(lda*n,sizeof(double),sixfour);
  for( jr=0; jr<lda; jr++)for(jc=0; jc<n; jc++)aa[jr+jc*lda]=a[jr][jc];

  double *sv;
  sv=(double*)mkl_calloc(MIN(m,n),sizeof(double),sixfour);
  double *u;
  u=(double*)mkl_calloc(ldu*m,sizeof(double),sixfour);
  double *vt;
  vt=(double*)mkl_calloc(ldvt*n,sizeof(double),sixfour);


  /* Query and allocate the optimal workspace */
  lwork = -1;
  iwork = (int*)mkl_malloc( 8*MIN(m,n)*sizeof(int),sixfour);

  dgesdd("A", &m, &n, aa, &lda, sv,u,&ldu,vt,&ldvt,&wkopt,&lwork,iwork,&info);

  t1=dsecnd();
  /* calculate SVD */
  lwork=(int)wkopt;
  work = (double*)mkl_malloc(lwork*sizeof(double),sixfour);
  dgesdd("A", &m, &n, aa, &lda, sv,u,&ldu,vt,&ldvt,work,&lwork,iwork,&info);
  t2=dsecnd();

  fprintf(stderr,"\n\ntime for dgesdd: %f\n",(float)((t2-t1)));

  svfl=fopen("singVals","w");
  fprintf(svfl,"%d\n",MIN(m,n));
  for(jr=0; jr<MIN(m,n); jr++)fprintf(svfl," %8.4lf ",sv[jr]);
  fclose(svfl);
 
  /* Check for convergence */
  if( info > 0 ) {
    fprintf(stderr,"The algorithm computing SVD failed to converge. info=%d;\n",info );
    fprintf(stderr,"the least squares solution could not be computed.\n" );

    mkl_free((double*)work);
    mkl_free((int*)iwork);
    mkl_free((double*)aa);
    mkl_free((double*)sv);
    mkl_free((double*)u);
    mkl_free((double*)vt);
    return( -1 );
  }


  /* Determine SV threshold */
  int dif_sp=2; /* spacing between points for smoothed derivative */
  int npts=MIN(m,n);

  double avg_diff=0;
  double *diffs;
  int nsig=0;
  diffs=calloc(npts,sizeof(double));
  for( jc=dif_sp; jc<npts-dif_sp; jc++){
    diffs[jc]=fabs(1/sv[jc+dif_sp]-1/sv[jc-dif_sp])/(double)(2*dif_sp);
    avg_diff+=diffs[jc];
  }

  
  avg_diff/=(double)jc;
  for( jc=npts/2; jc<npts; jc++ ){
    if(log(diffs[jc]/avg_diff) > (double)DIFF_THRESH) break;
    nsig=jc;
    /*    nsig=jc*.98; */ /*sometimes needed...like for the potential mapping */
  }
  free(diffs);

  /*
  double diffs;
  double av_diff=0;
  int nsig=0;
  for(jc=npts/2; jc<npts-1; jc++){
    av_diff+=fabs(1/sv[jc-1]-1/sv[jc+1]);
    if(log(fabs(1/sv[jc-1]-1/sv[jc+1])/(2*av_diff)) > 0) break;
    nsig=jc;
  }
  */
  /*  int nsig=0;
  for(jc=0; jc<MIN(m,n); jc++){
    if( sv[jc] < sv[0]*(double)SV_THRESH) break;
    nsig=jc;
  }
  */
  fprintf(stderr,"NSIG: %d\n",nsig);
  if( nsig <= 0 ){
      mkl_free((double*)work);
      mkl_free((int*)iwork);
      mkl_free((double*)aa);
      mkl_free((double*)sv);
      mkl_free((double*)u);
      mkl_free((double*)vt);
      return( -1);
  }

  double temp;
  int i1=1;
  for( jr=0; jr<MIN(nsig,ldu); jr++)if(sv[jr]!=0){
    temp=(double)(1/sv[jr]);
    dscal(&m,&temp,(u+jr*ldu),&i1);
  }
  for( jr=nsig; jr<ldu; jr++)if(sv[jr]!=0){
    temp=(double)(0.);
    dscal(&m,&temp,(u+jr*ldu),&i1);
  }

  t1=dsecnd();
  double *C;
  C=(double *)mkl_malloc(m*n*sizeof(double),sixfour);


  double one=1;
  double zero=0;
  char *notrans="N";
  char *trans="T";

  dgemm(trans,trans,&n,&ldu,&ldvt,&one,vt,&n,u,&m,&zero,C,&n);

  t2=dsecnd();
  fprintf(stderr,"time for cblas_dgemm: %f\n\n\n",(float)((t2-t1)));

  for( jr=0; jr<n; jr++)for(jc=0; jc<m; jc++)inv[jr][jc]=C[jr+jc*n];

  mkl_free((double*)work);
  mkl_free((int*)iwork);
  mkl_free((double*)aa);
  mkl_free((double*)sv);
  mkl_free((double*)u);
  mkl_free((double*)vt);
  mkl_free((double*)C);

  return(0);
}

