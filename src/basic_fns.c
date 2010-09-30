/*********************************************************************
 **
 ** file: basic_fns.c
 **
 ** Aim: Function used by lvs_rlm.c
 **
 ** Copyright (C) 2007 Stefano Calza 
 ** Mainly adapted from code by Ben M .Bolstad <bolstad@stat.berkeley.edu> 2003
 **
 ** 
 ** All code was written by Ben M .Bolstad <bolstad@stat.berkeley.edu> 2003-2007
 ** 
 ** created on: May 10, 2007
 **
 ** Last modified: May 10, 2007
 **
 ** The aim is to provide a set of functions to fit a 
 ** robust linear models to affy data in order to get arrays & probes variance components. 
 ** It is focused on huber regression. Code is inspired by rlm() method which is
 ** part of the MASS package bundle and from affyPLM set of functions.
 */

#include "basic_fns.h"

/* ok */
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h> /* for dgesv and dgecon*/
#include <R_ext/Applic.h> /* for dgemm */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif


void matprod(double *x, int nrx, int ncx,
	     double *y, int nry, int ncy, double *z)
{
  char *transa = "N", *transb = "N";
  int i,  j, k;
  double one = 1.0, zero = 0.0;
  LDOUBLE sum;
  Rboolean have_na = FALSE;
  
  if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
    /* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
     * The test is only O(n) here
     */
    for (i = 0; i < nrx*ncx; i++)
      if (ISNAN(x[i])) {have_na = TRUE; break;}
    if (!have_na)
      for (i = 0; i < nry*ncy; i++)
	if (ISNAN(y[i])) {have_na = TRUE; break;}
    if (have_na) {
      for (i = 0; i < nrx; i++)
	for (k = 0; k < ncy; k++) {
	  sum = 0.0;
	  for (j = 0; j < ncx; j++)
	    sum += x[i + j * nrx] * y[j + k * nry];
	  z[i + k * nrx] = sum;
	}
    } else
      F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
		      x, &nrx, y, &nry, &zero, z, &nrx);
  } else /* zero-extent operations should return zeroes */
    for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}


void crossprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z)
{
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
			x, &nrx, y, &nry, &zero, z, &ncx);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < ncx*ncy; i++) z[i] = 0;
    }
}

void tcrossprod(double *x, int nrx, int ncx,
		      double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "T";
    double one = 1.0, zero = 0.0;
    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	F77_CALL(dgemm)(transa, transb, &nrx, &nry, &ncx, &one,
			x, &nrx, y, &nry, &zero, z, &nrx);
    } else { /* zero-extent operations should return zeroes */
	int i;
	for(i = 0; i < nrx*nry; i++) z[i] = 0;
    }
}



void lvs_dgesv(double *A, double *Bin, double *B, int n, int p, double tol)
{
  int info, *ipiv;
  double *avals, anorm, rcond, *work;
  
  ipiv = (int *) R_alloc(n, sizeof(int));
  
  Memcpy(B, Bin, n * p);

  avals = (double *) R_alloc(n * n, sizeof(double));
  /* work on a copy of A */
  Memcpy(avals, A, n * n);
  F77_CALL(dgesv)(&n, &p, avals, &n, ipiv, B, &n, &info);
  if (info > 0)
    error("Lapack routine dgesv: system is exactly singular");
  anorm = F77_CALL(dlange)("1", &n, &n, A, &n, (double*) NULL);
  work = (double *) R_alloc(4*n, sizeof(double));
  F77_CALL(dgecon)("1", &n, avals, &n, &anorm, &rcond, work, ipiv, &info);
  if (rcond < tol)
    error("system is computationally singular: reciprocal condition number = %g",
	  rcond);
  
}



/**************************************************************************
 **
 ** double median(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x
 **
 *************************************************************************/

double  lvs_median(double *x, int length){
  int i;
  int half;
  double med;
  double *buffer = Calloc(length,double);
  
  memcpy(buffer,x,length*sizeof(double));

  half = (length + 1)/2;

  rPsort(buffer, length, half-1);
  med = buffer[half-1];
  if (length % 2 == 0){
    rPsort(buffer, length, half);
    med = (med + buffer[half])/2.0;
  }
  
  Free(buffer);
  return med;
}




/***************************************************************
 **
 ** double irls_delta(double *old, double *new, int length)
 **
 ** double *old - previous value of a vector
 ** double *new - new value of a vector
 ** int length - length of vector
 **
 ** this function computes the sum of the difference of two vectors
 ** divides this by the sum squared of the old datavector.
 **
 ** the aim of this function is compute something to test for 
 ** convergence in the iteratively reweighted least squares (IRLS)
 ** 
 **
 ** Author: Ben Bolstad
 **
 **************************************************************/

double lvs_irls_delta(double *old, double *new, int length){
  int i=0;
  double sum = 0.0;
  double sum2 =0.0;
  double divisor=1e-20;

  for (i=0; i < length; i++){
    sum = sum + (old[i] - new[i])*(old[i]-new[i]);
    sum2 = sum2 + old[i]*old[i];
  }
  
  if(sum2 >= divisor){
    divisor = sum2;
  }

  return sqrt(sum/divisor); 
} 


/**********************************************************************************
 **
 ** double med_abs(double *x, int length)
 **
 ** double *x - a vector of data
 ** int length - length of the vector.
 ** 
 ** returns the median of the absolute values.
 **
 ** computes the median of the absolute values of a given vector.
 **
 ** Author: Ben Bolstad
 **********************************************************************************/

double lvs_med_abs(double *x, int length){
  int i;
  double med_abs;
  double *buffer = Calloc(length,double);
  
  for (i = 0; i < length; i++)
    buffer[i] = fabs(x[i]);
  
  med_abs = lvs_median(buffer,length);
  
  Free(buffer);
  return(med_abs);
}


double lvs_psi_huber(double u, double k,int deriv)
{
  
  if (deriv == 0){
    if ( 1 < k/fabs(u)){
      return 1.0;
    } else {
      return  k/fabs(u);
    }
  } else if (deriv == 1){
    if (fabs(u) <= k){
      return 1.0;
    } else {
      return 0.0;
    }
  } else {
    if (fabs(u) <= k){
      return u;
    } else {
      if (u < 0){
	return -k;
      } else {
	return k;
      }
    }
  }
}

double lvs_psi_huber2(double u, double k, double sigma)
{
  
  if(fabs(u) <= k*sigma)
    {
      return(1/(sigma*sigma));
    }
  else
    {
      return(k/fabs(u)/sigma);
    }
  
}

double lvs_psi_huber3(double u, double k, double sigma, double mu)
{
  
  if(fabs(u) <= k*sigma)
    {
      return(1/mu);
    }
  else
    {
      return((k*sigma/(fabs(u)+0.01))/mu);
    }
  
}


double lvs_min(double *x, int length)
  {
    int i;
    double min = x[0];
    
    for (i = 1; i < length; i++)
      {
	if (x[i] > 0 && x[i] < min)
	  {
	    min = x[i];
	  }
      }
    
    return min;
  }
  
double lvs_check_conv(double *old, double *new, int len){
  
  int i;
  double max_b = 0.0;
  double tmp_b = 0.0;

  for (i = 0; i < len; i++)
    {
      tmp_b = fabs((new[i] - old[i])/old[i]);
      max_b = fmax2(max_b,tmp_b);
    }

  return max_b;
  
}

/**************************************************************************
 **
 ** double quartile3(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the 3rd quartile of *x. Mimics R function quantile type=7
 **
 *************************************************************************/

double  lvs_quartile3(double *x, int length){

  int hi,lo;
  double index,h;
  double q3;
  double *buffer = Calloc(length,double);
  
  memcpy(buffer,x,length*sizeof(double));

  index = 1 + (length - 1) * 0.75;
  lo = (int) index;
  hi = lo + 1;

  rPsort(buffer, length, lo - 1);

  q3 = buffer[lo - 1];

  h = (index - (double) lo);
  
  if(h > 0)
    {
      rPsort(buffer, length, hi - 1);
      q3 = (1-h)*q3 + h * buffer[hi - 1];}

  Free(buffer);
  return q3;
}

/* SEXP  test_quartile3(SEXP x){ */

/*   int i,hi,lo,len; */
/*   double index,h; */
/*   double q3; */
/*   SEXP Q3; */

/*   len = length(x); */

/*   double *buffer = Calloc(len,double); */
  
/*   memcpy(buffer,REAL(x),len*sizeof(double)); */


/*   index = 1 + (len - 1) * 0.75; */
/*   lo = (int) index; */
/*   hi = lo + 1; */

/*   rPsort(buffer, len, lo - 1); */

/*   q3 = buffer[lo - 1]; */

/*   h = (index - (double) lo); */
  
/*   if(h > 0) */
/*     { */
/*       rPsort(buffer, len, hi - 1); */
/*       q3 = (1-h)*q3 + h * buffer[hi - 1];} */

/*   Free(buffer); */

/*   PROTECT(Q3 = allocVector(REALSXP,1)); */

/*   REAL(Q3)[0]=q3; */
/*   UNPROTECT(1); */
/*   return Q3; */
/* } */
