/*********************************************************************
 **
 ** file: lvs_rlm.c
 **
 ** Aim: implement robust linear models to be used for LVS normalization
 **
 ** Stefano Calza <stefano.calza@med.unibs.it> 2010
 ** 
 */


#include "basic_fns.h"
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


#ifdef Win32
# include <fcntl.h>
#endif



SEXP rlm_fit(SEXP X, SEXP Y, SEXP Beta, SEXP max_iter, SEXP initialized,SEXP tolin)
{

  double sigma =0.0, conv, psi_k = 1.345, acc = 0.01;


  int maxit=20,init=0,iter,j,i,rows,cols,converged = 0;

  double *old_beta, *pred, *tmpBeta, *A, *B, *Xwt;

  SEXP ans, NewBeta,  names, Sigma, resids, wts, Xdims, Pred;


  maxit = INTEGER(max_iter)[0];
  init = INTEGER(initialized)[0];
  

  PROTECT(Sigma = allocVector(REALSXP,1));

  /*   List with results */
  PROTECT(ans = allocVector(VECSXP,7)); //24Mar2010: added converged & iteration
  
  Xdims = getAttrib(X, R_DimSymbol);
  
  rows = INTEGER(Xdims)[0];
  cols = INTEGER(Xdims)[1];

  tmpBeta = (double *) R_alloc(cols, sizeof(double));
  old_beta = (double *) R_alloc(cols, sizeof(double));
  pred = (double *) R_alloc(rows, sizeof(double));

  PROTECT(resids = allocVector(REALSXP,rows));
  PROTECT(wts = allocVector(REALSXP,rows));
  PROTECT(NewBeta = allocVector(REALSXP,cols));  
  PROTECT(Pred = allocVector(REALSXP,rows));  

  Xwt = (double *) R_alloc(cols*rows, sizeof(double));
  A = (double *) R_alloc(cols*rows, sizeof(double));
  B = (double *) R_alloc(cols, sizeof(double));


  if (!init)
    {
      
      crossprod(REAL(X),rows,cols,REAL(X),rows,cols,A);
      crossprod(REAL(X),rows,cols,REAL(Y),rows,1,B);

      lvs_dgesv(A, B, tmpBeta, cols, 1, asReal(tolin));

      Memcpy(old_beta,tmpBeta,cols);
    }
  else
    {
      matprod(REAL(X),rows,cols,REAL(Beta),cols,1,pred);
      Memcpy(old_beta,REAL(Beta),cols);
    }

  
  for (i=0; i < rows; i++){
    REAL(resids)[i] = REAL(Y)[i] - pred[i];
  }
    
  for (iter = 0; iter < maxit; iter++){
    
    sigma = lvs_med_abs(REAL(resids),rows)/0.6745;

    for (i=0; i < rows; i++){
      REAL(wts)[i] = lvs_psi_huber2(REAL(resids)[i],psi_k,sigma);
    }
    

    for (i=0; i < rows; i++){
      for (j = 0; j < cols; j++){
	Xwt[j*rows+i] = REAL(wts)[i]*REAL(X)[j*rows+i];
      }
    }

    crossprod(Xwt,rows,cols,REAL(X),rows,cols,A);
    crossprod(Xwt,rows,cols,REAL(Y),rows,1,B);


    lvs_dgesv(A, B, tmpBeta, cols, 1, asReal(tolin));


    matprod(REAL(X),rows,cols,tmpBeta,cols,1,pred);
    
    for (i=0; i < rows; i++){
      REAL(resids)[i] = REAL(Y)[i] - pred[i];
    }    
    
    /*check convergence  based on betas */
    
    conv = lvs_check_conv(old_beta,tmpBeta,cols);    

    if (conv < acc){
      converged = 1;
      iter++; /* to have a record of how many iterations needed */
      break; 
    }

    Memcpy(old_beta, tmpBeta, cols);

  }
  

  REAL(Sigma)[0] = sigma;


  PROTECT(names = allocVector(STRSXP,7));

  SET_STRING_ELT(names, 0, mkChar("Beta"));
  SET_STRING_ELT(names, 1, mkChar("resids"));
  SET_STRING_ELT(names, 2, mkChar("weights"));
  SET_STRING_ELT(names, 3, mkChar("sigma"));
  SET_STRING_ELT(names, 4, mkChar("predicted"));
  SET_STRING_ELT(names, 5, mkChar("converged"));
  SET_STRING_ELT(names, 6, mkChar("iteration"));

  Memcpy(REAL(NewBeta),tmpBeta,cols);
  Memcpy(REAL(Pred),pred,rows);

  SEXP Conv, Iter;
  PROTECT(Conv = allocVector(INTSXP,1));
  PROTECT(Iter = allocVector(INTSXP,1));

  INTEGER(Conv)[0] = converged;
  INTEGER(Iter)[0] = iter;

  
  SET_VECTOR_ELT(ans, 0, NewBeta);
  SET_VECTOR_ELT(ans, 1, resids);
  SET_VECTOR_ELT(ans, 2, wts);
  SET_VECTOR_ELT(ans, 3, Sigma);
  SET_VECTOR_ELT(ans, 4, Pred);
  SET_VECTOR_ELT(ans, 5, Conv);
  SET_VECTOR_ELT(ans, 6, Iter);
  
  setAttrib(ans,R_NamesSymbol,names);

  UNPROTECT(9);

  return ans; /* return a list with Beta, resids, weights, sigma */
  
}


/* Y = residuals */
/* X = design matrix */
void gamma_fit(SEXP X, SEXP Y, int max_iter, double tol, double *mu)
{
  
  int i,j, rows, cols, iter, converged=0;
  double  q3, conv, acc = 0.01;
  double *A, *B, *tmpBeta, *oldBeta, *resids, *abs_resids, *wts, *Xwt, *pred, *Y2, *buffer;
  /* double min_y = 1e-7; */
  double min_y = sqrt(tol);
  
  SEXP Xdims;


  Xdims = getAttrib(X, R_DimSymbol);
  
  rows = INTEGER(Xdims)[0];
  cols = INTEGER(Xdims)[1];

  tmpBeta = (double *) Calloc(cols, double);
  oldBeta = (double *) Calloc(cols, double);
  abs_resids = (double *) Calloc(rows, double);
  resids = (double *) Calloc(rows, double);
  pred = (double *) Calloc(rows, double);
  wts = (double *) Calloc(rows, double);
  Xwt = (double *) Calloc(cols*rows, double);
  A = (double *) Calloc(cols*rows, double);
  B = (double *) Calloc(cols, double);
  Y2 = (double *) Calloc(rows, double);
  buffer = (double *) Calloc(rows, double);

  for(i=0; i < rows; i++)
    {
      /* As input we need residuals^2 */
      Y2[i] = R_pow_di(REAL(Y)[i],2);

      if(Y2[i] <= min_y){
	buffer[i] = log(min_y);
      }
      else{
	buffer[i] = log(Y2[i]);
      }
    }
  
  /* log(y^2) ~ Beta*X */
  crossprod(REAL(X),rows,cols,REAL(X),rows,cols,A);
  crossprod(REAL(X),rows,cols,buffer,rows,1,B);
  
  /*  Starting values  */
  lvs_dgesv(A, B, tmpBeta, cols, 1, tol);

  /* Start iterations */
  for(iter = 1; iter <= max_iter; iter++)
    {

      Memcpy(oldBeta,tmpBeta,cols);
      
      matprod(REAL(X),rows,cols,tmpBeta,cols,1,pred); /* x %*% Beta = fitted values */

      for (i=0; i < rows; i++){
	mu[i] = exp(pred[i]);
	resids[i] = (Y2[i] - mu[i])/mu[i];
	abs_resids[i] = fabs(resids[i]);

	buffer[i] = pred[i] + resids[i];
      }

      q3 = lvs_quartile3(abs_resids,rows);

      for (i=0; i < rows; i++){
	if(abs_resids[i] < q3){
	  wts[i] = 1;}
	else{
	  wts[i] = q3/abs_resids[i];}
      }

      for (i=0; i < rows; i++){
	for (j = 0; j < cols; j++){
	  Xwt[j*rows+i] = wts[i]*REAL(X)[j*rows+i];
	}
      }
	
      crossprod(Xwt,rows,cols,REAL(X),rows,cols,A); /* t(Xw) %*% X */
      crossprod(Xwt,rows,cols,buffer,rows,1,B); /* t(Xw) %*% Y */
      
      /*  Starting values  */
      lvs_dgesv(A, B, tmpBeta, cols, 1, tol);

      conv = lvs_check_conv(oldBeta,tmpBeta,cols);    
    

      if (conv < acc){
	converged = 1; /* I don't output it for the moment */
	break; 
      }
    }

  Free(tmpBeta);
  Free(oldBeta);
  Free(abs_resids);
  Free(resids);
  Free(pred);
  Free(wts);
  Free(Xwt);
  Free(A);
  Free(B);
  Free(buffer);
  Free(Y2);

}



SEXP test_gamma_fit(SEXP X, SEXP Y, int max_iter, double tol)
{
  
  int i,j, rows, cols, iter, converged=0;
  double q3, conv, acc = 0.01;
  double *A, *B, *tmpBeta, *oldBeta, *resids, *abs_resids, *wts, *Xwt, *pred, *Y2, *buffer;
  
  SEXP Xdims, mu;
  double min_y = sqrt(tol);

  Xdims = getAttrib(X, R_DimSymbol);
  
  rows = INTEGER(Xdims)[0];
  cols = INTEGER(Xdims)[1];

  /* PROTECT(Y2 = allocVector(REALSXP,rows));   */
  /* PROTECT(YY = allocVector(REALSXP,rows));   */

  PROTECT(mu = allocVector(REALSXP,rows));  

  tmpBeta = (double *) Calloc(cols, double);
  oldBeta = (double *) Calloc(cols, double);
  abs_resids = (double *) Calloc(rows, double);
  resids = (double *) Calloc(rows, double);
  pred = (double *) Calloc(rows, double);
  wts = (double *) Calloc(rows, double);
  Xwt = (double *) Calloc(cols*rows, double);
  A = (double *) Calloc(cols*rows, double);
  B = (double *) Calloc(cols, double);
  Y2 = (double *) Calloc(rows, double);
  buffer = (double *) Calloc(rows, double);

  for(i=0; i < rows; i++)
    {
      /* As input we need residuals^2 */
      Y2[i] = R_pow_di(REAL(Y)[i],2);

      if(Y2[i] <= min_y){
	buffer[i] = log(min_y);
      }
      else{
	buffer[i] = log(Y2[i]);
      }
    }
  
  /* log(y^2) ~ Beta*X */
  crossprod(REAL(X),rows,cols,REAL(X),rows,cols,A);
  crossprod(REAL(X),rows,cols,buffer,rows,1,B);
  
  /*  Starting values  */
  lvs_dgesv(A, B, tmpBeta, cols, 1, tol);

  /* Start iterations */
  for(iter = 1; iter <= max_iter; iter++)
    {

      Memcpy(oldBeta,tmpBeta,cols);
      
      matprod(REAL(X),rows,cols,tmpBeta,cols,1,pred); /* x %*% Beta = fitted values */

      for (i=0; i < rows; i++){
	REAL(mu)[i] = exp(pred[i]);
	resids[i] = (Y2[i] - REAL(mu)[i])/REAL(mu)[i];
	abs_resids[i] = fabs(resids[i]);

	buffer[i] = pred[i] + resids[i];
      }

      q3 = lvs_quartile3(abs_resids,rows);

      for (i=0; i < rows; i++){
	if(abs_resids[i] < q3){
	  wts[i] = 1;}
	else{
	  wts[i] = q3/abs_resids[i];}
      }

      for (i=0; i < rows; i++){
	for (j = 0; j < cols; j++){
	  Xwt[j*rows+i] = wts[i]*REAL(X)[j*rows+i];
	}
      }
	
      crossprod(Xwt,rows,cols,REAL(X),rows,cols,A); /* t(Xw) %*% X */
      crossprod(Xwt,rows,cols,buffer,rows,1,B); /* t(Xw) %*% Y */
      
      /*  Starting values  */
      lvs_dgesv(A, B, tmpBeta, cols, 1, tol);

      conv = lvs_check_conv(oldBeta,tmpBeta,cols);    
    
      if (conv < acc){
	converged = 1; /* I don't output it for the moment */
	break; 
      }
    }

  Free(tmpBeta);
  Free(oldBeta);
  Free(abs_resids);
  Free(resids);
  Free(pred);
  Free(wts);
  Free(Xwt);
  Free(A);
  Free(B);
  Free(buffer);
  Free(Y2);


  UNPROTECT(1);
  return mu;
}



SEXP joint_fit(SEXP X, SEXP Y, SEXP Beta, SEXP max_iter, SEXP initialized,SEXP tolin)
{

  double sigma =0.0, conv, psi_k = 1.345, acc = 0.01;

  int maxit=0,init=0,iter,j,i,rows,cols,converged = 0;

  double *old_beta, *pred, *tmpBeta, *A, *B, *Xwt, *mu;

  SEXP ans, NewBeta,  names, Sigma, resids, wts, Xdims, Pred;

  PROTECT(Sigma = allocVector(REALSXP,1));

  /*   List with results */
  PROTECT(ans = allocVector(VECSXP,7));
  
  Xdims = getAttrib(X, R_DimSymbol);
  
  rows = INTEGER(Xdims)[0];
  cols = INTEGER(Xdims)[1];


  PROTECT(resids = allocVector(REALSXP,rows));
  PROTECT(wts = allocVector(REALSXP,rows));
  PROTECT(NewBeta = allocVector(REALSXP,cols));  
  PROTECT(Pred = allocVector(REALSXP,rows));

  maxit = INTEGER(max_iter)[0];
  init = INTEGER(initialized)[0];

/*   tmpBeta = (double *) R_alloc(cols, sizeof(double)); */
/*   old_beta = (double *) R_alloc(cols, sizeof(double)); */
/*   pred = (double *) R_alloc(rows, sizeof(double)); */
/*   mu = (double *) R_alloc(rows, sizeof(double)); */
/*   Xwt = (double *) R_alloc(cols*rows, sizeof(double)); */
/*   A = (double *) R_alloc(cols*rows, sizeof(double)); */
/*   B = (double *) R_alloc(cols, sizeof(double)); */
  

  tmpBeta = (double *) Calloc(cols, double);
  old_beta = (double *) Calloc(cols, double);
  pred = (double *) Calloc(rows, double);
  mu = (double *) Calloc(rows, double);
  Xwt = (double *) Calloc(cols*rows, double);
  A = (double *) Calloc(cols*rows, double);
  B = (double *) Calloc(cols, double);


  if (!init){
    crossprod(REAL(X),rows,cols,REAL(X),rows,cols,A);
    crossprod(REAL(X),rows,cols,REAL(Y),rows,1,B);
  
    lvs_dgesv(A, B, tmpBeta, cols, 1, asReal(tolin));
    matprod(REAL(X),rows,cols,tmpBeta,cols,1,pred);
    Memcpy(old_beta,tmpBeta,cols);
  }
  else{
    matprod(REAL(X),rows,cols,REAL(Beta),cols,1,pred);
    Memcpy(old_beta,REAL(Beta),cols);
  }
  

  for (i=0; i < rows; i++){
    REAL(resids)[i] = REAL(Y)[i] - pred[i];
  }
  
  for (iter = 1; iter <= maxit; iter++){
    
    sigma = lvs_med_abs(REAL(resids),rows)/0.6745;

    gamma_fit(X,resids,maxit,asReal(tolin),mu);

    for (i=0; i < rows; i++){
      REAL(wts)[i] = lvs_psi_huber3(REAL(resids)[i],psi_k,sigma,mu[i]);
    }
    

    for (i=0; i < rows; i++){
      for (j = 0; j < cols; j++){
	Xwt[j*rows+i] = REAL(wts)[i]*REAL(X)[j*rows+i];
      }
    }

    crossprod(Xwt,rows,cols,REAL(X),rows,cols,A);
    crossprod(Xwt,rows,cols,REAL(Y),rows,1,B);


    lvs_dgesv(A, B, tmpBeta, cols, 1, asReal(tolin));


    matprod(REAL(X),rows,cols,tmpBeta,cols,1,pred);
    
    for (i=0; i < rows; i++){
      REAL(resids)[i] = REAL(Y)[i] - pred[i];
    }    
    
    /*check convergence  based on betas */
    
    conv = lvs_check_conv(old_beta,tmpBeta,cols);    


    if (conv < acc){
      converged = 1; /* I don't output it for the moment */
      break; 
    }

    Memcpy(old_beta, tmpBeta, cols);

  }
  

/*   for(i=0; i< rows;i++) */
/*     { */
/*     printf("%.3f\n",mu[i]); */
/*     } */

  REAL(Sigma)[0] = sigma;
  
  PROTECT(names = allocVector(STRSXP,7));

  SET_STRING_ELT(names, 0, mkChar("Beta"));
  SET_STRING_ELT(names, 1, mkChar("resids"));
  SET_STRING_ELT(names, 2, mkChar("weights"));
  SET_STRING_ELT(names, 3, mkChar("sigma"));
  SET_STRING_ELT(names, 4, mkChar("predicted"));
  SET_STRING_ELT(names, 5, mkChar("converged"));
  SET_STRING_ELT(names, 6, mkChar("iteration"));

  Memcpy(REAL(NewBeta),tmpBeta,cols);
  Memcpy(REAL(Pred),pred,rows);

  SEXP Conv, Iter;
  PROTECT(Conv = allocVector(INTSXP,1));
  PROTECT(Iter = allocVector(INTSXP,1));

  INTEGER(Conv)[0] = converged;
  INTEGER(Iter)[0] = iter;

  
  SET_VECTOR_ELT(ans, 0, NewBeta);
  SET_VECTOR_ELT(ans, 1, resids);
  SET_VECTOR_ELT(ans, 2, wts);
  SET_VECTOR_ELT(ans, 3, Sigma);
  SET_VECTOR_ELT(ans, 4, Pred);
  SET_VECTOR_ELT(ans, 5, Conv);
  SET_VECTOR_ELT(ans, 6, Iter);

  setAttrib(ans,R_NamesSymbol,names);

  UNPROTECT(9);


  Free(tmpBeta);
  Free(old_beta);
  Free(pred);
  Free(mu);
  Free(Xwt);
  Free(A);
  Free(B);


  return ans; /* return a list with Beta, resids, weights, sigma */
  
}
