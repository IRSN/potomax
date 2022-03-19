#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#define NODEBUG

/* ===========================================================================
   AUTHOR Yves Deville <deville.yves@alpestat.com>
   
   GOAL Compute the probability functions related to the one parameter
   Exponential Distribution, possibly including the gradient and the
   Hessian for the density, the cumulative distribution and the
   quantile functions. The implementation in C is faster than a pure R
   implementation. NAs are returned when the shape is negative, which
   can be a desirable behaviour for when unconstrained optimisation is
   used to maximise the log-likelihood.
 
   NOTE This program is part of the 'potomax' R package.

   LICENCE Contact the author for details. The 'potomax' R package is
   still in a early stage of development.
   =========================================================================== */  

/* ==========================================================================
 * Density function
 * ========================================================================== */

SEXP Call_dexp1(SEXP x,             /*  double                          */
		SEXP scale,         /*  double                          */
		SEXP logFlag,       /*  integer                         */
		SEXP derivFlag,     /*  integer                         */ 
		SEXP hessianFlag) { /*  integer                         */
  
  int n, nx, nscale, i, ix, iscale,
    deriv = INTEGER(derivFlag)[0], hessian = INTEGER(hessianFlag)[0];
  
  double z;
  
  SEXP val;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));

  double *rx = REAL(x), *rscale = REAL(scale);

  nx = LENGTH(x);											
  nscale = LENGTH(scale);		
  
  if ((nx == 0) || (nscale == 0)) 			
    return(allocVector(REALSXP, 0));				
  
  n = nx;						       			
  if (n < nscale) n = nscale;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);
  
  if (deriv) {

    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n));
    double *rgrad = REAL(grad);

    /* =====================================================================
       NOTE: although some variables that are used only when the
       hessian is required and hence could be declared only in this
       case, it seems that this does not work. It seems that the
       compilator does not see declarations inside two nested 'if'.
       ===================================================================== */

    // if (hessian) { 

    SEXP hess;
    
    PROTECT(hess = allocVector(REALSXP, n));
    double *rhess = REAL(hess);
    
    // specific auxiliary variables for hessian
    double sigma;
    
    // }
    
    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

    for (i = ix = iscale = 0;  i < n; 
	 ix = (++ix == nx) ? 0 : ix, 	       
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ++i) {
      
      if (ISNA(rx[ix]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
	rgrad[i] = NA_REAL;

	if (hessian) {
	  // row 'sigma'
	  rhess[i] = NA_REAL;
	}
	
      } else {
	
	z = rx[ix] / rscale[iscale];
	sigma = rscale[iscale]; 
 
	// Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, sigma2);
 
	if (z < 0.0) {
	  
	  rval[i] = 0.0;
	  rgrad[i] = 0.0;
	  
	  if (hessian) {
	    rhess[i] = 0.0;
	  }

	} else {
	  
	  rval[i] = -log(sigma) - z;
	  rgrad[i] = -(1.0 - z) / sigma;
	  
	  if (hessian) {
	    rhess[i] =  (1 - 2 * z) / sigma / sigma;
	  }
	}
	 
	if (!INTEGER(logFlag)[0]) {
	  rval[i] = exp(rval[i]);
	  rgrad[i] *= rval[i];
	  
	  if (hessian) {
	    // multiply the hessian by the density and add the
	    // 'tcrossprod' of the gradient
	    rhess[i] *= rval[i];
	    rhess[i] += rgrad[i] * rgrad[i];
	  }
	  
	}
	
      }   /* non-NA case     */
      
    }

    SET_ATTR(val, attrNm, grad);

    if (hessian) {
      SET_STRING_ELT(attrNm, 0, mkChar("hessian"));
      SET_ATTR(val, attrNm, hess);
    } 

    UNPROTECT(6);

    return(val);
    
  } else {

    for (i = ix = iscale = 0;  i < n; 
	 ix = (++ix == nx) ? 0 : ix, 	        
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ++i) {
      
      if (ISNA(rx[ix]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
		
      } else {
	
	z = rx[ix] / rscale[iscale];
	// Rprintf("%d, %d, %6.3f, sigma = %6.3f\n", i, ishape, xi, rscale[iscale]);
       
	if (z < 0.0) {
	    rval[i] = 0.0;	  
	} else {
	  rval[i] = -log(rscale[iscale]) - z;	
	}
	
	if (!INTEGER(logFlag)[0]) {
	  rval[i] = exp(rval[i]);
	}
	
      }   /* non-NA case     */
      
    }
  
    UNPROTECT(3);
    return(val);
    
  }

}

 /* ==========================================================================
  * Distribution function
  * ========================================================================== */

SEXP Call_pexp1(SEXP q,               /*  double                          */
		SEXP scale,           /*  double                          */
		SEXP lowerTailFlag,   /*  integer                         */
		SEXP derivFlag,       /*  integer                         */
		SEXP hessianFlag) {   /*  integer                         */
  
  int n, nq, nscale, i, iq, iscale,
    deriv = INTEGER(derivFlag)[0], hessian = INTEGER(hessianFlag)[0];
  
  double z, S, sigma;
  
  SEXP val;
  
  PROTECT(q = coerceVector(q, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));

  double *rq = REAL(q), *rscale = REAL(scale);

  nq = LENGTH(q);					       						
  nscale = LENGTH(scale);		
  
  if ((nq == 0) || (nscale == 0)) 			
    return(allocVector(REALSXP, 0));				
  
  n = nq;						       						
  if (n < nscale) n = nscale;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);

  if (deriv) {

    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n));
    double *rgrad = REAL(grad);

    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
    
    // if (hessian) { 

    SEXP hess;

    PROTECT(hess = allocVector(REALSXP, n));
    double *rhess = REAL(hess);
    double z1;
    
    // }
   
    for (i = iq = iscale = 0;  i < n; 
	 iq = (++iq == nq) ? 0 : iq, 	       
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ++i) {
      
      // as a general rule, only 3 elements of the hessan can be non
      // zero
      if (hessian) {
	// row 'mu'
	rhess[i] = 0.0;
      }
      
      if (ISNA(rq[iq]) || (rscale[iscale] <= 0.0)) {

	rval[i] = NA_REAL;
	rgrad[i] = NA_REAL;

	if (hessian) {
	  rhess[i] = NA_REAL;
	}
	
      } else {
	
	z = rq[iq] / rscale[iscale];
	sigma = rscale[iscale];
       
	if (z < 0.0) {
	  S = 1.0;
	  rval[i] = S;
	  rgrad[i] = 0.0;
	} else {
	  S = exp(-z);
	  rval[i] = S;
	  rgrad[i] = S * z / sigma;
	}
	
	if (hessian) {
	  z1 = 2.0 - z;
	  rhess[i] = - z * z1 * S / sigma / sigma;
	}
	
	if (INTEGER(lowerTailFlag)[0]) {
	  rval[i] = 1.0 - rval[i];
	  rgrad[i] = -rgrad[i];
	  if (hessian) {
	    rhess[i] = -rhess[i];
	  }
	}
	
      }   /* non-NA case     */
      
    }
  
    SET_ATTR(val, attrNm, grad);

    if (hessian) {
      SET_STRING_ELT(attrNm, 0, mkChar("hessian"));
      SET_ATTR(val, attrNm, hess);
    } 

    UNPROTECT(6);
    return(val);
    
  } else {

    for (i = iq = iscale = 0;  i < n; 
	 iq = (++iq == nq) ? 0 : iq, 	        
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ++i) {
      
      if (ISNA(rq[iq]) || (rscale[iscale] <= 0.0)) {
	rval[i] = NA_REAL;	
      } else {
	z = rq[iq] / rscale[iscale];
	if (z < 0.0) {
	  rval[i] = 1.0;
	}	  
	else {
	  rval[i] = exp(-z);
	}
       
	if (INTEGER(lowerTailFlag)[0]) {
	  rval[i] = 1.0 - rval[i];
	}
	
      }   /* non-NA case     */
      
    }
  
    UNPROTECT(3);
    return(val);
    
  }

}

/* ==========================================================================
 * Quantile function
 * ========================================================================== */

SEXP Call_qexp1(SEXP p,               /*  double                          */
		SEXP scale,           /*  double                          */
		SEXP lowerTailFlag,   /*  integer                         */
		SEXP derivFlag,       /* integer                          */
		SEXP hessianFlag) {   /*  integer                         */
  
  int n, np, nscale, i, ip, iscale,
    deriv = INTEGER(derivFlag)[0], hessian = INTEGER(hessianFlag)[0] ;
  
  double q, lq, sigma, rpi;
  
  SEXP val;
  
  PROTECT(p = coerceVector(p, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));

  double *rp = REAL(p), *rscale = REAL(scale);

  np = LENGTH(p);	    	
  nscale = LENGTH(scale);		
  
  if ((np == 0) || (nscale == 0)) 			
    return(allocVector(REALSXP, 0));				
  
  n = np;		  			
  if (n < nscale) n = nscale;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);

  if (deriv) {

    SEXP grad, attrNm;
    
    PROTECT(grad = allocVector(REALSXP, n));
    double *rgrad = REAL(grad);
    
    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
    
    SEXP hess;
    
    PROTECT(hess = allocVector(REALSXP, n));
    double *rhess = REAL(hess);

    for (i = ip = iscale = 0;  i < n; 
	 ip = (++ip == np) ? 0 : ip, 	       
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ++i) {

      if (ISNA(rp[ip]) || (rscale[iscale] <= 0.0)) {

	rval[i] = NA_REAL;
	rgrad[i] = NA_REAL;
	rhess[i] = NA_REAL;
	
      } else {
	
	if (hessian) {
	  rhess[i] = 0.0;
	}

	rpi = rp[ip];
	
	if (!INTEGER(lowerTailFlag)[0]) {
	  rpi = 1.0 - rpi;
	}
	
	sigma = rscale[iscale];
	q = 1.0 - rpi;
	lq = log(q);
	rval[i] = - sigma * lq;
	rgrad[i] =  -lq;
     
      }   /* non-NA case     */
      
    }     /* loop            */

    SET_ATTR(val, attrNm, grad);
    
    if (hessian) {
      SET_STRING_ELT(attrNm, 0, mkChar("hessian"));
      SET_ATTR(val, attrNm, hess);
    } 

    UNPROTECT(6);
    return(val);
    
  } else {

    for (i = ip = iscale = 0;  i < n; 
	 ip = (++ip == np) ? 0 : ip, 	       
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ++i) {
      
      if (ISNA(rp[ip]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
		
      } else {

	rpi = rp[ip];
	
	if (!INTEGER(lowerTailFlag)[0]) {
	  rpi = 1.0 - rpi;
	}

	q = 1 - rpi;
	rval[i] = -rscale[iscale] * log(q);	  
	
      }   /* non-NA case     */
      
    }     /* loop            */
  
    UNPROTECT(3);
    return(val);
    
  }
}
