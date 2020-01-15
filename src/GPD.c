#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#define NODEBUG

/* ===========================================================================
   AUTHOR Yves Deville <deville.yves@alpestat.com>
   
   GOAL Compute the probability functions related to the two-parameter
   Generalized Pareto Distribution (GPD), possibly including the
   gradient and the Hessian for the density, the cumulative
   distribution and the quantile functions. The implementation in C is
   faster than a pure R implementation. Unlike some other
   implementations in existing R packages, NAs are returned when the
   shape is negative, which can be a desirable behaviour for when
   unconstrained optimisation is used to maximise the log-likelihood.
 
   NOTE This program is part of the 'potomax' R package.

   LICENCE Contact the author for details. The 'potomax' R package is
   still in a early stage of development.
   =========================================================================== */  


/* ==========================================================================
 * Density function
 * ========================================================================== */

SEXP Call_dGPD2(SEXP x,             /*  double                          */
		SEXP scale,         /*  double                          */
		SEXP shape,         /*  double                          */
		SEXP logFlag,       /*  integer                         */
		SEXP derivFlag,     /*  integer                         */ 
		SEXP hessianFlag) { /*  integer                         */
  
  int n, nx, nscale, nshape, i, ix, iscale, ishape,
    deriv = INTEGER(derivFlag)[0], hessian = INTEGER(hessianFlag)[0];
  
  double eps = 1e-6, z, V, xi;
  
  SEXP val;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));
  PROTECT(shape = coerceVector(shape, REALSXP));

  double *rx = REAL(x), *rscale = REAL(scale), 
    *rshape = REAL(shape);

  nx = LENGTH(x);											
  nscale = LENGTH(scale);		
  nshape = LENGTH(shape);
  
  if ((nx == 0) || (nscale == 0) || (nshape == 0)) 			
    return(allocVector(REALSXP, 0));				
  
  n = nx;						       			
  if (n < nscale) n = nscale;
  if (n < nshape) n = nshape;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);
  
  if (deriv) {

    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n * 2));
    double *rgrad = REAL(grad);

    /* =====================================================================
       NOTE: although some variables that are used only when the
       hessian is required and hence could be declared only in this
       case, it seems that this does not work. It seems that the
       compilator does not see declarations inside two nested 'if'.
       ===================================================================== */

    // if (hessian) { 

      SEXP hess;

      PROTECT(hess = allocVector(REALSXP, n * 2 * 2));
      double *rhess = REAL(hess);

      // specific auxiliary variables for hessian
      int j, k;
      double A, B, B1, sigma, sigma2;

      // }
    
    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

    for (i = ix = iscale = ishape = 0;  i < n; 
	 ix = (++ix == nx) ? 0 : ix, 	       
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rx[ix]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
	
	rgrad[i] = NA_REAL;
	rgrad[i + n] = NA_REAL;

	if (hessian) {
	  // row 'sigma'
	  rhess[i] = NA_REAL;
	  rhess[i + n] = NA_REAL;
	  rhess[i + 2 * n] = NA_REAL;
	  // row 'xi'
	  rhess[i + 3 * n] = NA_REAL;
	}
	
      } else {
	
	z = rx[ix] / rscale[iscale];
	xi = rshape[ishape];
	sigma = rscale[iscale]; 
	sigma2 = sigma * sigma;

	// Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, sigma2);
 
	if (fabs(xi) < eps) {

	  if (z < 0.0) {
	    
	    rgrad[i] = 0.0;
	    rgrad[i + n] = 0.0;
	    
	    if (hessian) {
	      
	      // row 'sigma'
	      rhess[i] = 0.0;
	      rhess[i + n] = 0.0;
	      rhess[i + 2 * n] = 0.0;
	      // row 'xi'
	      rhess[i + 3 * n] = 0.0;
	      
	    }

	  } else {
	    
	    rval[i] = -log(sigma) - z;
	    
	    rgrad[i] = -(1.0 - (xi + 1.0) * z) / sigma;
	    rgrad[i + n] = 0.5 * z * (z - 2.0);
	    
	    if (hessian) {

	      // row 'sigma'
	      rhess[i] =  (z + 1) * (z + 1) / sigma / sigma;
	      rhess[i + n] =  - z * (z - 1.0) / sigma;
	      rhess[i + 2 * n] = rhess[i + n];
	      // row 'xi'
	      rhess[i + 3 * n] = z * z * (-2.0 * z / 3.0 + 1);
	      
	    }
	  }
	  
	} else {
	  
	  V = 1.0 + xi * z;

	  // Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, V);
	  
	  if (V > 0.0) {
	
	    A = log(V);
	    B = z / V; 
	    rval[i] = -log(sigma) - (1.0 / xi + 1.0) * A;
	      
	    B1 = (xi + 1.0) * B;
       
	    rgrad[i] = -(1.0 - B1) / sigma;
	    rgrad[i + n] = (A / xi - B1) / xi;
	    
	    if (hessian) {

	      // row 'sigma'
	      rhess[i] = (1.0 - B1 * (2.0 - B * xi)) / sigma2 ;
	      rhess[i + n] =  B * (1.0 - B1) / sigma;
	      rhess[i + 2 * n] = rhess[i + n];
	      // row 'xi'
	      rhess[i + 3 * n] = (-2.0 * A + xi * B * (2.0 + xi * B1)) / xi / xi / xi;

	    }
	    
	  } else {

	    rval[i] = R_NegInf;
	    rgrad[i] = 0.0;
	    rgrad[i + n] = 0.0;

	    if (hessian) {
	      
	      // row 'sigma'
	      rhess[i] = 0.0;
	      rhess[i + n] = 0.0;
	      rhess[i + 2 * n] = 0.0;
	      // row 'xi'
	      rhess[i + 3 * n] = 0.0;
	     
	    }

	  }

	} /* non-exponential case */
	
	if (!INTEGER(logFlag)[0]) {
	  rval[i] = exp(rval[i]);
	  rgrad[i] *= rval[i];
	  rgrad[i + n] *= rval[i];

	  if (hessian) {
	    // multiply the hessian by the density and add the
	    // 'tcrossprod' of the gradient
	    for (j = 0; j < 2; j++) {
	      for (k = 0; k < 2; k++) {
		rhess[i + (j + k * 2) * n] *= rval[i];
		rhess[i + (j + k * 2) * n] += rgrad[i + j * n] * rgrad[i + k * n];
	      }
	    }
	  }
	  
	}
	
      }   /* non-NA case     */
      
    }

    SET_ATTR(val, attrNm, grad);

    if (hessian) {
      SET_STRING_ELT(attrNm, 0, mkChar("hessian"));
      SET_ATTR(val, attrNm, hess);
    } 

    UNPROTECT(7);

    return(val);
    
  } else {

    for (i = ix = iscale = ishape = 0;  i < n; 
	 ix = (++ix == nx) ? 0 : ix, 	        
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rx[ix]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
		
      } else {
	
	z = rx[ix] / rscale[iscale];
	xi = rshape[ishape];
	// Rprintf("%d, %d, %6.3f, sigma = %6.3f\n", i, ishape, xi, rscale[iscale]);
	
	if (fabs(xi) < eps) {
 
	  if (z < 0.0) {
	    rval[i] = 0.0;	  
	  } else {
	    rval[i] = -log(rscale[iscale]) - z;	
	  }
	  
	} else {
	  
	  if (z < 0.0) {
	    rval[i] = 0.0;
	  } else {
	    V = 1.0 + xi * z;
	    if (V > 0.0) {
	      rval[i] = -log(rscale[iscale]) - (1.0 / xi + 1.0) * log(V);
	    } else {
	      rval[i] = R_NegInf;
	    }
	  }

	  // Rprintf("V =  %6.3f, val = %6.3f\n", V, rval[i]);
	  
	} /* non-expoential case */
	
	if (!INTEGER(logFlag)[0]) {
	  rval[i] = exp(rval[i]);
	}

      }   /* non-NA case     */
      
    }
  
    UNPROTECT(4);
    return(val);
    
  }

}

 /* ==========================================================================
  * Distribution function
  * ========================================================================== */

SEXP Call_pGPD2(SEXP q,               /*  double                          */
		SEXP scale,           /*  double                          */
		SEXP shape,           /*  double                          */
		SEXP lowerTailFlag,   /*  integer                         */
		SEXP derivFlag,       /*  integer                         */
		SEXP hessianFlag) {   /*  integer                         */
 
  int n, nq, nscale, nshape, i, iq, iscale, ishape,
    deriv = INTEGER(derivFlag)[0], hessian = INTEGER(hessianFlag)[0];
  
  double eps = 1e-6, z, V, A, B, S, sigma, xi;
  
  SEXP val;
  
  PROTECT(q = coerceVector(q, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));
  PROTECT(shape = coerceVector(shape, REALSXP));

  double *rq = REAL(q), *rscale = REAL(scale), 
    *rshape = REAL(shape);

  nq = LENGTH(q);					       						
  nscale = LENGTH(scale);		
  nshape = LENGTH(shape);
  
  if ((nq == 0) || (nscale == 0) || (nshape == 0)) 			
    return(allocVector(REALSXP, 0));				
  
  n = nq;						       						
  if (n < nscale) n = nscale;
  if (n < nshape) n = nshape;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);

  if (deriv) {

    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n * 2));
    double *rgrad = REAL(grad);

    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
    
    // if (hessian) { 

    SEXP hess;

    PROTECT(hess = allocVector(REALSXP, n * 2 * 2));
    double *rhess = REAL(hess);
    double z1;
    
    // }
   
    for (i = iq = iscale = ishape = 0;  i < n; 
	 iq = (++iq == nq) ? 0 : iq, 	       
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {


      // as a general rule, only 3 elements of the hessan can be non
      // zero
      if (hessian) {
	// row 'mu'
	rhess[i] = 0.0;
	rhess[i + n] = 0.0;
	rhess[i + 2 * n] = 0.0;
	rhess[i + 3 * n] = 0.0;
      }
      
      if (ISNA(rq[iq]) || (rscale[iscale] <= 0.0)) {

	rval[i] = NA_REAL;
	rgrad[i] = NA_REAL;
	rgrad[i + n] = NA_REAL;

	if (hessian) {
	  // row 'sigma'
	  rhess[i ] = NA_REAL;
	  rhess[i + 2 * n] = NA_REAL;
	  // row 'xi'
	  rhess[i + n] = NA_REAL;
	  rhess[i + 3 * n] = NA_REAL;
	}
	
      } else {
	
	z = rq[iq] / rscale[iscale];
	xi = rshape[ishape];
	sigma = rscale[iscale];
	
	if (fabs(xi) < eps) {
   
	  if (z < 0.0) {
	    rval[i] = 1.0;
	    rgrad[i] = 0.0;
	    rgrad[i + n] = 0.0;
	  } else {
	    S = exp(-z);
	    rval[i] = S;
	    rgrad[i] = S * z / sigma;
	    rgrad[i + n] = S * z * z / 2;
	  }

	  // Rprintf("%4d, %4d, %6.3f, %7.2f, %6.3f, %6.3f", i, ishape, xi, z, Z, rscale[iscale]);
	  // Rprintf(" grad: %6.3f, %6.3f, %6.3f\n", rgrad[i], rgrad[i + n], rgrad[i + 2 * n]);

	  if (hessian) {
	    z1 = 2.0 - z;
	    // row 'sigma'
	    rhess[i] = - z * z1 * S / sigma / sigma;
	    rhess[i + 2 * n] = - z * z * z1 * S / 2.0 / sigma; 
	    // row 'xi'
	    rhess[i + n] = rhess[i + 2 * n];
	    rhess[i + 3 * n] = - z * z * z * (8.0 - 3.0 * z) * S / 12.0;
	  }
	 
	} else {
	  
	  if (z < 0.0) {
	    rval[i] = 1.0;
	    rgrad[i] = 0.0;
	    rgrad[i + n] = 0.0;	    
	  } else {
	    
	    V = 1.0 + xi * z;
	    // Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, V);
	    
	    if (V > 0.0) {
	      
	      A = log(V);
	      B = z / V;
	      S = pow(V, - 1.0 / xi);
	      rval[i] = S;
	      
	      rgrad[i] =  S * B / sigma;
	      rgrad[i + n] = S * (A / xi - B) / xi;

	      if (hessian) {
		// row 'sigma'
		rhess[i] = - B * (2.0 - xi * B) * S / sigma / sigma +  rgrad[i] *  rgrad[i] / S;
		rhess[i + 2 * n] = - B * B * S / sigma +  rgrad[i] * rgrad[i + n] / S;
		// row 'xi'
		rhess[i + n] = rhess[i + 2 * n];
		z1 = (- 2.0 * A / xi / xi + B * (2.0 / xi + B)) / xi;
		rhess[i + 3 * n] = z1 * S + rgrad[i + n] *  rgrad[i + n] / S;
	      }
	      
	    } else {
	      
	      if (xi > 0.0) {
		rval[i] = 1.0;
	      } else {
		rval[i] = 0.0;
	      }
	      rgrad[i] = 0.0;
	      rgrad[i + n] = 0.0;   

	      if (hessian) {
		// row 'sigma'
		rhess[i] = 0.0;
		rhess[i + 2 * n] = 0.0;
		// row 'xi'
		rhess[i + n] = 0.0;
		rhess[i + 3 * n] = 0.0;
	      }	      
	    }

	  }

	} /* non-exponential case */
	
	if (INTEGER(lowerTailFlag)[0]) {
	  rval[i] = 1.0 - rval[i];
	  rgrad[i] = -rgrad[i];
	  rgrad[i + n] = -rgrad[i + n];
	  if (hessian) {
	    rhess[i] = -rhess[i];
	    rhess[i + n] = -rhess[i + n];
	    rhess[i + 2 * n] = -rhess[i + 2 * n];
	    rhess[i + 3 * n] = - rhess[i + 3 * n];
	  }
	}
	
      }   /* non-NA case     */
      
    }
  
    SET_ATTR(val, attrNm, grad);

    if (hessian) {
      SET_STRING_ELT(attrNm, 0, mkChar("hessian"));
      SET_ATTR(val, attrNm, hess);
    } 

    
    UNPROTECT(7);
    return(val);
    
  } else {

    for (i = iq = iscale = ishape = 0;  i < n; 
	 iq = (++iq == nq) ? 0 : iq, 	        
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rq[iq]) || (rscale[iscale] <= 0.0)) {
	rval[i] = NA_REAL;	
      } else {
	z = rq[iq] / rscale[iscale];
	xi = rshape[ishape];
	if (fabs(xi) < eps) {
	  if (z < 0.0) {
	    rval[i] = 1.0;
	  }	  
	  else {
	    rval[i] = exp(-z);
	  }
	} else {
	  if (z < 0.0) {
	    rval[i] = 1.0;
	  }
	  else {
	    V = 1.0 + xi * z;
	    if (V > 0.0) { 
	      rval[i] = pow(V, - 1.0 / xi);
	    } else {
	      if (xi > 0.0) {
		rval[i] = 1.0;
	      } else {
		rval[i] = 0.0;
	      }
	    }
	  }
	} /* non-exponential case */
	
	if (INTEGER(lowerTailFlag)[0]) {
	  rval[i] = 1.0 - rval[i];
	}
	
      }   /* non-NA case     */
      
    }
  
    UNPROTECT(4);
    return(val);
    
  }

}

/* ==========================================================================
 * Quantile function
 * ========================================================================== */

SEXP Call_qGPD2(SEXP p,               /*  double                          */
		SEXP scale,           /*  double                          */
		SEXP shape,           /*  double                          */
		SEXP lowerTailFlag,   /*  integer                         */
		SEXP derivFlag,       /* integer                          */
		SEXP hessianFlag) {   /*  integer                         */
  
  int n, np, nscale, nshape, i, ip, iscale, ishape,
    deriv = INTEGER(derivFlag)[0], hessian = INTEGER(hessianFlag)[0] ;
  
  double eps = 1e-6, q, lq, q1, xi, sigma, V, W;
  
  SEXP val;
  
  PROTECT(p = coerceVector(p, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));
  PROTECT(shape = coerceVector(shape, REALSXP));

  double *rp = REAL(p), *rscale = REAL(scale), 
    *rshape = REAL(shape);

  np = LENGTH(p);	    	
  nscale = LENGTH(scale);		
  nshape = LENGTH(shape);
  
  if ((np == 0) || (nscale == 0) || (nshape == 0)) 			
    return(allocVector(REALSXP, 0));				
  
  n = np;		  			
  if (n < nscale) n = nscale;
  if (n < nshape) n = nshape;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);

  if (deriv) {

    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n * 2));
    double *rgrad = REAL(grad);

    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

 // if (hessian) { 

      SEXP hess;

      PROTECT(hess = allocVector(REALSXP, n * 2 * 2));
      double *rhess = REAL(hess);
      // double dV, d2V;

      // }

    for (i = ip = iscale = ishape = 0;  i < n; 
	 ip = (++ip == np) ? 0 : ip, 	       
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      // as a general rule, only 1 element of the hessian can be non
      // zero
      
      if (hessian) {
	
	// row 'sigma'
	rhess[i] = 0.0;
	rhess[i + n] = 0.0;
	rhess[i + 2 * n] = 0.0;
	// row 'xi'
	rhess[i + 3 * n] = 0.0;
	
      }

      if (ISNA(rp[ip]) || (rscale[iscale] <= 0.0)) {
	// Rprintf("NA case\n");

	rval[i] = NA_REAL;
	
	rgrad[i] = NA_REAL;
	rgrad[i + n] = NA_REAL;
	
	if (hessian) {
	  
	  // row 'sigma'
	  rhess[i] = 0.0;
	  rhess[i + n] = NA_REAL;
	  rhess[i + 2 * n] = NA_REAL;
	  // row 'xi'
	  rhess[i + 3 * n] = NA_REAL;
	}

      } else {

	if (!INTEGER(lowerTailFlag)[0]) {
	  rp[ip] = 1.0 - rp[ip];
	}

	xi = rshape[ishape];
	sigma = rscale[iscale];
	q = 1.0 - rp[ip];
	lq = log(q);

	if (fabs(xi) < eps) {

	  V = - lq;
	  W = lq * lq / 2.0;
	  rval[i] = sigma * V;
	  
	  rgrad[i] =  V;
	  rgrad[i + n] = sigma * W;

	  if (hessian) {
	    // row 'sigma'
	    rhess[i] = 0.0;
	    rhess[i + n] = W;
	    rhess[i + 2 * n] = W;
	    // row 'xi'
	    rhess[i + 3 * n] = - sigma * lq * lq * lq / 3.0;

	  }

	} else {

	  q1 =  pow(q, - xi);
	  V = (q1 - 1.0) / xi;
	  rval[i] = sigma * V;  
	  W = - (V + q1 * lq) / xi;
	  
	  rgrad[i] = V ;
	  rgrad[i + n] = sigma * W;

	  if (hessian) {
	    // row 'sigma'
	    rhess[i] = 0.0;
	    rhess[i + n] = W;
	    rhess[i + 2 * n] = W;
	    // row 'xi'
	    rhess[i + 3 * n] = sigma * (2 * V + q1 * lq * (2 + xi * lq)) / xi / xi;
	  }


	} /* non-exponential case */

      }   /* non-NA case     */
      
    }     /* loop            */

    SET_ATTR(val, attrNm, grad);
    
    if (hessian) {
      SET_STRING_ELT(attrNm, 0, mkChar("hessian"));
      SET_ATTR(val, attrNm, hess);
    } 

    UNPROTECT(7);
    return(val);
    
  } else {

    for (i = ip = iscale = ishape = 0;  i < n; 
	 ip = (++ip == np) ? 0 : ip, 	       
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rp[ip]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
		
      } else {
		
	if (!INTEGER(lowerTailFlag)[0]) {
	  rp[ip] = 1.0 - rp[ip];
	}

	xi = rshape[ishape];
	q = 1 - rp[ip];

	if (fabs(xi) < eps) {

	  rval[i] = -rscale[iscale] * log(q);	  
	  
	} else {

	  V = (pow(q, - xi) - 1.0) / xi;
	  rval[i] = rscale[iscale] * V;  
       	  	  
	} /* non-exp case */
	
      }   /* non-NA case     */
      
    }     /* loop            */
  
    UNPROTECT(4);
    return(val);
    
  }
}
