#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#define NODEBUG

/* ===========================================================================
   AUTHOR Yves Deville <deville.yves@alpestat.com>
   
   GOAL 

   Compute the transformation of Poisson-GP parameters into PP
   parameters and its Jacobian.

   INPUT 

   'lambda', 'loc', 'scale' and 'shape' are vectors coped with by
   using the recycling rule. They are used in this order for the
   differentiation.  'w' is assumed to be of length one and is not
   used in differentiation.  'lambda' is the Poisson rate, 'loc,
   'scale' and 'shape' are the GP parameters and 'w' is the block
   duration so that the product lambda * w is dimensionless.

   RESULT 

   A numeric vector to be coerced into an array of dim (n, 3). The
   columns correspond to the 3 GEV parameters 'locStar', 'scaleStar'
   and 'shapeStar'.

   When a derivative is required it is returned as a "gradient"
   attribute of the result. This is a vector to be coerced into an
   array of dimension (n, 3, 4) with indices : index in the input
   vectors, GEV parameter, then Poisson-GP parameter. The element [i1,
   i2, i3] corresponds to the derivative of the i-th GEV vector, j-th
   element and differentiation w.r.t. It corresponds to the element i1
   + n * i2 + (3 * n) * i3 in the vector indexation, according to the
   usual rule "the first index varies the fastest".

   NOTE 

   Since 'lambda' and 'w' are used only through their product, a
   variable 'w' could be used by changing the input: replace 'lambda'
   by 'lambda * w' and 'w' by 1.0.

   DETAILS 

   See 'Renext Computing Details'

   LICENCE 

   Contact the author for details. The 'potomax' R package is still in a
   early stage of development.
   =========================================================================== */  


SEXP Call_poisGP2PP(SEXP lambda,        /*  double                          */
		    SEXP loc,           /*  double                          */
		    SEXP scale,         /*  double                          */
		    SEXP shape,         /*  double                          */
		    SEXP w,             /*  double                          */  
		    SEXP derivFlag) {   /*  integer                         */ 

  int n, n3, nlambda, nloc, nscale, nshape, 
    i, i12, ilambda, iloc, iscale, ishape, 
    deriv = INTEGER(derivFlag)[0];
  
  double eps = 1e-6, w0, xi, sigma, A, L, E, BC;
  
  SEXP val;

  PROTECT(lambda = coerceVector(lambda, REALSXP));
  PROTECT(loc = coerceVector(loc, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));
  PROTECT(shape = coerceVector(shape, REALSXP));
  PROTECT(w = coerceVector(w, REALSXP));

  double *rlambda = REAL(lambda), *rloc = REAL(loc), *rscale = REAL(scale), 
    *rshape = REAL(shape);

  nlambda = LENGTH(lambda);
  nloc = LENGTH(loc);
  nscale = LENGTH(scale);		
  nshape = LENGTH(shape);

  w0 = REAL(w)[0];
  
  if ((nlambda == 0) || (nloc == 0) || (nscale == 0) || (nshape == 0)) 			
    return(allocVector(REALSXP, 0));				
  
  n = nlambda;
  if (n < nloc) n = nloc;					       			
  if (n < nscale) n = nscale;
  if (n < nshape) n = nshape;
  n3 = 3 * n;
  
  PROTECT(val = allocVector(REALSXP, n * 3));
  double *rval = REAL(val);
  
  if (deriv) {	  
	
    SEXP grad, attrNm;
    PROTECT(grad = allocVector(REALSXP, n * 3 * 4));
    double *rgrad = REAL(grad);	  
    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
 
    for (i = ilambda = iloc = iscale = ishape = 0;  i < n; 
	 ilambda = (++ilambda == nlambda) ? 0 : ilambda,
	   iloc = (++iloc == nloc) ? 0 : iloc,
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if ((rscale[iscale] <= 0.0)) {
	// do something here
      }

      
      sigma = rscale[iscale];
      A = rlambda[ilambda] * w0;
      xi = rshape[ishape];
      L = log(A);
      
      if (fabs(xi) < eps) {
	
	rval[i] = rloc[iloc] + L * sigma;  // locStar
	rval[i + n] = sigma;               // scaleStar
	rval[i + 2 * n] = xi;              // shapeStar
	
	
	// row #1      i1 = i, i2 = 0, i3 = 0 to 3
	i12 = i;
	rgrad[i12] = sigma / rlambda[ilambda];               // @locStar / @lambda
	rgrad[i12 + n3] = 1.0;                               // @locStar / @loc
	rgrad[i12 + n3 * 2] = L;                             // @locStar / @scale
	rgrad[i12 + n3 * 3] = sigma * L * L / 2.0;           // @locStar / @shape
	
	// row #2     i1 = i, i2 = 1, i3 = 0 to 3
	i12 = i + n;
	rgrad[i12] = 0.0;                                    // @scaleStar / @lambda
	rgrad[i12 + n3] = 0.0;                               // @scaleStar / @loc
	rgrad[i12 + n3 * 2] = 1.0;                           // @scaleStar / @scale
	rgrad[i12 + n3 * 3] = sigma * L;                     // @scaleStar / @shape
	
	// row #3     i2 = 2, i2 = 2, i3 = 0 to 3
	i12 = i + 2 * n;
	rgrad[i12] = 0.0;                                    // @shapeStar / @lambda
	rgrad[i12 + n3] = 0.0;                               // @hapeStar / @loc
	rgrad[i12 + n3 * 2] = 0.0;                           // @shapeStar / @scale
	rgrad[i12 + n3 * 3] = 1.0;                           // @shapeStar / @shape
	
	
	
      } else {
      
	E = exp(xi * L);          // (lambda * w) ^ xi
	BC = (E - 1.0) / xi;      // Box-Cox

	rval[i] = rloc[iloc] + BC * sigma;     // locStar
	rval[i + n] = E * sigma;              // scaleStar
	rval[i + 2 * n] = xi;                 // shapeStar
	  
	// row #1      i1 = i, i2 = 0, i3 = 0 to 3
	i12 = i;
	rgrad[i12] = E * sigma / rlambda[ilambda];           // @locStar / @lambda
	rgrad[i12 + n3] = 1.0;                               // @locStar / @loc
	rgrad[i12 + n3 * 2] = BC;                            // @locStar / @scale
	rgrad[i12 + n3 * 3] = sigma * (E * L - BC) / xi;     // @locStar / @shape
	
	// row #2     i1 = i, i2 = 1, i3 = 0 to 3
	i12 = i + n;
	rgrad[i12] = xi * rgrad[i];                          // @scaleStar / @lambda
	rgrad[i12 + n3] = 0.0;                               // @scaleStar / @loc
	rgrad[i12 + n3 * 2] = E;                             // @scaleStar / @scale
	rgrad[i12 + n3 * 3] = sigma * E * L;                 // @scaleStar / @shape
	
	// row #3     i2 = 2, i2 = 2, i3 = 0 to 3
	i12 = i + 2 * n;
	rgrad[i12] = 0.0;                                    // @shapeStar / @lambda
	rgrad[i12 + n3] = 0.0;                               // @hapeStar / @loc
	rgrad[i12 + n3 * 2] = 0.0;                           // @shapeStar / @scale
	rgrad[i12 + n3 * 3] = 1.0;                           // @shapeStar / @shape
	
      }
	
    } // i++ loop

    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
    SET_ATTR(val, attrNm, grad);
    UNPROTECT(8);
    return(val);
    
  } else {


    for (i = ilambda = iloc = iscale = ishape = 0;  i < n; 
	 ilambda = (++ilambda == nlambda) ? 0 : ilambda,
	   iloc = (++iloc == nloc) ? 0 : iloc,
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if ((rscale[iscale] <= 0.0)) {
	// do something here
      }
      
      sigma = rscale[iscale];
      A = rlambda[ilambda] * w0;
      xi = rshape[ishape];
      L = log(A);
      
      if (fabs(xi) < eps) {
	
	rval[i] = rloc[iloc] + L * sigma;  // locStar
	rval[i + n] = sigma;               // scaleStar
	rval[i + 2 * n] = xi;              // shapeStar	
	
      } else {
      
	E = exp(xi * L);          // (lambda * w) ^ xi
	BC = (E - 1.0) / xi;      // Box-Cox
	
	rval[i] = rloc[iloc] + BC * sigma;    // locStar
	rval[i + n] = E * sigma;              // scaleStar
	rval[i + 2 * n] = xi;                 // shapeStar
	 
      }
	
    } // i++ loop
    
    UNPROTECT(6);
    return(val);
    
  }

  
}

    
/* ===========================================================================
   AUTHOR Yves Deville <deville.yves@alpestat.com>
   
   GOAL 

   Compute the transformation of PP parameters into Poisson-GP
   parameters and its Jacobian.

   INPUT 

   'loc', 'scale', 'shape' and 'threshold' are vectors coped with by
   using the recycling rule. They are used in this order for the
   differentiation.  'w' is assumed to be of length one and is not
   used in differentiation.  'loc, 'scale' and 'shape' are the PP
   parameters, 'threshold' is the Poisson-GP threshold (or location
   parameter) and 'w' is the block duration so that the product lambda
   * w is dimensionless.

   RESULT 

   A numeric vector to be coerced into an array of dim (n, 3). The
   columns correspond to the 3 GEV parameters 'lambda', 'scale' and
   'shape'.

   When a derivative is required it is returned as a "gradient"
   attribute of the result. This is a vector to be coerced into an
   array say 'a' of dimension (n, 3, 3) with indices : index in the
   input vectors, Poisson-GP parameter, then PP parameter. So each
   slice a[i1, , ] is a Jacobian matrix.  The element [i1, i2, i3]
   corresponds to the derivative of the i1-th Poisson-GP vector, i2-th
   element and differentiation w.r.t the i3-th PP parameter. It
   corresponds to the element i1 + n * i2 + (3 * n) * i3 in the vector
   indexation, according to the usual rule "the first index varies the
   fastest".

   NOTE 

   Since 'lambda' and 'w' are used only through their product, a
   variable 'w' could be used by changing the input: replace 'lambda'
   by 'lambda * w' and 'w' by 1.0.

  =========================================================================== */

SEXP Call_PP2poisGP(SEXP locStar,           /*  double                          */
		    SEXP scaleStar,         /*  double                          */
		    SEXP shapeStar,         /*  double                          */
		    SEXP threshold,         /*  double                          */ 
		    SEXP w,                 /*  double                          */  
		    SEXP derivFlag) {       /*  integer                         */ 
  
  int n, n3, nlocStar, nscaleStar, nshapeStar, nthreshold,
    i, i12, ilocStar, iscaleStar, ishapeStar, ithreshold,
    deriv = INTEGER(derivFlag)[0];
  
  double eps = 1e-6, w0, xiStar, sigmaStar, C, L, z, lambda;
  
  SEXP val;
  
  PROTECT(locStar = coerceVector(locStar, REALSXP));
  PROTECT(scaleStar = coerceVector(scaleStar, REALSXP));
  PROTECT(shapeStar = coerceVector(shapeStar, REALSXP));
  PROTECT(threshold = coerceVector(threshold, REALSXP));
  PROTECT(w = coerceVector(w, REALSXP));
  
  double *rlocStar = REAL(locStar), *rscaleStar = REAL(scaleStar), 
    *rshapeStar = REAL(shapeStar), *rthreshold = REAL(threshold);
  
  nlocStar = LENGTH(locStar);
  nscaleStar = LENGTH(scaleStar);		
  nshapeStar = LENGTH(shapeStar);
  nthreshold = LENGTH(threshold);
  
  w0 = REAL(w)[0];
  
  if ((nlocStar == 0) || (nscaleStar == 0) || (nshapeStar == 0) || (nthreshold == 0)) 			
    return(allocVector(REALSXP, 0));				

  // 'n' is the max of 'nloc', 'nscale', 'nshape' and 'nthreshold'
  n = nlocStar;					       			
  if (n < nscaleStar) n = nscaleStar;
  if (n < nshapeStar) n = nshapeStar;
  if (n < nthreshold) n = nthreshold;
  n3 = 3 * n;
  
  PROTECT(val = allocVector(REALSXP, n * 3));
  double *rval = REAL(val);

  if (deriv) {	  
	
    SEXP grad, attrNm;
    PROTECT(grad = allocVector(REALSXP, n * 3 * 3));
    double *rgrad = REAL(grad);	  
    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
 
    for (i = ilocStar = iscaleStar = ishapeStar = ithreshold = 0; i < n; 
	 ilocStar = (++ilocStar == nlocStar) ? 0 : ilocStar,
	   iscaleStar = (++iscaleStar == nscaleStar) ? 0 : iscaleStar,
	   ishapeStar = (++ishapeStar == nshapeStar) ? 0 : ishapeStar,
	   ithreshold = (++ithreshold == nthreshold) ? 0 : ithreshold,
	   ++i) {
      
      if ((rscaleStar[iscaleStar] <= 0.0)) {
	// do something here
      }

      z = (rthreshold[ithreshold] - rlocStar[ilocStar]) / rscaleStar[iscaleStar];
      xiStar = rshapeStar[ishapeStar];
      // Rprintf("ishapeStar = %d, xiStar = %8.2f, z = %8.2f\n", ishapeStar, xiStar, z);
 
      if (fabs(xiStar) < eps) {

	// in this case we heve C = 1.0
	lambda = exp(-z) / w0;
	sigmaStar = rscaleStar[iscaleStar];
	rval[i] =  lambda;                   // lambda
	rval[i + n] = sigmaStar;             // scale
	rval[i + 2 * n] = xiStar;            // shape
	
	// row #1      i1 = i, i2 = 0, i3 = 0 to 2
	i12 = i;
	rgrad[i12] = lambda / sigmaStar;                     // @lambda / @locStar
	rgrad[i12 + n3] = z * rgrad[i12];                    // @lambda / @scaleStar
	rgrad[i12 + n3 * 2] = lambda * z * z / 2.0 ;         // @lambda / @shapeStar
	
	// row #2     i1 = i, i2 = 1, i3 = 0 to 2
	i12 = i + n;
	rgrad[i12] = 0.0;                                    // @scale / @locStar
	rgrad[i12 + n3] = 1.0;                               // @scale / @scaleStar
	rgrad[i12 + n3 * 2] = sigmaStar * z;                 // @scale / @shapeStar
	
	// row #3     i2 = 2, i2 = 2, i3 = 0 to 2
	i12 = i + 2 * n;
	rgrad[i12] = 0.0;                                    // @shape / @locStar
	rgrad[i12 + n3] = 0.0;                               // @shape / @scaleStar
	rgrad[i12 + n3 * 2] = 1.0;                           // @shape / @shapeStar

	
      } else {

	C = 1.0 + xiStar * z;
	L = log(C);
	
	if (C <= 0.0) {
	  error("each 'threshold' must be in the support of the GEV distribution");
	}

	lambda = pow(C, - 1.0 / xiStar) / w0;
	sigmaStar = rscaleStar[iscaleStar];
	rval[i] =  lambda;                   // lambda
	rval[i + n] = C * sigmaStar;         // scale
	rval[i + 2 * n] = xiStar;            // shape
	
	// row #1      i1 = i, i2 = 0, i3 = 0 to 2
	i12 = i;
	rgrad[i12] = lambda / C / sigmaStar;                                    // @lambda / @locStar
	rgrad[i12 + n3] = z * rgrad[i12];                                       // @lambda / @scaleStar
	rgrad[i12 + n3 * 2] = lambda * (L - xiStar * z / C) / xiStar / xiStar ; // @lambda / @shapeStar
	
	// row #2     i1 = i, i2 = 1, i3 = 0 to 2
	i12 = i + n;
	rgrad[i12] = -xiStar;                                // @scale / @locStar
	rgrad[i12 + n3] = 1.0;                               // @scale / @scaleStar
	rgrad[i12 + n3 * 2] = sigmaStar * z;                 // @scale / @shapeStar
	
	// row #3     i2 = 2, i2 = 2, i3 = 0 to 2
	i12 = i + 2 * n;
	rgrad[i12] = 0.0;                                    // @shape / @locStar
	rgrad[i12 + n3] = 0.0;                               // @shape / @scaleStar
	rgrad[i12 + n3 * 2] = 1.0;                           // @shape / @shapeStar
	
      }
     
    } // i++ loop
    
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
    SET_ATTR(val, attrNm, grad);
    UNPROTECT(8);
    return(val);

  } else {
    
    for (i = ilocStar = iscaleStar = ishapeStar = ithreshold = 0; i < n; 
	 ilocStar = (++ilocStar == nlocStar) ? 0 : ilocStar,
	   iscaleStar = (++iscaleStar == nscaleStar) ? 0 : iscaleStar,
	   ishapeStar = (++ishapeStar == nshapeStar) ? 0 : ishapeStar,
	   ithreshold = (++ithreshold == nthreshold) ? 0 : ithreshold,
	   ++i) {

      
      if ((rscaleStar[iscaleStar] <= 0.0)) {
	// do something here
      }
      
      z = (rthreshold[ithreshold] - rlocStar[ilocStar]) / rscaleStar[iscaleStar];
      xiStar = rshapeStar[ishapeStar];

      if (fabs(xiStar) < eps) {
	
	// in this case we heve C = 1.0
	lambda = exp(-z) / w0;
	sigmaStar = rscaleStar[iscaleStar];
	rval[i] =  lambda;                   // lambda
	rval[i + n] = sigmaStar;             // scale
	rval[i + 2 * n] = xiStar;            // shape
	
      } else {

	C = 1.0 + xiStar * z;
	L = log(C);
	
	if (C <= 0.0) {
	  error("each 'threshold' must be in the support of the GEV distribution");
	}

	lambda = pow(C, - 1.0 / xiStar) / w0;
	sigmaStar = rscaleStar[iscaleStar];
	rval[i] =  lambda;                   // lambda
	rval[i + n] = C * sigmaStar;         // scale
	rval[i + 2 * n] = xiStar;            // shape
	
      }
     
    } // i++ loop
        
    UNPROTECT(6);
    return(val);
    
  } // else deriv != 0
  
  
}
