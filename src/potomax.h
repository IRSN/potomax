#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP Call_dGPD2(SEXP x,
		SEXP scale, SEXP shape,
		SEXP logFlag, SEXP derivFlag, SEXP hessianFlag);

SEXP Call_pGPD2(SEXP q, 
		SEXP scale, SEXP shape, 
		SEXP lowerTailFlag, SEXP derivFlag, SEXP hessianFlag);
 
SEXP Call_qGPD2(SEXP p,
		SEXP scale, SEXP shape,
		SEXP lowerTailFlag, SEXP derivFlag, SEXP hessianFlag);

SEXP Call_poisGP2PP(SEXP lambda, SEXP loc, SEXP scale, SEXP shape,
		    SEXP w,
		    SEXP derivFlag);

SEXP Call_PP2poisGP(SEXP locStar, SEXP scaleStar, SEXP shapeStar,
		    SEXP threshold, SEXP w,
		    SEXP derivFlag);
