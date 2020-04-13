#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP Call_dGPD2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Call_pGPD2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Call_qGPD2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Call_poisGP2PP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Call_PP2poisGP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"Call_dGPD2", (DL_FUNC) &Call_dGPD2, 6},
  {"Call_pGPD2", (DL_FUNC) &Call_pGPD2, 6},
  {"Call_qGPD2", (DL_FUNC) &Call_qGPD2, 6},
  {"Call_poisGP2PP", (DL_FUNC) &Call_poisGP2PP, 6},
  {"Call_PP2poisGP", (DL_FUNC) &Call_PP2poisGP, 6},
  {NULL, NULL, 0}
};

void R_init_potomax(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
