#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP graphicalVAR_Beta_C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP graphicalVAR_beta_ridge_C(SEXP, SEXP, SEXP);
extern SEXP graphicalVAR_LogLik_and_BIC(SEXP, SEXP, SEXP);
extern SEXP graphicalVAR_VAR_logLik_C(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"graphicalVAR_Beta_C",         (DL_FUNC) &graphicalVAR_Beta_C,         8},
  {"graphicalVAR_beta_ridge_C",   (DL_FUNC) &graphicalVAR_beta_ridge_C,   3},
  {"graphicalVAR_LogLik_and_BIC", (DL_FUNC) &graphicalVAR_LogLik_and_BIC, 3},
  {"graphicalVAR_VAR_logLik_C",   (DL_FUNC) &graphicalVAR_VAR_logLik_C,   4},
  {NULL, NULL, 0}
};

void R_init_graphicalVAR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}