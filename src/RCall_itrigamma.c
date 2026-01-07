#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include "itrigamma.h"


/* take a function pointer to `f` and return `f(value)` */
static inline double callOnValue(double (*f)(double), double v) {
  return f(v);
}


SEXP itrigamma_trigamma(SEXP x) {
  if (!isReal(x)) x = PROTECT(coerceVector(x, REALSXP));
  else PROTECT(x);

  R_xlen_t n = XLENGTH(x);
  SEXP out = PROTECT(allocVector(REALSXP, n));
  const double *ptrTo_x = REAL(x);
  double *ptrTo_out = REAL(out);

  for (R_xlen_t i = 0; i < n; i++) {
    double v = ptrTo_x[i];
    if (R_IsNA(v)) {
      ptrTo_out[i] = NA_REAL;
    } else {
      ptrTo_out[i] = callOnValue(trigamma, v);
    }
  }

  UNPROTECT(2);
  return out;
}

SEXP itrigamma_tetragamma(SEXP x) {
  if (!isReal(x)) x = PROTECT(coerceVector(x, REALSXP));
  else PROTECT(x);

  R_xlen_t n = XLENGTH(x);
  SEXP out = PROTECT(allocVector(REALSXP, n));

  const double *ptrTo_x = REAL(x);
  double *ptrTo_out = REAL(out);

  for (R_xlen_t i = 0; i < n; i++) {
    double v = ptrTo_x[i];
    if (R_IsNA(v)) ptrTo_out[i] = NA_REAL;
    else {
      ptrTo_out[i] = tetragamma(v);
    }
  }

  UNPROTECT(2);
  return out;
}

SEXP itrigamma_itrigamma(SEXP y) {
  if (!isReal(y)) y = PROTECT(coerceVector(y, REALSXP));
  else PROTECT(y);
  R_xlen_t n = XLENGTH(y);
  /* allocate n*sizeof(double) bytes for output*/
  SEXP out = PROTECT(allocVector(REALSXP, n));
  const double *ptrTo_y = REAL(y);
  double *ptrTo_out = REAL(out);

  for (R_xlen_t i = 0; i < n; i++) {
    double v = ptrTo_y[i];
    if (R_IsNA(v)) ptrTo_out[i] = NA_REAL;
    else {
        ptrTo_out[i] = itrigamma(v);
    }
  }

  UNPROTECT(2);
  return out;
}
