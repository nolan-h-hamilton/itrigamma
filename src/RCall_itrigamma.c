#include "itrigamma.h"
#include <R.h>
#include <R_ext/Arith.h>
#include <Rinternals.h>

SEXP itrigamma_pos_trigamma(SEXP x)
{
    int protectCount = 0;
    // block GC only if x is not real/needs coercion
    if (!isReal(x)) {
      x = PROTECT(coerceVector(x, REALSXP)); protectCount++;
    }

    R_xlen_t n = XLENGTH(x);
    SEXP out = PROTECT(allocVector(REALSXP, n)); protectCount++;

    const double *ptrTo_x = REAL(x);
    double *ptrTo_out = REAL(out);

    for (R_xlen_t i = 0; i < n; i++) {
        double v = ptrTo_x[i];
        if (R_IsNA(v)) {
            ptrTo_out[i] = NA_REAL;
        } else if (R_IsNaN(v)) {
            ptrTo_out[i] = R_NaN;
        } else {
            ptrTo_out[i] = pos_trigamma(v);
        }
    }

    // allow gc for x, out
    UNPROTECT(protectCount);
    return out;
}

SEXP itrigamma_pos_tetragamma(SEXP x)
{
    int protectCount = 0;
    if (!isReal(x)) {
      x = PROTECT(coerceVector(x, REALSXP)); protectCount++;
    }

    R_xlen_t n = XLENGTH(x);
    SEXP out = PROTECT(allocVector(REALSXP, n)); protectCount++;

    const double *ptrTo_x = REAL(x);
    double *ptrTo_out = REAL(out);
    for (R_xlen_t i = 0; i < n; i++) {
        double v = ptrTo_x[i];

        if (R_IsNA(v)) {
            ptrTo_out[i] = NA_REAL;
        } else if (R_IsNaN(v)) {
            ptrTo_out[i] = R_NaN;
        } else {
            ptrTo_out[i] = pos_tetragamma(v);
        }
    }
    // allow gc for x, out
    UNPROTECT(protectCount);
    return out;
}

SEXP itrigamma_itrigamma(SEXP y)
{
    int protectCount = 0;
    if (!isReal(y)) {
        y = PROTECT(coerceVector(y, REALSXP));
        protectCount++;
    }

    R_xlen_t n = XLENGTH(y);
    SEXP out = PROTECT(allocVector(REALSXP, n)); protectCount++;

    const double *ptrTo_y = REAL(y);
    double *ptrTo_out = REAL(out);
    for (R_xlen_t i = 0; i < n; i++) {
        double v = ptrTo_y[i];
        if (R_IsNA(v)) {
            ptrTo_out[i] = NA_REAL;
        } else if (R_IsNaN(v)) {
            ptrTo_out[i] = R_NaN;
        } else {
            ptrTo_out[i] = itrigamma(v);
        }
    }

    UNPROTECT(protectCount);
    return out;
}
