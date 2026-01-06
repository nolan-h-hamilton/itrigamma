#ifndef ITRIGAMMA_H
#define ITRIGAMMA_H

/* in case we call from Rcpp */
#ifdef __cplusplus
extern "C" {
#endif

double trigamma(double x);
double tetragamma(double x);
double itrigamma(double y);

#ifdef __cplusplus
}
#endif

#endif
