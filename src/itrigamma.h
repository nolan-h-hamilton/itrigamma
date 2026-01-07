#ifndef ITRIGAMMA_H
#define ITRIGAMMA_H

#include <stddef.h>

double trigamma(double x);
double tetragamma(double x);
double itrigamma(double y);

/* array-level versions for each */
void trigamma_vec(const double *x, double *out, size_t n);
void tetragamma_vec(const double *x, double *out, size_t n);
void itrigamma_vec(const double *y, double *out, size_t n);

#endif

