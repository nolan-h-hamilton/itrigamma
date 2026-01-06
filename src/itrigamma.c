#include <math.h>
#include <float.h>

/*


The following from chapter6, Abramowitz and Stegun 1964, "Handbook of Mathematical Functions" are used to implement
the trigamma, tetragamma, and inverse trigamma functions:


log(Gamma(x)) = \int_0^\infty [ (x-1) e^{-t} - (e^{-t} - e^{-x t}) / (1 - e^{-t}) ] dt / t

- Digamma(x) (psi) is the first derivative of the log-gamma function.
- Trigamma(x) (psi') is the second derivative of the log-gamma function.
- Tetragamma(x) (psi'') is the third derivative of the log-gamma function.
- Polygamma(n, x) is the general nth derivative of the log-gamma function.

Recurrences

---digamma (psi)---
    psi(x+1) = d/dz log(gamma(x+1)) = d/dz log(x * gamma(x)) = d/dz log(x) + d/dz log(gamma(x))
    --> psi(x+1) = (1/x) + psi(x)


---trigamma (psi')---

    psi'(x+1) = d/dz psi(x+1) = d/dz [ (1/x) + psi(x) ] = -1/x^2 + psi'(x)
        --> psi'(x+1) = (1/x^2) + psi'(x)


---tetragamma (psi'')---

    psi''(x+1) = d/dz psi'(x+1) = d/dz [ (1/x^2) + psi'(x) ] = -2/x^3 + psi''(x)
        --> psi''(x+1) = -2/(x^3) + psi''(x)


Asymptotic expansions

---trigamma asymptotic expansion (6.4.12)---

    psi'(z) ~ 1/z + 1/(2 z^2) + 1/(6 z^3) - 1/(30 z^5) + 1/(42 z^7) - 1/(30 z^9) + 5/(66 z^11) - 691/(2730 z^13) + 7/(6 z^15) + ...


---tetragamma asymptotic (6.4.13)---

    psi''(z) ~ -1/z^2 - 1/z^3 - 1/(2 z^4) + 1/(6 z^6) - 1/(6 z^8) + 3/(10 z^10) - 5/(6 z^12) + 691/(210 z^14) - 35/(2 z^16) + ...

*/

static const double SHIFT_TO = 8.0;
static const double INVERSE_TRIGAMMA_Y_BIG = 1e7;
static const double INVERSE_TRIGAMMA_Y_SMALL = 1e-6;
static const double INVERSE_TRIGAMMA_REL_TOL = 1e-12;
static const int INVERSE_TRIGAMMA_NEWTON_MAX_ITER = 50;
static const int INVERSE_TRIGAMMA_BISECT_MAX_ITER = 200;
static const double TRIGAMMA_AE_1 = 1.0;
static const double TRIGAMMA_AE_2 = 0.5;
static const double TRIGAMMA_AE_3 = 1.0 / 6.0;
static const double TRIGAMMA_AE_5 = -1.0 / 30.0;
static const double TRIGAMMA_AE_7 = 1.0 / 42.0;
static const double TRIGAMMA_AE_9 = -1.0 / 30.0;
static const double TRIGAMMA_AE_11 = 5.0 / 66.0;
static const double TRIGAMMA_AE_13 = -691.0 / 2730.0;
static const double TRIGAMMA_AE_15 = 7.0 / 6.0;


static inline double minPositiveXNoOverflow(void)
{
    return 1.0 / sqrt(DBL_MAX);
}

static double trigamma(double x)
{
    double result;
    double reciprocalX;
    double reciprocalX2;
    double reciprocalPower;

    if (!(x > 0.0))
        return NAN;
    if (!isfinite(x))
        return NAN;

    if (x <= minPositiveXNoOverflow())
        return INFINITY;

    result = 0.0;

    while (x < SHIFT_TO) {
        double x2_;

        x2_ = x*x;

        /* check overflow */
        if (!(x2_ > 0.0))
            return INFINITY;
        if (!isfinite(x2_))
            return INFINITY;

        result += 1.0 / x2_;
        x += 1.0;
    }

    /* build up the asymptotic expansion */
    reciprocalX = 1.0 / x;
    reciprocalX2 = reciprocalX * reciprocalX;

    reciprocalPower = reciprocalX;
    result += TRIGAMMA_AE_1 * reciprocalPower;

    reciprocalPower *= reciprocalX;
    result += TRIGAMMA_AE_2 * reciprocalPower;

    reciprocalPower *= reciprocalX;
    result += TRIGAMMA_AE_3 * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += TRIGAMMA_AE_5 * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += TRIGAMMA_AE_7 * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += TRIGAMMA_AE_9 * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += TRIGAMMA_AE_11 * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += TRIGAMMA_AE_13 * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += TRIGAMMA_AE_15 * reciprocalPower;

    return result;
}

static double tetragamma(double x)
{
    double result;
    double reciprocalX;
    double reciprocalX2;
    double reciprocalPower;

    if (!(x > 0.0))
        return NAN;
    if (!isfinite(x))
        return NAN;

    if (x <= minPositiveXNoOverflow())
        return -INFINITY;

    result = 0.0;

    /*
     * tetragamma recurrence
     * psi''(x) = psi''(x+1) - 2/x^3
     */
    while (x < SHIFT_TO) {
        double x3_;

        x3_ = x * x * x;

        if (!(x3_ > 0.0))
            return -INFINITY;
        if (!isfinite(x3_))
            return -INFINITY;

        result -= 2.0 / x3_;
        x += 1.0;
    }

    /*
     * truncated tetragamma asymptotic expansion (6.4.13)
     */
    reciprocalX = 1.0 / x;
    reciprocalX2 = reciprocalX * reciprocalX;

    reciprocalPower = reciprocalX2;
    result += -1.0 * reciprocalPower;

    reciprocalPower *= reciprocalX;
    result += -1.0 * reciprocalPower;

    reciprocalPower *= reciprocalX;
    result += -0.5 * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += (1.0 / 6.0) * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += (-1.0 / 6.0) * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += (3.0 / 10.0) * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += (-5.0 / 6.0) * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += (691.0 / 210.0) * reciprocalPower;

    reciprocalPower *= reciprocalX2;
    result += (-35.0 / 2.0) * reciprocalPower;

    return result;
}

static double inverseTrigamma(double y)
{

    double minX;
    double currentX;
    double lowerBound;
    double upperBound;
    double lowerResidual;
    double upperResidual;

    if (!(y > 0.0))
        return NAN;
    if (!isfinite(y))
        return NAN;

    /*
     * early-exits -- use expansions for large and small
     */
    if (y > INVERSE_TRIGAMMA_Y_BIG)
        return 1.0 / sqrt(y);
    if (y < INVERSE_TRIGAMMA_Y_SMALL)
        return 1.0 / y;

    minX = minPositiveXNoOverflow();
    if (y < 1.0)
        currentX = 1.0 / y + 0.5;
    else
        currentX = 1.0 / sqrt(y);

    if (!(currentX > minX))
        currentX = 1.0;
    if (!isfinite(currentX))
        currentX = 1.0;


    lowerBound = fmax(minX, 0.5 * currentX); /* avoid zero or negative */
    upperBound = fmax(lowerBound * 2.0, currentX); /* ensure upper >= curr, lower */
    lowerResidual = trigamma(lowerBound) - y;

    if (!isfinite(lowerResidual))
        lowerResidual = INFINITY;

    upperResidual = trigamma(upperBound) - y;

    if (!isfinite(upperResidual))
        upperResidual = -y;

    for (int k = 0; k < 256 && upperResidual > 0.0; k++) {
        upperBound *= 2.0;

        if (!isfinite(upperBound))
            return NAN;
        if (upperBound > 1e100)
            return NAN;

        upperResidual = trigamma(upperBound) - y;

        if (!isfinite(upperResidual))
            upperResidual = -y;
    }

    for (int k = 0; k < 256 && lowerResidual < 0.0 && lowerBound > minX; k++) {
        lowerBound *= 0.5;

        if (lowerBound < minX)
            lowerBound = minX;

        lowerResidual = trigamma(lowerBound) - y;

        if (!isfinite(lowerResidual))
            lowerResidual = INFINITY;
    }

    if (!(lowerResidual > 0.0 && upperResidual < 0.0))
        return NAN;

    /*
     * Newton-raphson iteration to find the root: trigamma(x) - y == 0
     */
    for (int iter = 0; iter < INVERSE_TRIGAMMA_NEWTON_MAX_ITER; iter++) {
        double functionValue;
        double grad_;
        double nextX;
        double boundedInterval;
        double xTolerance;

        functionValue = trigamma(currentX) - y;

        /* we're optimizing over [minX, +infinity)*/
        if (!isfinite(functionValue)) {
            if (currentX <= minX)
                functionValue = INFINITY;
            else
                functionValue = -y;
        }

        if (functionValue > 0.0) {
            lowerBound = currentX;
            lowerResidual = functionValue;
        } else if (functionValue < 0.0) {
            upperBound = currentX;
            upperResidual = functionValue;
        } else {
            return currentX;
        }

        boundedInterval = upperBound - lowerBound;

        xTolerance = INVERSE_TRIGAMMA_REL_TOL * fmax(1.0, fabs(currentX));

        if (boundedInterval <= xTolerance)
            return 0.5 * (lowerBound + upperBound);

        grad_ = tetragamma(currentX);
        if (!isfinite(grad_))
            grad_ = 0.0;

        if (grad_ == 0.0) {
            nextX = 0.5 * (lowerBound + upperBound);
        } else {
            nextX = currentX - functionValue / grad_;

            if (!(nextX > lowerBound && nextX < upperBound))
                nextX = 0.5 * (lowerBound + upperBound);

            if (!isfinite(nextX))
                nextX = 0.5 * (lowerBound + upperBound);

            if (!(nextX > 0.0))
                nextX = 0.5 * (lowerBound + upperBound);
        }

        currentX = nextX;
    }

    /* bisect */
    for (int iter = 0; iter < INVERSE_TRIGAMMA_BISECT_MAX_ITER; iter++) {
        double midpoint;
        double midpointResidual;
        double boundedInterval;
        double midpointTolerance;

        midpoint = 0.5 * (lowerBound + upperBound);
        midpointResidual = trigamma(midpoint) - y;

        if (!isfinite(midpointResidual))
            return NAN;

        boundedInterval = upperBound - lowerBound;
        midpointTolerance = INVERSE_TRIGAMMA_REL_TOL * fmax(1.0, fabs(midpoint));

        if (boundedInterval <= midpointTolerance)
            return midpoint;

        if (midpointResidual > 0.0)
            lowerBound = midpoint;
        else
            upperBound = midpoint;
    }

    return 0.5 * (lowerBound + upperBound);
}
