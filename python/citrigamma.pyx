r"""Cython wrapper for `itrigamma` C library

This module provides the python bindings to the `itrigamma` C library. See `itrigamma.c` for the implementation details of the numerical routines.

"""

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, nonecheck=False, initializedcheck=False
import numpy as np
cimport numpy as cnp
from libc.stddef cimport size_t
cnp.import_array()


cdef extern from "itrigamma.h":
    void pos_digamma_vec(const double *x, double *out, size_t n)
    void pos_trigamma_vec(const double *x, double *out, size_t n)
    void pos_tetragamma_vec(const double *x, double *out, size_t n)
    void itrigamma_vec(const double *y, double *out, size_t n)


cdef inline cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] ravelContig_F64(object a):
    cdef object arr = np.asarray(a)
    cdef object ret

    if (<object>arr).dtype == np.dtype(np.float64):
        ret = (<object>arr).ravel(order="K")
        if (<object>ret).flags["C_CONTIGUOUS"]:
            return <cnp.ndarray[cnp.float64_t, ndim=1, mode="c"]> ret

    ret = np.ascontiguousarray((<object>arr).ravel(order="K"), dtype=np.float64)
    return <cnp.ndarray[cnp.float64_t, ndim=1, mode="c"]> ret


cdef inline object applyVec_F64(object x, void (*fn)(const double*, double*, size_t)):
    cdef object arr = np.asarray(x)
    cdef tuple shape = (<object>arr).shape
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] x1
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] out1
    cdef size_t n
    cdef object outObj

    if shape == ():
        # case: scalar input, so skip ravel and the contig check
        x1 = <cnp.ndarray[cnp.float64_t, ndim=1, mode="c"]> np.empty(1, dtype=np.float64)
        out1 = <cnp.ndarray[cnp.float64_t, ndim=1, mode="c"]> np.empty(1, dtype=np.float64)
        (<double*>x1.data)[0] = <double>np.float64((<object>arr)[()])
        fn(<const double*>x1.data, <double*>out1.data, <size_t>1)
        return float((<double*>out1.data)[0])

    x1 = ravelContig_F64(arr)
    outObj = np.empty_like(x1)
    out1 = <cnp.ndarray[cnp.float64_t, ndim=1, mode="c"]> outObj
    n = <size_t>x1.size
    fn(<const double*>x1.data, <double*>out1.data, n)
    return (<object>out1).reshape(shape)


cpdef digamma(x):
    r"""Digamma function :math:`\psi(x)` for positive real arguments.

    .. math::

        \psi(x) = \frac{d}{dx}\log\Gamma(x)

    Numerical Routine
    ------------------

    For :math:`x>0`, the backend combines:

    **Recurrence (shift small x upward)**

    .. math::

        \psi(x) = \psi(x+1) - \frac{1}{x}

    applied repeatedly until :math:`x \geq x_0` (a fixed switching threshold).

    **Asymptotic expansion (Euler–Maclaurin)**

    Euler–Maclaurin summation applied to the tail of the defining series for
    :math:`\psi(x)` yields the large-:math:`x` asymptotic expansion (Bernoulli
    corrections). A truncated form is evaluated once :math:`x \geq x_0`, e.g.

    .. math::

        \psi(x) \sim \log x - \frac{1}{2x} - \frac{1}{12x^2} + \frac{1}{120x^4}
        - \cdots

    The switching threshold :math:`x_0` is chosen conservatively so the truncated
    expansion is near machine precision while keeping the recurrence cheap.

    Parameters
    -----------
    x : array_like
        Positive input values.

    Returns
    --------
    out : float or ndarray
        :math:`\psi(x)` with the same shape as `x`.
    """
    return applyVec_F64(x, pos_digamma_vec)


cpdef trigamma(x):
    r"""Trigamma function :math:`\psi^{(1)}(x)` for positive real arguments.

    .. math::

        \psi^{(1)}(x) = \frac{d}{dx}\psi(x)

    Numerical Routine
    -----------------------------

    **Recurrence (shift small x upward)**

    .. math::

        \psi^{(1)}(x) = \psi^{(1)}(x+1) + \frac{1}{x^2}

    applied until :math:`x \geq x_0`.

    **Asymptotic expansion (Euler–Maclaurin)**

    Using Euler–Maclaurin on the defining series

    .. math::

        \psi^{(1)}(x) = \sum_{k=0}^{\infty}\frac{1}{(x+k)^2}

    gives the large-:math:`x` expansion with Bernoulli corrections; a truncated
    form is evaluated for :math:`x \geq x_0`, e.g.

    .. math::

        \psi^{(1)}(x) \sim \frac{1}{x} + \frac{1}{2x^2} + \frac{1}{6x^3}
        - \frac{1}{30x^5} + \cdots

    Parameters
    -----------
    x : array_like
        Positive input values.

    Returns
    --------
    out : float or ndarray
        :math:`\psi^{(1)}(x)` with the same shape as `x`.

    """

    return applyVec_F64(x, pos_trigamma_vec)


cpdef tetragamma(x):
    r"""Tetragamma function :math:`\psi^{(2)}(x)` for positive real arguments.

    .. math::

        \psi^{(2)}(x) = \frac{d}{dx}\psi^{(1)}(x)


    Numerical Routine
    -----------------------------

    **Recurrence (shift small x upward)**

    .. math::

        \psi^{(2)}(x) = \psi^{(2)}(x+1) - \frac{2}{x^3}

    applied until :math:`x \geq x_0`.

    **Asymptotic expansion (Euler–Maclaurin)**

    Euler–Maclaurin summation yields a large-:math:`x` expansion; a truncated
    form is evaluated for :math:`x \geq x_0`, e.g.

    .. math::

        \psi^{(2)}(x) \sim -\frac{1}{x^2} - \frac{1}{x^3} - \frac{1}{2x^4}
        + \frac{1}{6x^6} - \cdots

    Parameters
    ----------
    x : array_like
        Positive input values.

    Returns
    -------
    out : float or ndarray
        :math:`\psi^{(2)}(x)` with the same shape as `x`.
    """
    return applyVec_F64(x, pos_tetragamma_vec)


cpdef itrigamma(y):
    r"""Inverse trigamma :math:`(\psi^{(1)})^{-1}(y)` for :math:`y>0`.

    Solves for :math:`x>0` in

    .. math::

        \psi^{(1)}(x) = y.

    Numerical Routine
    ----------------------------

    In the C backend, first find the root of :math:`f(x) = \psi^{(1)}(x) - y` on
    :math:`(0,\infty)` using a safeguarded Newton method:

    #. Find an interval :math:`[a,b]` s.t. :math:`f(a)>0` and :math:`f(b)<0`
       (note trigamma is strictly decreasing on :math:`(0,\infty)`).

    #. Newton update using tetragamma as the derivative:

       .. math::

           x_{n+1} = x_n - \frac{\psi^{(1)}(x_n) - y}{\psi^{(2)}(x_n)}.

    #. If the Newton step exits the bracket or produces non-finite values --> go to bisection :math:`x_{n+1}=(a+b)/2`.

    #. Stop when the bracket width is small or :math:`|f(x)|` is below a relative tolerance.

    The quality of Newton steps depends on accurate evaluation of
    :math:`\psi^{(1)}` and :math:`\psi^{(2)}`. See :func:`trigamma` and :func:`tetragamma`.

    Parameters
    ----------
    y : array_like
        Positive target values.

    Returns
    -------
    out : float or ndarray
        Solution :math:`x` with the same shape as `y`.

    """

    cdef object arr = np.asarray(y)
    cdef tuple shape = (<object>arr).shape
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] y1
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] out1
    cdef size_t n
    cdef object outObj

    if shape == ():
        y1 = <cnp.ndarray[cnp.float64_t, ndim=1, mode="c"]> np.empty(1, dtype=np.float64)
        out1 = <cnp.ndarray[cnp.float64_t, ndim=1, mode="c"]> np.empty(1, dtype=np.float64)
        (<double*>y1.data)[0] = <double>np.float64((<object>arr)[()])
        itrigamma_vec(<const double*>y1.data, <double*>out1.data, <size_t>1)
        return float((<double*>out1.data)[0])

    y1 = ravelContig_F64(arr)
    outObj = np.empty_like(y1)
    out1 = <cnp.ndarray[cnp.float64_t, ndim=1, mode="c"]> outObj
    n = <size_t>y1.size
    itrigamma_vec(<const double*>y1.data, <double*>out1.data, n)
    return (<object>out1).reshape(shape)
