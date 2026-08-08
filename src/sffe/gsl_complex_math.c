/* complex/math.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jorma Olavi Tähtinen, Brian
 * Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */

/* Basic complex arithmetic functions

 * Original version by Jorma Olavi Tähtinen <jotahtin@cc.hut.fi>
 *
 * Modified for GSL by Brian Gough, 3/2000
 */

/* The following references describe the methods used in these
 * functions,
 *
 *   T. E. Hull and Thomas F. Fairgrieve and Ping Tak Peter Tang,
 *   "Implementing Complex Elementary Functions Using Exception
 *   Handling", ACM Transactions on Mathematical Software, Volume 20
 *   (1994), pp 215-244, Corrigenda, p553
 *
 *   Hull et al, "Implementing the complex arcsin and arccosine
 *   functions using exception handling", ACM Transactions on
 *   Mathematical Software, Volume 23 (1997) pp 299-335
 *
 *   Abramowitz and Stegun, Handbook of Mathematical Functions, "Inverse
 *   Circular Functions in Terms of Real and Imaginary Parts", Formulas
 *   4.4.37, 4.4.38, 4.4.39
 */

#include "config.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "number_math.h"

/**********************************************************************
 * Complex numbers
 **********************************************************************/

gsl_complex gsl_complex_polar(number_t r, number_t theta)
{ /* return z = r nexp(i theta) */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, r * ncos(theta), r * nsin(theta));
    return z;
}

/**********************************************************************
 * Properties of complex numbers
 **********************************************************************/

number_t gsl_complex_arg(gsl_complex z)
{ /* return arg(z),  -pi < arg(z) <= +pi */
    number_t x = GSL_REAL(z);
    number_t y = GSL_IMAG(z);

    if (x == 0.0 && y == 0.0) {
        return 0;
    }

    return natan2(y, x);
}

number_t gsl_complex_abs(gsl_complex z)
{ /* return |z| */
    return nhypot(GSL_REAL(z), GSL_IMAG(z));
}

number_t gsl_complex_abs2(gsl_complex z)
{ /* return |z|^2 */
    number_t x = GSL_REAL(z);
    number_t y = GSL_IMAG(z);

    return (x * x + y * y);
}

number_t gsl_complex_logabs(gsl_complex z)
{ /* return log|z| */
    number_t xabs = nfabs(GSL_REAL(z));
    number_t yabs = nfabs(GSL_IMAG(z));
    number_t max, u;

    if (xabs >= yabs) {
        max = xabs;
        u = yabs / xabs;
    } else {
        max = yabs;
        u = xabs / yabs;
    }

    /* Handle underflow when u is close to 0 */

    return nlog(max) + 0.5 * nlog1p(u * u);
}

/***********************************************************************
 * Complex arithmetic operators
 ***********************************************************************/

gsl_complex gsl_complex_add(gsl_complex a, gsl_complex b)
{ /* z=a+b */
    number_t ar = GSL_REAL(a), ai = GSL_IMAG(a);
    number_t br = GSL_REAL(b), bi = GSL_IMAG(b);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, ar + br, ai + bi);
    return z;
}

gsl_complex gsl_complex_add_real(gsl_complex a, number_t x)
{ /* z=a+x */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, GSL_REAL(a) + x, GSL_IMAG(a));
    return z;
}

gsl_complex gsl_complex_add_imag(gsl_complex a, number_t y)
{ /* z=a+iy */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, GSL_REAL(a), GSL_IMAG(a) + y);
    return z;
}

gsl_complex gsl_complex_sub(gsl_complex a, gsl_complex b)
{ /* z=a-b */
    number_t ar = GSL_REAL(a), ai = GSL_IMAG(a);
    number_t br = GSL_REAL(b), bi = GSL_IMAG(b);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, ar - br, ai - bi);
    return z;
}

gsl_complex gsl_complex_sub_real(gsl_complex a, number_t x)
{ /* z=a-x */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, GSL_REAL(a) - x, GSL_IMAG(a));
    return z;
}

gsl_complex gsl_complex_sub_imag(gsl_complex a, number_t y)
{ /* z=a-iy */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, GSL_REAL(a), GSL_IMAG(a) - y);
    return z;
}

gsl_complex gsl_complex_mul(gsl_complex a, gsl_complex b)
{ /* z=a*b */
    number_t ar = GSL_REAL(a), ai = GSL_IMAG(a);
    number_t br = GSL_REAL(b), bi = GSL_IMAG(b);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, ar * br - ai * bi, ar * bi + ai * br);
    return z;
}

gsl_complex gsl_complex_mul_real(gsl_complex a, number_t x)
{ /* z=a*x */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, x * GSL_REAL(a), x * GSL_IMAG(a));
    return z;
}

gsl_complex gsl_complex_mul_imag(gsl_complex a, number_t y)
{ /* z=a*iy */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, -y * GSL_IMAG(a), y * GSL_REAL(a));
    return z;
}

gsl_complex gsl_complex_div(gsl_complex a, gsl_complex b)
{ /* z=a/b */
    number_t ar = GSL_REAL(a), ai = GSL_IMAG(a);
    number_t br = GSL_REAL(b), bi = GSL_IMAG(b);

    number_t s = 1.0 / gsl_complex_abs(b);

    number_t sbr = s * br;
    number_t sbi = s * bi;

    number_t zr = (ar * sbr + ai * sbi) * s;
    number_t zi = (ai * sbr - ar * sbi) * s;

    gsl_complex z;
    GSL_SET_COMPLEX(&z, zr, zi);
    return z;
}

gsl_complex gsl_complex_div_real(gsl_complex a, number_t x)
{ /* z=a/x */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, GSL_REAL(a) / x, GSL_IMAG(a) / x);
    return z;
}

gsl_complex gsl_complex_div_imag(gsl_complex a, number_t y)
{ /* z=a/(iy) */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, GSL_IMAG(a) / y, -GSL_REAL(a) / y);
    return z;
}

gsl_complex gsl_complex_conjugate(gsl_complex a)
{ /* z=conj(a) */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, GSL_REAL(a), -GSL_IMAG(a));
    return z;
}

gsl_complex gsl_complex_negative(gsl_complex a)
{ /* z=-a */
    gsl_complex z;
    GSL_SET_COMPLEX(&z, -GSL_REAL(a), -GSL_IMAG(a));
    return z;
}

gsl_complex gsl_complex_inverse(gsl_complex a)
{ /* z=1/a */
    number_t s = 1.0 / gsl_complex_abs(a);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, (GSL_REAL(a) * s) * s, -(GSL_IMAG(a) * s) * s);
    return z;
}

/**********************************************************************
 * Elementary complex functions
 **********************************************************************/

gsl_complex gsl_complex_sqrt(gsl_complex a)
{ /* z=nsqrt(a) */
    gsl_complex z;

    if (GSL_REAL(a) == 0.0 && GSL_IMAG(a) == 0.0) {
        GSL_SET_COMPLEX(&z, 0, 0);
    } else {
        number_t x = nfabs(GSL_REAL(a));
        number_t y = nfabs(GSL_IMAG(a));
        number_t w;

        if (x >= y) {
            number_t t = y / x;
            w = nsqrt(x) * nsqrt(0.5 * (1.0 + nsqrt(1.0 + t * t)));
        } else {
            number_t t = x / y;
            w = nsqrt(y) * nsqrt(0.5 * (t + nsqrt(1.0 + t * t)));
        }

        if (GSL_REAL(a) >= 0.0) {
            number_t ai = GSL_IMAG(a);
            GSL_SET_COMPLEX(&z, w, ai / (2.0 * w));
        } else {
            number_t ai = GSL_IMAG(a);
            number_t vi = (ai >= 0) ? w : -w;
            GSL_SET_COMPLEX(&z, ai / (2.0 * vi), vi);
        }
    }

    return z;
}

gsl_complex gsl_complex_sqrt_real(number_t x)
{ /* z=nsqrt(x) */
    gsl_complex z;

    if (x >= 0) {
        GSL_SET_COMPLEX(&z, nsqrt(x), 0.0);
    } else {
        GSL_SET_COMPLEX(&z, 0.0, nsqrt(-x));
    }

    return z;
}

gsl_complex gsl_complex_exp(gsl_complex a)
{ /* z=nexp(a) */
    number_t rho = nexp(GSL_REAL(a));
    number_t theta = GSL_IMAG(a);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, rho * ncos(theta), rho * nsin(theta));
    return z;
}

gsl_complex gsl_complex_pow(gsl_complex a, gsl_complex b)
{ /* z=a^b */
    gsl_complex z;

    if (GSL_REAL(a) == 0 && GSL_IMAG(a) == 0.0) {
        if (GSL_REAL(b) == 0 && GSL_IMAG(b) == 0.0) {
            GSL_SET_COMPLEX(&z, 1.0, 0.0);
        } else {
            GSL_SET_COMPLEX(&z, 0.0, 0.0);
        }
    } else if (GSL_REAL(b) == 1.0 && GSL_IMAG(b) == 0.0) {
        return a;
    } else if (GSL_REAL(b) == -1.0 && GSL_IMAG(b) == 0.0) {
        return gsl_complex_inverse(a);
    } else {
        number_t logr = gsl_complex_logabs(a);
        number_t theta = gsl_complex_arg(a);

        number_t br = GSL_REAL(b), bi = GSL_IMAG(b);

        number_t rho = nexp(logr * br - bi * theta);
        number_t beta = theta * br + bi * logr;

        GSL_SET_COMPLEX(&z, rho * ncos(beta), rho * nsin(beta));
    }

    return z;
}

gsl_complex gsl_complex_pow_real(gsl_complex a, number_t b)
{ /* z=a^b */
    gsl_complex z;

    if (GSL_REAL(a) == 0 && GSL_IMAG(a) == 0) {
        if (b == 0) {
            GSL_SET_COMPLEX(&z, 1, 0);
        } else {
            GSL_SET_COMPLEX(&z, 0, 0);
        }
    } else {
        number_t logr = gsl_complex_logabs(a);
        number_t theta = gsl_complex_arg(a);
        number_t rho = nexp(logr * b);
        number_t beta = theta * b;
        GSL_SET_COMPLEX(&z, rho * ncos(beta), rho * nsin(beta));
    }

    return z;
}

gsl_complex gsl_complex_log(gsl_complex a)
{ /* z=nlog(a) */
    number_t logr = gsl_complex_logabs(a);
    number_t theta = gsl_complex_arg(a);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, logr, theta);
    return z;
}

gsl_complex gsl_complex_log10(gsl_complex a)
{ /* z = nlog10(a) */
    return gsl_complex_mul_real(gsl_complex_log(a), 1 / nlog(10.));
}

gsl_complex gsl_complex_log_b(gsl_complex a, gsl_complex b)
{
    return gsl_complex_div(gsl_complex_log(a), gsl_complex_log(b));
}

/***********************************************************************
 * Complex trigonometric functions
 ***********************************************************************/

gsl_complex gsl_complex_sin(gsl_complex a)
{ /* z = nsin(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);

    gsl_complex z;

    if (I == 0.0) {
        /* avoid returning negative zero (-0.0) for the imaginary part  */

        GSL_SET_COMPLEX(&z, nsin(R), 0.0);
    } else {
        GSL_SET_COMPLEX(&z, nsin(R) * ncosh(I), ncos(R) * nsinh(I));
    }

    return z;
}

gsl_complex gsl_complex_cos(gsl_complex a)
{ /* z = ncos(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);

    gsl_complex z;

    if (I == 0.0) {
        /* avoid returning negative zero (-0.0) for the imaginary part  */

        GSL_SET_COMPLEX(&z, ncos(R), 0.0);
    } else {
        GSL_SET_COMPLEX(&z, ncos(R) * ncosh(I), nsin(R) * nsinh(-I));
    }

    return z;
}

gsl_complex gsl_complex_tan(gsl_complex a)
{ /* z = ntan(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);

    gsl_complex z;

    if (nfabs(I) < 1) {
        number_t D = npow(ncos(R), 2.0) + npow(nsinh(I), 2.0);

        GSL_SET_COMPLEX(&z, 0.5 * nsin(2 * R) / D, 0.5 * nsinh(2 * I) / D);
    } else {
        number_t D = npow(ncos(R), 2.0) + npow(nsinh(I), 2.0);
        number_t F = 1 + npow(ncos(R) / nsinh(I), 2.0);

        GSL_SET_COMPLEX(&z, 0.5 * nsin(2 * R) / D, 1 / (ntanh(I) * F));
    }

    return z;
}

gsl_complex gsl_complex_sec(gsl_complex a)
{ /* z = sec(a) */
    gsl_complex z = gsl_complex_cos(a);
    return gsl_complex_inverse(z);
}

gsl_complex gsl_complex_csc(gsl_complex a)
{ /* z = csc(a) */
    gsl_complex z = gsl_complex_sin(a);
    return gsl_complex_inverse(z);
}

gsl_complex gsl_complex_cot(gsl_complex a)
{ /* z = cot(a) */
    gsl_complex z = gsl_complex_tan(a);
    return gsl_complex_inverse(z);
}

/**********************************************************************
 * Inverse Complex Trigonometric Functions
 **********************************************************************/

gsl_complex gsl_complex_arcsin(gsl_complex a)
{ /* z = arcsin(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);
    gsl_complex z;

    if (I == 0) {
        z = gsl_complex_arcsin_real(R);
    } else {
        number_t x = nfabs(R), y = nfabs(I);
        number_t r = nhypot(x + 1, y), s = nhypot(x - 1, y);
        number_t A = 0.5 * (r + s);
        number_t B = x / A;
        number_t y2 = y * y;

        number_t real, imag;

        const number_t A_crossover = 1.5, B_crossover = 0.6417;

        if (B <= B_crossover) {
            real = nasin(B);
        } else {
            if (x <= 1) {
                number_t D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
                real = natan(x / nsqrt(D));
            } else {
                number_t Apx = A + x;
                number_t D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
                real = natan(x / (y * nsqrt(D)));
            }
        }

        if (A <= A_crossover) {
            number_t Am1;

            if (x < 1) {
                Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            } else {
                Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

            imag = nlog1p(Am1 + nsqrt(Am1 * (A + 1)));
        } else {
            imag = nlog(A + nsqrt(A * A - 1));
        }

        GSL_SET_COMPLEX(&z, (R >= 0) ? real : -real, (I >= 0) ? imag : -imag);
    }

    return z;
}

gsl_complex gsl_complex_arcsin_real(number_t a)
{ /* z = arcsin(a) */
    gsl_complex z;

    if (nfabs(a) <= 1.0) {
        GSL_SET_COMPLEX(&z, nasin(a), 0.0);
    } else {
        if (a < 0.0) {
            GSL_SET_COMPLEX(&z, -N_PI_2, nacosh(-a));
        } else {
            GSL_SET_COMPLEX(&z, N_PI_2, -nacosh(a));
        }
    }

    return z;
}

gsl_complex gsl_complex_arccos(gsl_complex a)
{ /* z = arccos(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);
    gsl_complex z;

    if (I == 0) {
        z = gsl_complex_arccos_real(R);
    } else {
        number_t x = nfabs(R), y = nfabs(I);
        number_t r = nhypot(x + 1, y), s = nhypot(x - 1, y);
        number_t A = 0.5 * (r + s);
        number_t B = x / A;
        number_t y2 = y * y;

        number_t real, imag;

        const number_t A_crossover = 1.5, B_crossover = 0.6417;

        if (B <= B_crossover) {
            real = nacos(B);
        } else {
            if (x <= 1) {
                number_t D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
                real = natan(nsqrt(D) / x);
            } else {
                number_t Apx = A + x;
                number_t D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
                real = natan((y * nsqrt(D)) / x);
            }
        }

        if (A <= A_crossover) {
            number_t Am1;

            if (x < 1) {
                Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            } else {
                Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

            imag = nlog1p(Am1 + nsqrt(Am1 * (A + 1)));
        } else {
            imag = nlog(A + nsqrt(A * A - 1));
        }

        GSL_SET_COMPLEX(&z, (R >= 0) ? real : N_PI - real,
                        (I >= 0) ? -imag : imag);
    }

    return z;
}

gsl_complex gsl_complex_arccos_real(number_t a)
{ /* z = arccos(a) */
    gsl_complex z;

    if (nfabs(a) <= 1.0) {
        GSL_SET_COMPLEX(&z, nacos(a), 0);
    } else {
        if (a < 0.0) {
            GSL_SET_COMPLEX(&z, N_PI, -nacosh(-a));
        } else {
            GSL_SET_COMPLEX(&z, 0, nacosh(a));
        }
    }

    return z;
}

gsl_complex gsl_complex_arctan(gsl_complex a)
{ /* z = arctan(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);
    gsl_complex z;

    if (I == 0) {
        GSL_SET_COMPLEX(&z, natan(R), 0);
    } else {
        /* FIXME: This is a naive implementation which does not fully
           take into account cancellation errors, overflow, underflow
           etc.  It would benefit from the Hull et al treatment. */

        number_t r = nhypot(R, I);

        number_t imag;

        number_t u = 2 * I / (1 + r * r);

        /* FIXME: the following cross-over should be optimized but 0.1
           seems to work ok */

        if (nfabs(u) < 0.1) {
            imag = 0.25 * (nlog1p(u) - nlog1p(-u));
        } else {
            number_t A = nhypot(R, I + 1);
            number_t B = nhypot(R, I - 1);
            imag = 0.5 * nlog(A / B);
        }

        if (R == 0) {
            if (I > 1) {
                GSL_SET_COMPLEX(&z, N_PI_2, imag);
            } else if (I < -1) {
                GSL_SET_COMPLEX(&z, -N_PI_2, imag);
            } else {
                GSL_SET_COMPLEX(&z, 0, imag);
            };
        } else {
            GSL_SET_COMPLEX(&z, 0.5 * natan2(2 * R, ((1 + r) * (1 - r))), imag);
        }
    }

    return z;
}

gsl_complex gsl_complex_arcsec(gsl_complex a)
{ /* z = arcsec(a) */
    gsl_complex z = gsl_complex_inverse(a);
    return gsl_complex_arccos(z);
}

gsl_complex gsl_complex_arcsec_real(number_t a)
{ /* z = arcsec(a) */
    gsl_complex z;

    if (a <= -1.0 || a >= 1.0) {
        GSL_SET_COMPLEX(&z, nacos(1 / a), 0.0);
    } else {
        if (a >= 0.0) {
            GSL_SET_COMPLEX(&z, 0, nacosh(1 / a));
        } else {
            GSL_SET_COMPLEX(&z, N_PI, -nacosh(-1 / a));
        }
    }

    return z;
}

gsl_complex gsl_complex_arccsc(gsl_complex a)
{ /* z = arccsc(a) */
    gsl_complex z = gsl_complex_inverse(a);
    return gsl_complex_arcsin(z);
}

gsl_complex gsl_complex_arccsc_real(number_t a)
{ /* z = arccsc(a) */
    gsl_complex z;

    if (a <= -1.0 || a >= 1.0) {
        GSL_SET_COMPLEX(&z, nasin(1 / a), 0.0);
    } else {
        if (a >= 0.0) {
            GSL_SET_COMPLEX(&z, N_PI_2, -nacosh(1 / a));
        } else {
            GSL_SET_COMPLEX(&z, -N_PI_2, nacosh(-1 / a));
        }
    }

    return z;
}

gsl_complex gsl_complex_arccot(gsl_complex a)
{ /* z = arccot(a) */
    gsl_complex z;

    if (GSL_REAL(a) == 0.0 && GSL_IMAG(a) == 0.0) {
        GSL_SET_COMPLEX(&z, N_PI_2, 0);
    } else {
        z = gsl_complex_inverse(a);
        z = gsl_complex_arctan(z);
    }

    return z;
}

/**********************************************************************
 * Complex Hyperbolic Functions
 **********************************************************************/

gsl_complex gsl_complex_sinh(gsl_complex a)
{ /* z = nsinh(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, nsinh(R) * ncos(I), ncosh(R) * nsin(I));
    return z;
}

gsl_complex gsl_complex_cosh(gsl_complex a)
{ /* z = ncosh(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, ncosh(R) * ncos(I), nsinh(R) * nsin(I));
    return z;
}

gsl_complex gsl_complex_tanh(gsl_complex a)
{ /* z = ntanh(a) */
    number_t R = GSL_REAL(a), I = GSL_IMAG(a);

    gsl_complex z;

    if (nfabs(R) < 1.0) {
        number_t D = npow(ncos(I), 2.0) + npow(nsinh(R), 2.0);

        GSL_SET_COMPLEX(&z, nsinh(R) * ncosh(R) / D, 0.5 * nsin(2 * I) / D);
    } else {
        number_t D = npow(ncos(I), 2.0) + npow(nsinh(R), 2.0);
        number_t F = 1 + npow(ncos(I) / nsinh(R), 2.0);

        GSL_SET_COMPLEX(&z, 1.0 / (ntanh(R) * F), 0.5 * nsin(2 * I) / D);
    }

    return z;
}

gsl_complex gsl_complex_sech(gsl_complex a)
{ /* z = sech(a) */
    gsl_complex z = gsl_complex_cosh(a);
    return gsl_complex_inverse(z);
}

gsl_complex gsl_complex_csch(gsl_complex a)
{ /* z = csch(a) */
    gsl_complex z = gsl_complex_sinh(a);
    return gsl_complex_inverse(z);
}

gsl_complex gsl_complex_coth(gsl_complex a)
{ /* z = coth(a) */
    gsl_complex z = gsl_complex_tanh(a);
    return gsl_complex_inverse(z);
}

/**********************************************************************
 * Inverse Complex Hyperbolic Functions
 **********************************************************************/

gsl_complex gsl_complex_arcsinh(gsl_complex a)
{ /* z = arcsinh(a) */
    gsl_complex z = gsl_complex_mul_imag(a, 1.0);
    z = gsl_complex_arcsin(z);
    z = gsl_complex_mul_imag(z, -1.0);
    return z;
}

gsl_complex gsl_complex_arccosh(gsl_complex a)
{ /* z = arccosh(a) */
    gsl_complex z = gsl_complex_arccos(a);
    z = gsl_complex_mul_imag(z, GSL_IMAG(z) > 0 ? -1.0 : 1.0);
    return z;
}

gsl_complex gsl_complex_arccosh_real(number_t a)
{ /* z = arccosh(a) */
    gsl_complex z;

    if (a >= 1) {
        GSL_SET_COMPLEX(&z, nacosh(a), 0);
    } else {
        if (a >= -1.0) {
            GSL_SET_COMPLEX(&z, 0, nacos(a));
        } else {
            GSL_SET_COMPLEX(&z, nacosh(-a), N_PI);
        }
    }

    return z;
}

gsl_complex gsl_complex_arctanh(gsl_complex a)
{ /* z = arctanh(a) */
    if (GSL_IMAG(a) == 0.0) {
        return gsl_complex_arctanh_real(GSL_REAL(a));
    } else {
        gsl_complex z = gsl_complex_mul_imag(a, 1.0);
        z = gsl_complex_arctan(z);
        z = gsl_complex_mul_imag(z, -1.0);
        return z;
    }
}

gsl_complex gsl_complex_arctanh_real(number_t a)
{ /* z = arctanh(a) */
    gsl_complex z;

    if (a > -1.0 && a < 1.0) {
        GSL_SET_COMPLEX(&z, natanh(a), 0);
    } else {
        GSL_SET_COMPLEX(&z, natanh(1 / a), (a < 0) ? N_PI_2 : -N_PI_2);
    }

    return z;
}

gsl_complex gsl_complex_arcsech(gsl_complex a)
{ /* z = arcsech(a); */
    gsl_complex t = gsl_complex_inverse(a);
    return gsl_complex_arccosh(t);
}

gsl_complex gsl_complex_arccsch(gsl_complex a)
{ /* z = arccsch(a) */
    gsl_complex t = gsl_complex_inverse(a);
    return gsl_complex_arcsinh(t);
}

gsl_complex gsl_complex_arccoth(gsl_complex a)
{ /* z = arccoth(a) */
    gsl_complex t = gsl_complex_inverse(a);
    return gsl_complex_arctanh(t);
}
