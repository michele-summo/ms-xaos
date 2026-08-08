#ifndef NUMBER_MATH_H
#define NUMBER_MATH_H

/* The real maths behind number_t.
 *
 * number_t is chosen at compile time (see config.h), and the library function
 * that operates on it changes with it: sqrt for double, sqrtl for long double,
 * sqrtq for __float128. C has no overloading to do that for us -- and the
 * complex maths under the formula parser is compiled as C, so <cmath> is not
 * available either -- so the choice is made once, here.
 *
 * Calling plain sqrt() on a number_t still compiles: the argument is converted
 * to double, the result comes back with 53 bits of mantissa, and nothing warns.
 * That is how every user-defined formula came to be evaluated in double
 * regardless of the build it ran in, which put a floor under how far one could
 * be zoomed into -- around 1e-15 -- while the built-in formulas, written in
 * number_t throughout, went as deep as the type allowed.
 */

#include "config.h"

#ifdef USE_FLOAT128

#include <quadmath.h>

#define nsqrt sqrtq
#define nexp expq
#define nlog logq
#define nlog10 log10q
#define nlog1p log1pq
#define nexpm1 expm1q
#define nsin sinq
#define ncos cosq
#define ntan tanq
#define nasin asinq
#define nacos acosq
#define natan atanq
#define natan2 atan2q
#define nsinh sinhq
#define ncosh coshq
#define ntanh tanhq
#define nasinh asinhq
#define nacosh acoshq
#define natanh atanhq
#define npow powq
#define nfabs fabsq
#define nfloor floorq
#define nceil ceilq
#define nhypot hypotq
#define nfmod fmodq
#define nround roundq
#define ntrunc truncq
#define ncopysign copysignq
#define nfrexp frexpq
#define nldexp ldexpq
#define nisfinite finiteq

#define N_PI M_PIq
#define N_PI_2 (M_PIq / 2)
#define N_E M_Eq

#else /* !USE_FLOAT128 */

#include <math.h>

#ifdef USE_LONG_DOUBLE

#define nsqrt sqrtl
#define nexp expl
#define nlog logl
#define nlog10 log10l
#define nlog1p log1pl
#define nexpm1 expm1l
#define nsin sinl
#define ncos cosl
#define ntan tanl
#define nasin asinl
#define nacos acosl
#define natan atanl
#define natan2 atan2l
#define nsinh sinhl
#define ncosh coshl
#define ntanh tanhl
#define nasinh asinhl
#define nacosh acoshl
#define natanh atanhl
#define npow powl
#define nfabs fabsl
#define nfloor floorl
#define nceil ceill
#define nhypot hypotl
#define nfmod fmodl
#define nround roundl
#define ntrunc truncl
#define ncopysign copysignl
#define nfrexp frexpl
#define nldexp ldexpl

#else /* plain double */

#define nsqrt sqrt
#define nexp exp
#define nlog log
#define nlog10 log10
#define nlog1p log1p
#define nexpm1 expm1
#define nsin sin
#define ncos cos
#define ntan tan
#define nasin asin
#define nacos acos
#define natan atan
#define natan2 atan2
#define nsinh sinh
#define ncosh cosh
#define ntanh tanh
#define nasinh asinh
#define nacosh acosh
#define natanh atanh
#define npow pow
#define nfabs fabs
#define nfloor floor
#define nceil ceil
#define nhypot hypot
#define nfmod fmod
#define nround round
#define ntrunc trunc
#define ncopysign copysign
#define nfrexp frexp
#define nldexp ldexp

#endif /* USE_LONG_DOUBLE */

#define nisfinite isfinite

/* Written out rather than taken from M_PIl and friends, which are a glibc
 * extension the toolchains XaoS builds with do not all provide. The digits
 * beyond what long double holds are harmless padding. */
#define N_PI 3.14159265358979323846264338327950288L
#define N_PI_2 1.57079632679489661923132169163975144L
#define N_E 2.71828182845904523536028747135266250L

#endif /* USE_FLOAT128 */

#define N_2PI (2 * N_PI)
#define N_1_2PI (1 / N_2PI)

#endif /* NUMBER_MATH_H */
