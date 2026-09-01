/*/////////////////////////////////////////////////////////////////////////////////////
// project : sFFe ( SegFault (or Segmentation Fault :) ) formula evalutaor )
// author  : Mateusz Malczak ( mateusz@malczak.info )
// wpage   : www.segfaultlabs.com/projects/sffe
///////////////////////////////////////////////////////////////////////////////////////
// special build for XaoS, for more info visit
// http://www.segfaultlabs.com/projects/sfXaos
/////////////////////////////////////////////////////////////////////////////////////*/

#ifdef SFFE_CMPLX_GSL

#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>

#include <cstdint>

#include "number_math.h"

/* Every entry is {implementation, argument count, name}, optionally followed by
 * a selector that marks the arguments as lazily evaluated.
 *
 * The first sffnctsfirst entries are the operators. sffe_function only searches
 * from there on, so an operator's spelling can never be reached as a name.
 *
 * Arguments are named below in the order they are written in a formula, f(a, b,
 * c) -- note that the sfaramN macros number them the other way round.
 *
 * A few entries do not compute what their name suggests. They are marked as
 * such: the behaviour stays because saved position files depend on it. */
const sffunction sfcmplxfunc[sffnctscount] = {
    /* --- operators, reached through sffe_operator rather than by name --- */
    {sfpow, 2, "^\0"}, /* a ^ b */
    {sfadd, 2, "+\0"}, /* a + b */
    {sfsub, 2, "-\0"}, /* a - b */
    {sfmul, 2, "*\0"}, /* a * b */
    {sfdiv, 2, "/\0"}, /* a / b */
    /* prefix minus; shares the '-' spelling with sfsub but takes one operand.
     * Not reachable by name: sffe_function only scans from sffnctsfirst. */
    {sfneg, 1, "-\0"}, /* -a */

    /* --- trigonometry over the complex plane --- */
    {sfsin, 1, "sin\0"},   /* nsin(a) */
    {sfcos, 1, "cos\0"},   /* ncos(a) */
    {sftan, 1, "tan\0"},   /* ntan(a) */
    {sfcot, 1, "cot\0"},   /* cot(a) */
    {sfasin, 1, "asin\0"}, /* arcsin(a) */
    {sfacos, 1, "acos\0"}, /* arccos(a) */
    {sfatan, 1, "atan\0"}, /* arctan(a) */
    {sfacot, 1, "acot\0"}, /* arccot(a) */
    /* natan2(y, x): angle of the real parts, plus i times the angle of the
     * imaginary parts. Two real arguments give the ordinary atan2. */
    {sfatan2, 2, "atan2\0"},

    /* --- hyperbolic --- */
    {sfsinh, 1, "sinh\0"}, /* nsinh(a) */
    {sfcosh, 1, "cosh\0"}, /* ncosh(a) */
    {sftanh, 1, "tanh\0"}, /* ntanh(a) */
    {sfcoth, 1, "coth\0"}, /* coth(a) */

    /* --- exponential and logarithms --- */
    {sfexp, 1, "exp\0"},     /* e ^ a */
    {sflog, 1, "log\0"},     /* natural logarithm of a */
    {sflog10, 1, "log10\0"}, /* logarithm of a in base 10 */
    {sflog2, 1, "log2\0"},   /* logarithm of a in base 2 */
    {sflogN, 2, "logn\0"},   /* logarithm of b in base a -- the base comes first */

    /* --- powers and roots --- */
    {sfpow, 2, "pow\0"},   /* a ^ b */
    {sfpowd, 2, "powd\0"}, /* a ^ real(b); the imaginary part of b is ignored */
    {sfsqr, 1, "sqr\0"}, /* a * a. This used to raise a to itself: sqr(3) was 27 */
    {sfsqrt, 1, "sqrt\0"}, /* principal square root of a */
    /* the c-th of the b b-th roots of a. This used to store that root over its
     * own first argument and evaluate to -1 rather than returning it. */
    {sfrtni, 3, "rtni\0"},
    /* 1 / a. The name used to end in \n rather than \0, which made it four
     * characters long, so it never matched a three-character lookup and the
     * function could not be called at all. */
    {sfinv, 1, "inv\0"},

    /* --- rounding, and pulling a value apart --- */
    {sfceil, 1, "ceil\0"},   /* ceiling of each component */
    {sffloor, 1, "floor\0"}, /* floor of each component */
    {sfabs, 1, "abs\0"},     /* |a|, as a real */
    {sfrabs, 1, "rabs\0"},   /* |real(a)|, as a real; the imaginary part is dropped */
    {sfre, 1, "re\0"},       /* real(a), as a real */
    {sfim, 1, "im\0"},       /* imag(a), as a real */
    {sfcarg, 1, "arg\0"},    /* angle of a, as a real */
    {sfmod, 1, "mod\0"},     /* fractional part of each component */
    {sfconj, 1, "conj\0"},   /* real(a) - i*imag(a); this used to swap them */

    /* --- burning ship variants --- */
    {sfbship, 1, "bship\0"},   /* |real(a)| + i*|imag(a)| */
    {sfbshipr, 1, "bshipr\0"}, /* |real(a)| + i*imag(a) */
    {sfbshipi, 1, "bshipi\0"}, /* real(a) + i*|imag(a)| */

    /* --- assembling one value out of two --- */
    {sfrect, 2, "rect\0"},   /* real(a) + i*imag(b) */
    {sfpolar, 2, "polar\0"}, /* |a| * e^(i*arg(b)) */

    /* --- smaller of two; the suffix says which part is compared --- */
    {sfmin, 2, "min\0"},   /* smaller real part and smaller imaginary part */
    {sfminr, 2, "minr\0"}, /* smaller real part; imaginary part taken from a */
    {sfmini, 2, "mini\0"}, /* smaller imaginary part; real part taken from a */
    {sfminm, 2, "minm\0"}, /* smaller modulus, held at the angle of a */

    /* --- larger of two, same convention --- */
    {sfmax, 2, "max\0"},   /* larger real part and larger imaginary part */
    {sfmaxr, 2, "maxr\0"}, /* larger real part; imaginary part taken from a */
    {sfmaxi, 2, "maxi\0"}, /* larger imaginary part; real part taken from a */
    {sfmaxm, 2, "maxm\0"}, /* larger modulus, held at the angle of a */

    /* --- mid(a, b, c) confines a to the range b..c.
     * With b < c it is a plain clamp. With b > c the range is inverted and a
     * value outside it is sent to the opposite end: mid(0,10,1) is 10 while
     * mid(99,10,1) is 1. --- */
    {sfmid, 3, "mid\0"},   /* both components confined */
    {sfmidr, 3, "midr\0"}, /* real part confined; imaginary part taken from a */
    {sfmidi, 3, "midi\0"}, /* imaginary part confined; real part taken from a */
    {sfmidm, 3, "midm\0"}, /* modulus confined, held at the angle of a */

    /* --- real trigonometry applied to each component separately --- */
    {sfsincos, 1, "sincos\0"}, /* nsin(real(a)) + i*ncos(imag(a)) */
    {sfcossin, 1, "cossin\0"}, /* ncos(real(a)) + i*nsin(imag(a)) */
    {sfsinr, 1, "sinr\0"},     /* nsin(real(a)); imaginary part passes through */
    {sfcosr, 1, "cosr\0"},     /* ncos(real(a)); imaginary part passes through */
    {sfsini, 1, "sini\0"},     /* nsin(imag(a)); real part passes through */
    {sfcosi, 1, "cosi\0"},     /* ncos(imag(a)); real part passes through */

    {sftancot, 1, "tancot\0"}, /* ntan(real(a)) + i*cot(imag(a)) */
    {sfcottan, 1, "cottan\0"}, /* cot(real(a)) + i*ntan(imag(a)) */
    {sftanr, 1, "tanr\0"},     /* ntan(real(a)); imaginary part passes through */
    {sfcotr, 1, "cotr\0"},     /* cot(real(a)); imaginary part passes through */
    {sftani, 1, "tani\0"},     /* ntan(imag(a)); real part passes through */
    {sfcoti, 1, "coti\0"},     /* cot(imag(a)); real part passes through */

    /* --- waveforms, again component by component --- */
    {sftrunc, 1, "trunc\0"},       /* each component truncated towards zero */
    {sfsawtooth, 1, "sawtooth\0"}, /* x - nfloor(x), a ramp in [0, 1) */
    {sftwave, 1, "twave\0"},       /* triangle wave of period 2, in [-1, 1] */

    /* --- assorted --- */
    {sfjulian, 3, "julian\0"}, /* |a|^b * e^(i*c*arg(a)) */
    /* inveps(a, b): real(a)/(|a|^2 + real(b)) - i*imag(a)/(|a|^2 + imag(b)),
     * an inverse softened by b so that it stays finite at the origin */
    {sfinveps, 2, "inveps\0"},
    /* atan2s(y, x): atan2 of each pair of components. Differs from atan2 only
     * on real arguments, where negating one gives a negative zero and
     * natan2(-0, -0) is -pi rather than 0. */
    {sfatan2s, 2, "atan2s\0"},

    /* ngon(a, b, c, d): folds a about the centre b onto a c-sided polygon,
     * the corner radius raised to the power d */
    {sfngon, 4, "ngon\0"},
    /* parchment(a, b): quantises the angle of a into |b| sectors, keeping |a| */
    {sfparchment, 2, "parchment\0"},
    /* parchmenta(a, b): as parchment, but mirroring alternate half sectors */
    {sfparchmenta, 2, "parchmenta\0"},

    /* --- snapping to a grid of step 1/n; n == 0 leaves the value alone --- */
    {sftruncv, 2, "truncv\0"},   /* both components, step 1/|b| */
    {sftruncc, 2, "truncc\0"},   /* real step 1/real(b), imaginary step 1/imag(b) */
    {sftruncvr, 2, "truncvr\0"}, /* real component only, step 1/|b| */
    {sftruncvi, 2, "truncvi\0"}, /* imaginary component only, step 1/|b| */
    {sftruncvm, 2, "truncvm\0"}, /* the modulus, angle kept */
    {sftruncva, 2, "truncva\0"}, /* the angle, modulus kept */

    /* gamma(a): Lanczos approximation of the complex gamma function.
     * Known defect: the series is scaled by nlog(nsqrt(2*pi)) where the formula
     * calls for nsqrt(2*pi), so every result is 0.3666 times the true gamma --
     * gamma(5) gives 8.798 instead of 24. */
    {sferf, 1, "erf\0"}, /* error function over the complex plane */
    {sfgamma, 1, "gamma\0"},
    {sflambertw, 1, "lambertw\0"}, /* principal branch of the Lambert W of a */

    /* iteration-dependent selection; -1 parameters means variadic, and the
     * selector marks the arguments as lazily evaluated */
    {sfifiter, SFFE_VARIADIC, "ifiter\0", sfifiter_sel},
    {sfifiterl, SFFE_VARIADIC, "ifiterl\0", sfifiterl_sel},

    /* Names with no implementation behind them. sffe_parse turns a call to one
     * of these into an unknown-function error rather than jumping through a
     * null pointer, which is what it used to do. */

    {sfrand, 1, "rand\0"}, /* real(a) times a random number in [0, 1) */
    /* Coherent noise over the position, 1 to 3 arguments and so
     * variadic. randsc interpolates and gives blobs; the rest do not and
     * give a mosaic of flat cells, differing only in how they cut the
     * plane up. See each of them for the rest. */
    {sfrandsc, SFFE_VARIADIC, "randsc\0"},
    {sfrandscq, SFFE_VARIADIC, "randscq\0"},
    {sfrandscp, SFFE_VARIADIC, "randscp\0"},
    {sfrandsch, SFFE_VARIADIC, "randsch\0"},
    {sfrandsct, SFFE_VARIADIC, "randsct\0"}};

const char sfcnames[sfvarscount][6] = {"pi\0", "pi_2\0", "pi2\0",
                                       "e\0",  "i\0",    "rnd\0"};

const cfptr sfcvals[sfvarscount] = {sfcPI, sfcPI2, sfc2PI, sfcE, sfcI, sfcRND};

sfarg *sfadd(sfarg *const p)
{ /* + */
    sfvalue(p) = gsl_complex_add(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfsub(sfarg *const p)
{ /* - */
    sfvalue(p) = gsl_complex_sub(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfneg(sfarg *const p)
{ /* unary - */
    sfvalue(p) = gsl_complex_negative(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

/* Current iteration of the fractal loop, 0 on the first one. Maintained by the
 * engine (see formulas.cpp); 0 for anyone evaluating a formula outside it. */
thread_local unsigned int sffe_iteration = 0;

/* The selectors run before any argument has been evaluated, to tell the
 * evaluator which one to bother computing. They must therefore depend only on
 * the iteration, never on the arguments -- which is exactly what these two do.
 * The functions below then read the value that was actually produced. */
unsigned int sfifiter_sel(unsigned int argc)
{ /* ifiter: cycles through its arguments */
    return sffe_iteration % argc;
}

unsigned int sfifiterl_sel(unsigned int argc)
{ /* ifiterl: stays on the last argument once past it */
    return sffe_iteration < argc ? sffe_iteration : argc - 1;
}

/* args are held right to left, so source index i sits at args[argc - 1 - i] */
sfarg *sfifiter(sfarg *const p)
{
    sfvalue(p) = sfvalue(p->args[p->argc - 1 - sfifiter_sel(p->argc)]);
    return p;
}

sfarg *sfifiterl(sfarg *const p)
{
    sfvalue(p) = sfvalue(p->args[p->argc - 1 - sfifiterl_sel(p->argc)]);
    return p;
}

sfarg *sfmul(sfarg *const p)
{ /* *  */
    sfvalue(p) = gsl_complex_mul(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfdiv(sfarg *const p)
{ /*  /   */
    sfvalue(p) = gsl_complex_div(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfsin(sfarg *const p)
{ /* sin */
    sfvalue(p) = gsl_complex_sin(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcos(sfarg *const p)
{ /* cos */
    sfvalue(p) = gsl_complex_cos(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sftan(sfarg *const p)
{ /* tan */
    sfvalue(p) = gsl_complex_tan(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcot(sfarg *const p)
{ /* ctan */
    sfvalue(p) = gsl_complex_cot(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfasin(sfarg *const p)
{ /* asin */
    sfvalue(p) = gsl_complex_arcsin(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfacos(sfarg *const p)
{ /* acos */
    sfvalue(p) = gsl_complex_arccos(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfatan(sfarg *const p)
{ /* atan */
    sfvalue(p) = gsl_complex_arctan(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfacot(sfarg *const p)
{ /* actan */
    sfvalue(p) = gsl_complex_arccot(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfatan2(sfarg *const p)
{ /* natan2(y, x) = angle of the real parts + i * angle of the imaginary parts */
    sfNumber y = sfvalue(sfaram2(p));
    sfNumber x = sfvalue(sfaram1(p));

    number_t hor = natan2(GSL_REAL(y), GSL_REAL(x));

    /* With real arguments the answer has to be the plain atan2, so a pair of
     * zeros contributes nothing. Left to atan2 they would not: negating a real
     * number gives a negative zero, and natan2(-0, -0) is -pi, not 0. */
    number_t ver = 0.0;
    if (GSL_IMAG(y) != 0 || GSL_IMAG(x) != 0) {
        ver = natan2(GSL_IMAG(y), GSL_IMAG(x));
    }

    cmplxset(sfvalue(p), hor, ver);
    return p;
}

sfarg *sfsinh(sfarg *const p)
{ /* sinh */
    sfvalue(p) = gsl_complex_sinh(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcosh(sfarg *const p)
{ /* cosh */
    sfvalue(p) = gsl_complex_cosh(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sftanh(sfarg *const p)
{ /* tanh */
    sfvalue(p) = gsl_complex_tanh(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcoth(sfarg *const p)
{ /* ctanh */
    sfvalue(p) = gsl_complex_coth(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfexp(sfarg *const p)
{ /* exp */
    sfvalue(p) = gsl_complex_exp(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sflog(sfarg *const p)
{ /* log */
    sfvalue(p) = gsl_complex_log(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sflog10(sfarg *const p)
{ /* log10 */
    sfvalue(p) = gsl_complex_log10(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sflog2(sfarg *const p)
{ /* log2 */
    sfNumber base;
    real(base) = 2;
    imag(base) = 0;
    sfvalue(p) = gsl_complex_log_b(sfvalue(sfaram1(p)), base);
    return sfaram1(p);
}

sfarg *sflogN(sfarg *const p)
{ /* logN */
    sfvalue(p) = gsl_complex_log_b(sfvalue(sfaram1(p)), sfvalue(sfaram2(p)));
    return sfaram2(p);
}

sfarg *sfpow(sfarg *const p)
{ /* cmplx pow */
    sfvalue(p) = gsl_complex_pow(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfpowd(sfarg *const p)
{ /* int pow */
    sfvalue(p) = gsl_complex_pow_real(sfvalue(sfaram2(p)),
                                      GSL_REAL(sfvalue(sfaram1(p))));
    return sfaram2(p);
}

sfarg *sfsqr(sfarg *const p)
{ /* sqr: a squared.
   *
   * This used to compute gsl_complex_pow(a, a), a raised to itself, so sqr(3)
   * answered 27. A multiplication is both the right answer and far cheaper
   * than a complex power, which goes through a log and an exp. */
    sfvalue(p) = gsl_complex_mul(sfvalue(sfaram1(p)), sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfsqrt(sfarg *const p)
{ /* sqrt */
    sfvalue(p) = gsl_complex_sqrt(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfrtni(sfarg *const p)
{ /* rtni(a, b, c): the c-th of the b b-th roots of a.
   *
   * Arguments are (z, n, i) as written, so sfaram3 is z and sfaram1 is i.
   *
   * This used to store the root over its own first argument and evaluate to
   * -1, rather than returning it. Writing through that argument reached the
   * caller's variable, so rtni(z,...) silently redefined z for the rest of the
   * formula: rtni(z,12,6)+z answered -2.0595 instead of 0.9405, because the
   * second z read the root as well. */
    number_t n = (number_t)(int)real(sfvalue(sfaram2(p)));
    number_t nrz = npow(gsl_complex_abs(sfvalue(sfaram3(p))), 1.0 / n);
    number_t alfi = (gsl_complex_arg(sfvalue(sfaram3(p))) +
                   8 * natan(1.0) * (number_t)(int)real(sfvalue(sfaram1(p)))) /
                  n;

    cmplxset(sfvalue(p), nrz * ncos(alfi), nrz * nsin(alfi));
    return p;
}

sfarg *sfinv(sfarg *const p)
{ /* cinv */
    sfvalue(p) = gsl_complex_inverse(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfceil(sfarg *const p)
{ /* ceil */
    // sfvalue(p) = nceil( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = nceil(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = nceil(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sffloor(sfarg *const p)
{ /* floor */
    // sfvalue(p) = nfloor( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = nfloor(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = nfloor(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcarg(sfarg *const p)
{ /* floor */
    // sfvalue(p) = nfloor( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = gsl_complex_arg(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = 0.0;
    return sfaram1(p);
}

sfarg *sfmod(sfarg *const p)
{ /* floor */
    // sfvalue(p) = nfloor( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = nfmod(GSL_REAL(sfvalue(sfaram1(p))), 1);
    GSL_IMAG(sfvalue(p)) = nfmod(GSL_IMAG(sfvalue(sfaram1(p))), 1);
    return sfaram1(p);
}

sfarg *sfconj(sfarg *const p)
{ /* conj: real(a) - i*imag(a).
   *
   * This used to swap the two components instead of negating the imaginary
   * one, so conj(3+4i) answered 4+3i. */
    sfvalue(p) = gsl_complex_conjugate(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfabs(sfarg *const p)
{ /* abs - |z| */
    GSL_REAL(sfvalue(p)) = gsl_complex_abs(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = 0.0;
    return sfaram1(p);
}

sfarg *sfrabs(sfarg *const p)
{ /* abs - real numbers */
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    if (GSL_REAL(sfvalue(p)) < 0)
        GSL_REAL(sfvalue(p)) = -GSL_REAL(sfvalue(p));
    GSL_IMAG(sfvalue(p)) = 0;
    return sfaram1(p);
}

sfarg *sfre(sfarg *const p)
{ /* RE */
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = 0.0;
    return sfaram1(p);
}

sfarg *sfim(sfarg *const p)
{ /* IM */
    GSL_REAL(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = 0.0;
    return sfaram1(p);
}

sfarg *sfrand(sfarg *const p)
{ /* rand */
    GSL_REAL(sfvalue(p)) =
        GSL_REAL(sfvalue(sfaram1(p))) * (number_t)rand() / (number_t)RAND_MAX;
    GSL_IMAG(sfvalue(p)) = 0;
    return sfaram1(p);
}

sfarg *sfbship(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = abs(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = abs(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfbshipr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = abs(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfbshipi(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = abs(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

/* rect(a, b) and polar(a, b) each build one value from two, taking a part from
 * each. Both used to read them the wrong way round -- rect answered
 * real(b) + i*imag(a) -- which is the opposite of what their declarations in
 * sffe_cmplx_gsl.h have always said. Remember that sfaram1 is the argument
 * written last, so a is sfaram2 and b is sfaram1. */

sfarg *sfrect(sfarg *const p)
{ /* rect(a, b) = real(a) + i*imag(b) */
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfpolar(sfarg *const p)
{ /* polar(a, b) = |a| * e^(i*arg(b)) */
    number_t radius = gsl_complex_abs(sfvalue(sfaram2(p)));
    number_t theta = gsl_complex_arg(sfvalue(sfaram1(p)));
    sfvalue(p) = gsl_complex_polar(radius, theta);
    return sfaram2(p);
}

sfarg *sfmax(sfarg *const p)
{
    number_t r1 = GSL_REAL(sfvalue(sfaram2(p)));
    number_t r2 = GSL_REAL(sfvalue(sfaram1(p)));

    number_t i1 = GSL_IMAG(sfvalue(sfaram2(p)));
    number_t i2 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = r1 < r2 ? r2 : r1;
    GSL_IMAG(sfvalue(p)) = i1 < i2 ? i2 : i1;
    return sfaram2(p);
}

sfarg *sfmaxr(sfarg *const p)
{
    number_t r1 = GSL_REAL(sfvalue(sfaram2(p)));
    number_t r2 = GSL_REAL(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = r1 < r2 ? r2 : r1;
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    return sfaram2(p);
}

sfarg *sfmaxi(sfarg *const p)
{
    number_t i1 = GSL_IMAG(sfvalue(sfaram2(p)));
    number_t i2 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
    GSL_IMAG(sfvalue(p)) = i1 < i2 ? i2 : i1;
    return sfaram2(p);
}

sfarg *sfmaxm(sfarg *const p)
{
    number_t r1 = gsl_complex_abs(sfvalue(sfaram2(p)));
    number_t r2 = gsl_complex_abs(sfvalue(sfaram1(p)));
    number_t theta = gsl_complex_arg(sfvalue(sfaram2(p)));

    sfvalue(p) = gsl_complex_polar(r1 < r2 ? r2 : r1, theta);
    return sfaram2(p);
}

sfarg *sfmin(sfarg *const p)
{
    number_t r1 = GSL_REAL(sfvalue(sfaram2(p)));
    number_t r2 = GSL_REAL(sfvalue(sfaram1(p)));

    number_t i1 = GSL_IMAG(sfvalue(sfaram2(p)));
    number_t i2 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = r1 < r2 ? r1 : r2;
    GSL_IMAG(sfvalue(p)) = i1 < i2 ? i1 : i2;
    return sfaram2(p);
}

sfarg *sfminr(sfarg *const p)
{
    number_t r1 = GSL_REAL(sfvalue(sfaram2(p)));
    number_t r2 = GSL_REAL(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = r1 < r2 ? r1 : r2;
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    return sfaram2(p);
}

sfarg *sfmini(sfarg *const p)
{
    number_t i1 = GSL_IMAG(sfvalue(sfaram2(p)));
    number_t i2 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
    GSL_IMAG(sfvalue(p)) = i1 < i2 ? i1 : i2;
    return sfaram2(p);
}

sfarg *sfminm(sfarg *const p)
{
    number_t r1 = gsl_complex_abs(sfvalue(sfaram2(p)));
    number_t r2 = gsl_complex_abs(sfvalue(sfaram1(p)));
    number_t theta = gsl_complex_arg(sfvalue(sfaram2(p)));

    sfvalue(p) = gsl_complex_polar(r1 < r2 ? r1 : r2, theta);
    return sfaram2(p);
}

number_t calc_mid(number_t v1, number_t v2, number_t v3) {
    if (v2 < v3) {
        if (v1 < v2) {
            return v2;
        } else if (v1 > v3) {
            return v3;
        } else {
            return v1;
        }
    } else {
        if (v1 < v2 && v1 < v3) {
            return v2;
        } else if (v1 > v3 && v1 > v2) {
            return v3;
        } else {
            return v1;
        }
    }
}

sfarg *sfmid(sfarg *const p)
{
    number_t r1 = GSL_REAL(sfvalue(sfaram3(p)));
    number_t r2 = GSL_REAL(sfvalue(sfaram2(p)));
    number_t r3 = GSL_REAL(sfvalue(sfaram1(p)));

    number_t i1 = GSL_IMAG(sfvalue(sfaram3(p)));
    number_t i2 = GSL_IMAG(sfvalue(sfaram2(p)));
    number_t i3 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = calc_mid(r1, r2, r3);
    GSL_IMAG(sfvalue(p)) = calc_mid(i1, i2, i3);
    return sfaram3(p);
}

sfarg *sfmidr(sfarg *const p)
{
    number_t r1 = GSL_REAL(sfvalue(sfaram3(p)));
    number_t r2 = GSL_REAL(sfvalue(sfaram2(p)));
    number_t r3 = GSL_REAL(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = calc_mid(r1, r2, r3);
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram3(p)));
    return sfaram3(p);
}

sfarg *sfmidi(sfarg *const p)
{
    number_t i1 = GSL_IMAG(sfvalue(sfaram3(p)));
    number_t i2 = GSL_IMAG(sfvalue(sfaram2(p)));
    number_t i3 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram3(p)));
    GSL_IMAG(sfvalue(p)) = calc_mid(i1, i2, i3);
    return sfaram3(p);
}

sfarg *sfmidm(sfarg *const p)
{
    number_t r1 = gsl_complex_abs(sfvalue(sfaram3(p)));
    number_t r2 = gsl_complex_abs(sfvalue(sfaram2(p)));
    number_t r3 = gsl_complex_abs(sfvalue(sfaram1(p)));
    number_t theta = gsl_complex_arg(sfvalue(sfaram3(p)));

    sfvalue(p) = gsl_complex_polar(calc_mid(r1, r2, r3), theta);
    return sfaram3(p);
}

sfarg *sfsincos(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = nsin(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = ncos(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcossin(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = ncos(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = nsin(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfsinr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = nsin(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcosr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = ncos(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfsini(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = nsin(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcosi(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = ncos(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

number_t cot(number_t x) {
    return ncos(x)/nsin(x);
}

sfarg *sftancot(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = ntan(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = cot(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcottan(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = cot(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = ntan(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sftanr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = ntan(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcotr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = cot(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sftani(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = ntan(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcoti(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = cot(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sftrunc(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = ntrunc(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = ntrunc(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

number_t sawtooth(number_t x) {
    return x - nfloor(x);
}

sfarg *sfsawtooth(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = sawtooth(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = sawtooth(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

number_t twave(number_t x) {
    number_t xf = x/2.0;
    return 2.0*abs(2.0*(xf-nfloor(xf+0.5)))-1.0;
}

sfarg *sftwave(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = twave(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = twave(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfjulian(sfarg *const p)
{
    gsl_complex z = sfvalue(sfaram3(p));
    gsl_complex m;
    GSL_SET_COMPLEX(&m, gsl_complex_abs(z), 0);
    m = gsl_complex_pow(m, sfvalue(sfaram2(p)));
    gsl_complex b = sfvalue(sfaram1(p));
    number_t mx = GSL_REAL(m);
    number_t my = GSL_IMAG(m);
    number_t arg = gsl_complex_arg(z);
    number_t byg = nexp(-GSL_IMAG(b)*arg);
    number_t bxg = arg * GSL_REAL(b);
    number_t cosbxg = ncos(bxg);
    number_t sinbxg = nsin(bxg);

    GSL_REAL(sfvalue(p)) = byg*(mx*cosbxg - my*sinbxg);
    GSL_IMAG(sfvalue(p)) = byg*(my*cosbxg + mx*sinbxg);
    return sfaram3(p);
}

sfarg *sfinveps(sfarg *const p)
{ /* cinv */
    number_t x = GSL_REAL(sfvalue(sfaram2(p)));
    number_t y = GSL_IMAG(sfvalue(sfaram2(p)));
    number_t delta = (x*x + y*y);
    GSL_REAL(sfvalue(p)) = x/(delta + GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = -y/(delta + GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram2(p);
}

sfarg *sfatan2s(sfarg *const p)
{ /* cinv */
    GSL_REAL(sfvalue(p)) = natan2(GSL_REAL(sfvalue(sfaram2(p))), GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = natan2(GSL_IMAG(sfvalue(sfaram2(p))), GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram2(p);
}


sfarg *sfngon(sfarg *const p)
{
    gsl_complex i;
    GSL_SET_COMPLEX(&i, 0.0, 1.0);

    gsl_complex n = sfvalue(sfaram2(p));
    gsl_complex zc = gsl_complex_sub(sfvalue(sfaram4(p)), sfvalue(sfaram3(p)));
    number_t t = gsl_complex_arg(zc);
    gsl_complex tn = gsl_complex_mul_real(n, t * N_1_2PI);
    tn = gsl_complex_add_real(tn, 0.5);
    GSL_REAL(tn) = nfloor(GSL_REAL(tn));
    GSL_IMAG(tn) = nfloor(GSL_IMAG(tn));
    tn = gsl_complex_mul_real(tn, N_2PI);
    tn = gsl_complex_div(tn, n);
    number_t cr = ncos(t);
    number_t sr = nsin(t);
    gsl_complex ccn = gsl_complex_cos(tn);
    gsl_complex scn = gsl_complex_sin(tn);
    gsl_complex rn = gsl_complex_add(gsl_complex_mul_real(ccn, cr),
                                     gsl_complex_mul_real(scn, sr));
    rn = gsl_complex_mul_real(gsl_complex_pow(rn, sfvalue(sfaram1(p))),
                              gsl_complex_abs(zc));
    gsl_complex argexp = gsl_complex_exp(gsl_complex_mul_real(i, t));
    sfvalue(p) =  gsl_complex_add(gsl_complex_mul(rn, argexp), sfvalue(sfaram3(p)));

    return sfaram4(p);
}

sfarg *sfparchment(sfarg *const p)
{
    gsl_complex z = sfvalue(sfaram2(p));
    number_t n = gsl_complex_abs(sfvalue(sfaram1(p)));
    //if (n == 2 && abs(GSL_REAL(z) - (-1.5)) < 0.1  && abs(GSL_IMAG(z) - (-1)) < 0.1) {
    //    int vb = 1;
    //}

    number_t t = gsl_complex_arg(z);
    number_t dN = n * N_1_2PI;
    number_t nN = 1/dN;

    number_t trc = nceil(t * dN) * nN;

    number_t trm = t - trc + nN;

    sfvalue(p) = gsl_complex_polar(gsl_complex_abs(z), trm);
    return sfaram2(p);
}

sfarg *sfparchmenta(sfarg *const p)
{
    gsl_complex z = sfvalue(sfaram2(p));
    number_t n = gsl_complex_abs(sfvalue(sfaram1(p)));
    //if (n == 5 && abs(GSL_REAL(z) - (-1.5)) < 0.1  && abs(GSL_IMAG(z) - (-1)) < 0.1) {
    //    int vb = 1;
    //}

    number_t t = gsl_complex_arg(z);
    number_t dN = n * N_1_2PI;
    number_t nN = 1/dN;
    number_t trc = nceil(t * dN) * nN;

    dN = dN*2;
    nN = 1/dN;
    number_t trc2 = nceil(t * dN) * nN;

    number_t trm = trc < trc2 + 0.1 / n ? trc2 - t : t + nN - trc2;

    sfvalue(p) = gsl_complex_polar(gsl_complex_abs(z), trm);
    return sfaram2(p);
}

sfarg *sftruncv(sfarg *const p)
{
    number_t n = gsl_complex_abs(sfvalue(sfaram1(p)));

    if (n != 0) {
        GSL_REAL(sfvalue(p)) = ntrunc(GSL_REAL(sfvalue(sfaram2(p))) * n) / n;
        GSL_IMAG(sfvalue(p)) = ntrunc(GSL_IMAG(sfvalue(sfaram2(p))) * n) / n;
    } else {
        GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
        GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    }
    return sfaram2(p);
}

sfarg *sftruncc(sfarg *const p)
{
    number_t nr = GSL_REAL(sfvalue(sfaram1(p)));
    number_t ni = GSL_IMAG(sfvalue(sfaram1(p)));

    if (nr != 0) {
        GSL_REAL(sfvalue(p)) = ntrunc(GSL_REAL(sfvalue(sfaram2(p))) * nr) / nr;
    } else {
        GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
    }

    if (ni != 0) {
        GSL_IMAG(sfvalue(p)) = ntrunc(GSL_IMAG(sfvalue(sfaram2(p))) * ni) / ni;
    } else {
        GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    }
    return sfaram2(p);
}

sfarg *sftruncvr(sfarg *const p)
{
    number_t n = gsl_complex_abs(sfvalue(sfaram1(p)));

    if (n != 0) {
        GSL_REAL(sfvalue(p)) = ntrunc(GSL_REAL(sfvalue(sfaram2(p))) * n) / n;
    } else {
        GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
    }
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    return sfaram2(p);
}

sfarg *sftruncvi(sfarg *const p)
{
    number_t n = gsl_complex_abs(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
    if (n != 0) {
        GSL_IMAG(sfvalue(p)) = ntrunc(GSL_IMAG(sfvalue(sfaram2(p))) * n) / n;
    } else {
        GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    }
    return sfaram2(p);
}

sfarg *sftruncvm(sfarg *const p)
{
    number_t n = gsl_complex_abs(sfvalue(sfaram1(p)));
    number_t m;
    number_t a;

    if (n != 0) {
        m = ntrunc(gsl_complex_abs(sfvalue(sfaram2(p))) * n) / n;
    } else {
        m = gsl_complex_abs(sfvalue(sfaram2(p)));
    }
    a = gsl_complex_arg(sfvalue(sfaram2(p)));
    sfvalue(p) = gsl_complex_polar(m, a);
    return sfaram2(p);
}

sfarg *sftruncva(sfarg *const p)
{
    number_t n = gsl_complex_abs(sfvalue(sfaram1(p)));
    number_t m;
    number_t a;

    m = gsl_complex_abs(sfvalue(sfaram2(p)));
    if (n != 0) {
        a = ntrunc(gsl_complex_arg(sfvalue(sfaram2(p))) * n) / n;
    } else {
        a = gsl_complex_arg(sfvalue(sfaram2(p)));
    }
    sfvalue(p) = gsl_complex_polar(m, a);
    return sfaram2(p);
}

/* Lanczos coefficients for g = 7 with 9 terms, good to about 15 digits.
 * Kept here rather than in the header: nothing else needs them, and a static
 * array in a header gets a private copy in every translation unit. */
static const int LANCZOS_G = 7;
static const number_t LANCZOS_P[9] = {0.99999999999980993,
                                    676.5203681218851,
                                    -1259.1392167224028,
                                    771.32342877765313,
                                    -176.61502916214059,
                                    12.507343278686905,
                                    -0.13857109526572012,
                                    9.9843695780195716e-6,
                                    1.5056327351493116e-7};

#define SQRT_2PI 2.5066282746310005024157652848110

/**
 * @brief Complex gamma function via the Lanczos approximation.
 * @details Poles at the non-positive integers give NaN. For Re(z) < 0.5, where
 * the series does not apply, the reflection formula is used instead; that
 * recurses exactly once, since 1 - z then has real part above 0.5.
 * @param z The argument.
 * @return Gamma(z).
 */
gsl_complex complex_gamma_lanczos(gsl_complex z)
{
    number_t zr = GSL_REAL(z);
    number_t zi = GSL_IMAG(z);
    gsl_complex temp;

    /* poles at 0, -1, -2, ... */
    if (zi == 0.0 && zr <= 0.0 && nfloor(zr) == zr) {
        GSL_SET_COMPLEX(&temp, NAN, NAN);
        return temp;
    }

    /* Gamma(z) = pi / (nsin(pi z) * Gamma(1 - z)) */
    if (zr < 0.5) {
        gsl_complex one_minus_z;
        GSL_SET_COMPLEX(&one_minus_z, 1.0 - zr, -zi);
        gsl_complex den =
            gsl_complex_mul(gsl_complex_sin(gsl_complex_mul_real(z, N_PI)),
                            complex_gamma_lanczos(one_minus_z));
        GSL_SET_COMPLEX(&temp, N_PI, 0.0);
        return gsl_complex_div(temp, den);
    }

    z = gsl_complex_sub_real(z, 1.0);
    zr = GSL_REAL(z);
    zi = GSL_IMAG(z);

    /* x = p[0] + sum p[i] / (z + i).
     *
     * Every numerator is real, so each term is p[i] * conj(z + i) / |z + i|^2:
     * one division and a handful of multiplications, instead of the general
     * complex division this used to go through. */
    number_t xr = LANCZOS_P[0];
    number_t xi = 0.0;
    for (int i = 1; i < 9; ++i) {
        number_t dr = zr + (number_t)i;
        number_t s = LANCZOS_P[i] / (dr * dr + zi * zi);
        xr += s * dr;
        xi -= s * zi;
    }
    gsl_complex x;
    GSL_SET_COMPLEX(&x, xr, xi);

    /* Gamma(z + 1) = nsqrt(2 pi) * t^(z + 0.5) * e^-t * x, with t = z + g + 0.5.
     *
     * Folding the power and the exponential into a single nexp((z + 0.5) *
     * nlog(t) - t) saves one complex exponential: t^(z + 0.5) is itself an
     * nexp(nlog(...)) underneath.
     *
     * The scale factor is nsqrt(2 pi). It used to be nlog(nsqrt(2 pi)), which
     * left every result multiplied by 0.3666 -- gamma(5) came out as 8.798
     * rather than 24. */
    gsl_complex t = gsl_complex_add_real(z, LANCZOS_G + 0.5);
    gsl_complex e = gsl_complex_sub(
        gsl_complex_mul(gsl_complex_add_real(z, 0.5), gsl_complex_log(t)), t);

    return gsl_complex_mul_real(gsl_complex_mul(gsl_complex_exp(e), x),
                                SQRT_2PI);
}

/**
 * @brief sFFe wrapper for the complex gamma function.
 * @param p The call; its argument is sfaram1(p).
 * @return Pointer to the input argument, per the sffe convention.
 */
#ifdef USE_FLOAT128
/* The Lanczos approximation above is a fixed set of coefficients, chosen for
 * about fifteen significant digits, and no amount of arithmetic behind it does
 * better: measured against Gamma(n+1) = n!, which is exact, it holds 2.2e-15
 * whether it runs at 64 or 113 bits of mantissa. That is the whole of a
 * double's precision and roughly a five-thousandth of a quad's.
 *
 * So the quad build uses Stirling instead, which has no such ceiling: shift
 * the argument up by the recurrence until it is large enough for the
 * asymptotic series to converge, then take as many Bernoulli terms as the type
 * needs. The series error is bounded by its first omitted term, and at
 * Re(w) >= 50 the term after the last one kept here is below 1e-37.
 *
 * The coefficients are B(2k) / (2k (2k-1)), written as the exact integer
 * ratios they are so that the division happens at whatever precision number_t
 * has, rather than being transcribed as decimals that would reintroduce the
 * ceiling this exists to remove. Every numerator fits a double exactly, so
 * none of them is rounded on the way in. */
static const number_t STIRLING_NUM[15] = {
    1, -1, 1, -1, 1, -691, 1, -3617, 43867, -174611,
    854513, -236364091, 8553103, -23749461029, 8615841276005};
static const number_t STIRLING_DEN[15] = {
    12, 360, 1260, 1680, 1188, 360360, 156, 122400, 244188, 125400,
    63756, 1506960, 3900, 657720, 12460140};

static gsl_complex complex_gamma_stirling(gsl_complex z)
{
    /* The series needs Re(w) large and positive, so the left half plane comes
     * back through the reflection formula. One level deep only: 1 - z has real
     * part at least 0.5 whenever z's is below it. */
    /* GSL_COMPLEX_ONE and gsl_complex_rect are only defined under HAVE_INLINE,
     * which this build does not set, so constants are built the way the rest
     * of this file builds them. */
    gsl_complex one, pi;
    GSL_SET_COMPLEX(&one, 1, 0);
    GSL_SET_COMPLEX(&pi, N_PI, 0);

    if (GSL_REAL(z) < 0.5) {
        gsl_complex s = gsl_complex_sin(gsl_complex_mul_real(z, N_PI));
        gsl_complex g = complex_gamma_stirling(gsl_complex_sub(one, z));
        return gsl_complex_div(pi, gsl_complex_mul(s, g));
    }

    /* Gamma(z) = Gamma(z + n) / (z (z+1) ... (z+n-1)) */
    gsl_complex w = z;
    gsl_complex denominator = one;
    while (GSL_REAL(w) < 50) {
        denominator = gsl_complex_mul(denominator, w);
        w = gsl_complex_add_real(w, 1.0);
    }

    /* log Gamma(w) = (w - 1/2) log w - w + log(2 pi)/2 + sum B(2k) / (2k(2k-1)
     * w^(2k-1)) */
    gsl_complex lg = gsl_complex_sub(
        gsl_complex_mul(gsl_complex_add_real(w, -0.5), gsl_complex_log(w)), w);
    lg = gsl_complex_add_real(lg, nlog((number_t)2 * N_PI) / 2);

    gsl_complex term = gsl_complex_inverse(w);
    gsl_complex step = gsl_complex_mul(term, term);
    for (int k = 0; k < 15; ++k) {
        lg = gsl_complex_add(lg, gsl_complex_mul_real(term, STIRLING_NUM[k] /
                                                                STIRLING_DEN[k]));
        term = gsl_complex_mul(term, step);
    }

    return gsl_complex_div(gsl_complex_exp(lg), denominator);
}
#endif

/* 2/sqrt(pi), derived rather than written out: a decimal literal is a double
 * whatever it is assigned to, which would cap erf at sixteen digits in the
 * quad build for no reason at all. */
static const number_t ERF_2_SQRTPI = (number_t)2 / nsqrt(N_PI);

#ifdef USE_FLOAT128
#define ERF_TOLERANCE ((number_t)1e-35)
#define ERF_MAX_TERMS 200
#else
#define ERF_TOLERANCE ((number_t)1e-21)
#define ERF_MAX_TERMS 120
#endif

/**
 * @brief Error function over the complex plane.
 * @details The Maclaurin series
 *
 *     erf(z) = 2/sqrt(pi) * sum (-1)^n z^(2n+1) / (n! (2n+1))
 *
 * evaluated by recurrence: each term is the one before it times -z^2/n, so a
 * term costs one complex multiplication and one real division, with no powers
 * and no factorials computed at all. It stops as soon as a term can no longer
 * move the sum, which for the values a fractal actually iterates -- inside a
 * bailout of two or four -- happens after a dozen or so terms. The odd
 * symmetry erf(-z) = -erf(z) folds the left half plane onto the right one,
 * which costs a sign test and halves the ground to cover.
 *
 * One expansion and not two, deliberately. The usual companion is the
 * continued fraction for erfc, which is excellent far out along the real axis
 * and poor near the imaginary one, where erf grows like e^(|z|^2); measured
 * across every crossover tried, the two disagreed by between 1e-5 and 1 at the
 * boundary. A seam of that size is a visible edge in a rendered fractal, which
 * is a worse fault than the one it would fix -- the series holds to a few ulp
 * wherever a fractal goes, and only loses digits out past a modulus of four or
 * five, which is beyond any bailout the iteration would have stopped at.
 *
 * @param p The call; its argument is sfaram1(p).
 * @return Pointer to the input argument, per the sffe convention.
 */
sfarg *sferf(sfarg *const p)
{
    gsl_complex z = sfvalue(sfaram1(p));
    int negate = 0;

    if (GSL_REAL(z) < 0) {
        z = gsl_complex_negative(z);
        negate = 1;
    }

    gsl_complex minus_z2 = gsl_complex_negative(gsl_complex_mul(z, z));
    gsl_complex term = z;
    gsl_complex r = z;
    for (int n = 1; n < ERF_MAX_TERMS; n++) {
        term =
            gsl_complex_div_real(gsl_complex_mul(term, minus_z2), (number_t)n);
        gsl_complex add = gsl_complex_div_real(term, (number_t)(2 * n + 1));
        r = gsl_complex_add(r, add);
        if (gsl_complex_abs2(add) <=
            ERF_TOLERANCE * ERF_TOLERANCE * gsl_complex_abs2(r))
            break;
    }
    r = gsl_complex_mul_real(r, ERF_2_SQRTPI);

    if (negate)
        r = gsl_complex_negative(r);
    sfvalue(p) = r;
    return sfaram1(p);
}

/* The point being iterated, set by the engine before each formula is
 * evaluated. It lives here rather than being read from formulas.cpp so that
 * the parser stays linkable on its own -- the test binaries build it without
 * the engine -- and so that a test can place a point directly.
 *
 * The position and not z: z diverges between a long double and a quad build
 * after enough iterations, by construction, so anything hashed from it differs
 * between the two. The position is computed from the view in a few operations
 * that do not amplify, and agrees to the last bit or two. */
thread_local cmplx sffe_position = {{0, 0}};

/* A hash of integers, which is the same at any precision: the alternative,
 * hashing a float through sin(), depends on the exact sin() of the build and
 * would give the two binaries different pictures. */
static uint64_t randsc_hash(int64_t i, int64_t j, uint64_t seed)
{
    uint64_t x = (uint64_t)i * 0x9E3779B97F4A7C15ULL ^
                 (uint64_t)j * 0xC2B2AE3D27D4EB4FULL ^ seed;
    x ^= x >> 33;
    x *= 0xFF51AFD7ED558CCDULL;
    x ^= x >> 33;
    x *= 0xC4CEB9FE1A85EC53ULL;
    x ^= x >> 33;
    return x;
}

static number_t randsc_unit(uint64_t x)
{
    return (number_t)(x >> 11) / (number_t)((uint64_t)1 << 53);
}

/* A real seed has to survive being written once and read by two builds:
 * "0.525" lands just below the exact value at long double and just above it at
 * quad, so the two differ around 1e-20 and a hash of them shares nothing.
 * Keeping the leading 40 bits puts both in the same bucket -- the quantum,
 * about 9e-13, is seven orders coarser than the disagreement. An integer seed
 * needs none of this and is exact, which is the reason to prefer one. */
/* What randsc_setup returns.
 *
 * RANDSC_BEYOND says the cells have been degraded finer than the grid can
 * address, so the cell handed back is the one saturated index rather than a
 * position in the plane. Each function then hashes it with its own salt and
 * stops: there is nothing for a tiling to do with a single cell, and the
 * transforms would only run the same numbers off the scale again.
 *
 * The field is flat there, and that is not a choice made here -- a cell a
 * millionth of a pixel wide has no picture left to give. What matters is that
 * the five functions still differ from one another, since a formula
 * subtracting one from another would otherwise reach exactly zero and iterate
 * to the limit for nothing. */
#define RANDSC_STOP 0
#define RANDSC_OK 1
#define RANDSC_BEYOND 2

/* How far a grid coordinate may go before it stops being one.
 *
 * Every one of these functions ends by hashing a pair of integers, so a
 * position divided by a cell size has to fit in one. It stops fitting when the
 * degradation has shrunk the cells far enough -- 0.3 * 0.5^n against a
 * position of order one runs out at about the sixty-fifth pass -- and what
 * happened then was the worst of both: the conversion overflowed, every point
 * in the plane landed on the same saturated index, and the field quietly went
 * flat while costing five times as much to compute, since arithmetic on
 * astronomical values is the slow path of the library floor.
 *
 * So the size running out of integers is treated as the size running out
 * altogether, which the code already had a case for. The limit is a quarter of
 * what an int64 holds, which leaves room for the coordinate transforms the
 * hexagons and the triangles apply on top -- at most a factor of 1.6 -- and
 * costs nothing: a cell that fine is far below any pixel that could show it.
 */
#define RANDSC_INDEX_LIMIT ((number_t)2.0e18)

/* floor(x) as an integer index, with the fraction left over, or 0 if the
 * value is not one an index can hold. NaN fails the comparison and is
 * rejected with the rest.
 *
 * This is floorl followed by a conversion, without the library call: the
 * conversion truncates toward zero, which is floor for anything not negative
 * and one too many otherwise. x - (number_t)c is exact, the two being within
 * one of each other. Measured at a fifth of what floorl costs, and a
 * twentieth of what floorl costs on the values the guard now refuses. */
static int randsc_cell(number_t x, int64_t *cell, number_t *frac)
{
    if (!(x > -RANDSC_INDEX_LIMIT && x < RANDSC_INDEX_LIMIT))
        return 0;
    int64_t c = (int64_t)x;
    number_t t = (number_t)c;
    if (t > x) {
        c--;
        t -= 1;
    }
    *cell = c;
    *frac = x - t;
    return 1;
}

/* roundl without the library call, half away from zero, exact for anything
 * the guard above lets through: the truncation and the remainder are both
 * exact there, so the comparison against a half decides it. */
static number_t randsc_round(number_t x)
{
    number_t t = (number_t)(int64_t)x;
    number_t d = x - t;
    if (d >= (number_t)0.5)
        return t + 1;
    if (d <= (number_t)-0.5)
        return t - 1;
    return t;
}

static const number_t RANDSC_SEED_SCALE = nldexp((number_t)1, 40);

static uint64_t randsc_seed(cmplx seed)
{
    int64_t a = 0, b = 0;
    number_t unused;
    randsc_cell(GSL_REAL(seed) * RANDSC_SEED_SCALE, &a, &unused);
    randsc_cell(GSL_IMAG(seed) * RANDSC_SEED_SCALE, &b, &unused);
    return randsc_hash(a, b, 0x5DEECE66DULL);
}


/* Shared by randsc and randscq: reads the arguments, applies the degradation
 * for the iteration reached, and returns the cell the position falls in along
 * with where inside it. Zero means the caller should not compute -- a zero in
 * either component of either argument divides by zero once the degradation
 * gets there.
 *
 * The arguments are read right to left, so which is which depends on how many
 * were given; see sfaramN. */
/* Forced, not suggested. This is the preamble of all five functions and the
 * compiler kept it out of line, which meant writing six values to memory --
 * two cells, two fractions, a hash and a state -- and reading them straight
 * back. Inlined it costs the formula a twentieth less. */
#if defined(__GNUC__)
#define RANDSC_INLINE inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#define RANDSC_INLINE __forceinline
#else
#define RANDSC_INLINE inline
#endif

static RANDSC_INLINE int randsc_setup(sfarg *const p, int64_t *cx, int64_t *cy, number_t *u,
                        number_t *v, uint64_t *hash)
{
    cmplx size, degradation, seed;
    GSL_SET_COMPLEX(&size, 1, 1);
    GSL_SET_COMPLEX(&degradation, 1, 1);
    GSL_SET_COMPLEX(&seed, 0, 0);

    switch (p->argc) {
        case 3:
            degradation = sfvalue(sfaram1(p));
            size = sfvalue(sfaram2(p));
            seed = sfvalue(sfaram3(p));
            break;
        case 2:
            size = sfvalue(sfaram1(p));
            seed = sfvalue(sfaram2(p));
            break;
        case 1:
            seed = sfvalue(sfaram1(p));
            break;
        default:
            return RANDSC_STOP;
    }

    if (GSL_REAL(size) == 0 || GSL_IMAG(size) == 0 ||
        GSL_REAL(degradation) == 0 || GSL_IMAG(degradation) == 0)
        return RANDSC_STOP;

    /* size * degradation^n, per component, carried from one pass to the next
     * rather than worked out again.
     *
     * This is what degradation says it does: the size is multiplied by it at
     * every pass. Computing the power instead grew with the iteration limit,
     * since squaring takes a step per bit of n -- 5.6 ns at pass sixteen and
     * 48 at sixteen thousand, against a flat 3.5 for one multiplication -- and
     * a formula with three noise calls paid it three times a pass. The plain
     * loop before that was worse still, being a step per pass.
     *
     * The running product lives on the call site (see sfarg), so two calls in
     * one formula keep their own, and a thread cannot disturb another. A pass
     * earlier than the one carried means a new pixel has started, which is
     * where the count goes back to nothing.
     *
     * Reaching a pass by multiplying once per pass gives the same answer
     * whatever route was taken to get there, so the picture does not depend on
     * the order the pixels were computed in. It is not the same answer that
     * squaring gives -- the two associate the multiplications differently, and
     * differ by some three parts in 10^18. */
    if (p->carried == 0 || p->carried > sffe_iteration) {
        GSL_SET_COMPLEX(&p->carry, 1, 1);
        p->carried = 0;
    }
    if (GSL_REAL(degradation) == 1 && GSL_IMAG(degradation) == 1) {
        /* The default, and multiplying by one leaves the product where it is,
         * so the passes can be counted off without doing any of them. */
        p->carried = sffe_iteration;
    } else {
        while (p->carried < sffe_iteration) {
            GSL_SET_COMPLEX(&p->carry,
                            GSL_REAL(p->carry) * GSL_REAL(degradation),
                            GSL_IMAG(p->carry) * GSL_IMAG(degradation));
            p->carried++;
        }
    }
    number_t wr = GSL_REAL(size) * GSL_REAL(p->carry);
    number_t wi = GSL_IMAG(size) * GSL_IMAG(p->carry);
    if (wr == 0 || wi == 0) /* shrunk past what the type can hold */
        return RANDSC_STOP;

    /* The iteration goes into the hash, not only into the size above.
     * Without it the field is fixed once the size is: a degradation of one
     * never changes the size, so every pass returned the very same value at
     * the same point, and a degradation near one -- 0.99, say -- moved the
     * grid by a percent and returned very nearly it. Hashing the iteration
     * gives each pass a field of its own, which is what a formula asks for
     * when it calls this once per iteration. Space is untouched by it: for
     * a fixed pass the hash is a constant, so the noise is as coherent from
     * point to point as it ever was. */
    *hash = randsc_hash((int64_t)sffe_iteration, 0, randsc_seed(seed));
    if (!randsc_cell(GSL_REAL(sffe_position) / wr, cx, u) ||
        !randsc_cell(GSL_IMAG(sffe_position) / wi, cy, v)) {
        /* Past the resolution of the grid. Not an error and not a refusal:
         * the caller gets one flat cell over the whole plane, which is what
         * it got before by accident, and gets it for the price of a
         * comparison. See RANDSC_BEYOND. */
        *cx = *cy = INT64_MIN;
        *u = *v = 0;
        return RANDSC_BEYOND;
    }
    return RANDSC_OK;
}

/**
 * @brief Coherent noise over the position, seeded and reproducible.
 * @details randsc(seed), randsc(seed; size), randsc(seed; size; degradation).
 *
 * Value noise: the plane is cut into cells, each corner is hashed to a number,
 * and the value between them is interpolated with a smooth curve. Nearby
 * points therefore give nearby values -- blobs rather than the per-pixel snow
 * a plain hash gives -- and that continuity is also what makes the result
 * stable. A difference in the input produces a difference of the same order in
 * the output, so the two precisions agree to about 1e-19 where a raw hash
 * would agree not at all. It is also why the cell boundary is not visible:
 * leaving one cell with weight 1 gives the same corner value as entering the
 * next with weight 0.
 *
 * size, default 1+i, is the average width of a blob along the real axis and
 * its height along the imaginary one.
 *
 * degradation, default 1+i, shrinks the blobs as the iteration proceeds: the
 * size in force is size * degradation^n, taken component by component, so
 * degradation 0.5+0.2i with size 1+i gives 1+i on the first pass, then
 * 0.5+0.2i, then 0.25+0.04i. Component by component and not as a complex
 * power, which would give 0.21+0.2i for the third.
 *
 * A zero in either component of either argument would divide by zero once the
 * degradation reached it, so the function returns zero instead of computing.
 *
 * The result is real, in [0, 1), with the imaginary part left at zero, as rand
 * does. Two independent fields are two calls with different seeds.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
sfarg *sfrandsc(sfarg *const p)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p),
                        randsc_unit(randsc_hash(cx, cy, h)), 0);
        return sfaram1(p);
    }

    u = u * u * (3 - 2 * u); /* smoothstep: flat at both ends, so the value */
    v = v * v * (3 - 2 * v); /* meets its neighbour without a crease */

    number_t a = randsc_unit(randsc_hash(cx, cy, h));
    number_t b = randsc_unit(randsc_hash(cx + 1, cy, h));
    number_t c = randsc_unit(randsc_hash(cx, cy + 1, h));
    number_t d = randsc_unit(randsc_hash(cx + 1, cy + 1, h));
    number_t lo = a + (b - a) * u;
    number_t hi = c + (d - c) * u;

    GSL_SET_COMPLEX(&sfvalue(p), lo + (hi - lo) * v, 0);
    return sfaram1(p);
}

/* Where the seed of a cell may sit, as a fraction of the cell: the middle
 * 60%. See sfrandscp for why it is not the whole cell. */
#define RANDSCP_JITTER_LOW ((number_t)0.2)
#define RANDSCP_JITTER_SPAN ((number_t)0.6)

/* A 32-bit field turned into [0, 1), for the two halves of one hash. */
static number_t randsc_unit32(uint32_t x)
{
    return (number_t)x / (number_t)((uint64_t)1 << 32);
}

/* The splitmix64 finalizer, to get a value out of a hash already spent on
 * placing a seed without the two being related. */
static uint64_t randsc_remix(uint64_t x)
{
    x ^= x >> 30;
    x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27;
    x *= 0x94D049BB133111EBULL;
    x ^= x >> 31;
    return x;
}

/**
 * @brief The same field again with straight edges: irregular flat polygons.
 * @details randscp takes the arguments randsc takes and means the same by
 * them. Where randsc smooths between the corners of a cell and randscq leaves
 * the cell itself flat, this scatters one point inside each cell and gives
 * every position the value of the nearest of them.
 *
 * What that draws is a Voronoi diagram: the boundary between two neighbouring
 * points is the perpendicular bisector of the segment joining them, so every
 * cell is a convex polygon with straight sides, and no two are the same shape.
 * It is the same underlying grid as randsc -- so size and degradation stretch
 * and shrink it identically -- with the curves taken out and the regularity of
 * randscq's mosaic taken out with them.
 *
 * The nearest point is looked for in the nine cells around the one the point
 * falls in. That is enough only because the scatter is held to the middle 60%
 * of each cell: the furthest a point can be from its own cell's seed is then
 * sqrt(2) * 0.8 = 1.13 cells, while the closest seed two cells away is 1.2, so
 * nothing outside the nine can win. Letting the seeds reach the cell edges
 * would be more irregular and occasionally wrong.
 *
 * One hash per cell serves for both coordinates of its seed -- the two halves
 * of a 64-bit mix are independent enough -- and the winner is mixed once more
 * for the value, so that a cell's colour does not follow where its seed sits.
 *
 * Being a step function it shares randscq's caveat: two builds that place a
 * point on opposite sides of an edge return unrelated values, so they agree
 * everywhere but a hairline along the edges.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
sfarg *sfrandscp(sfarg *const p)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p),
                        randsc_unit(randsc_remix(randsc_hash(cx, cy, h))), 0);
        return sfaram1(p);
    }

    /* Larger than any distance the nine cells can produce. */
    number_t bestd = 16;
    uint64_t besth = 0;

    for (int j = -1; j <= 1; j++)
        for (int i = -1; i <= 1; i++) {
            uint64_t hh = randsc_hash(cx + i, cy + j, h);
            number_t fx = RANDSCP_JITTER_LOW +
                          RANDSCP_JITTER_SPAN * randsc_unit32((uint32_t)hh);
            number_t fy =
                RANDSCP_JITTER_LOW +
                RANDSCP_JITTER_SPAN * randsc_unit32((uint32_t)(hh >> 32));
            number_t dx = (number_t)i + fx - u;
            number_t dy = (number_t)j + fy - v;
            number_t d = dx * dx + dy * dy;
            if (d < bestd) {
                bestd = d;
                besth = hh;
            }
        }

    GSL_SET_COMPLEX(&sfvalue(p), randsc_unit(randsc_remix(besth)), 0);
    return sfaram1(p);
}

/**
 * @brief The same field without the interpolation: a mosaic of flat cells.
 * @details randscq takes the arguments randsc takes and means the same by
 * them, but returns the value of the cell the point falls in rather than a
 * blend of the four around it. Each cell is therefore one flat colour and the
 * result is a grid of squares, size wide and size high, instead of blobs.
 *
 * What is gained in look is paid for in stability, and it is worth being plain
 * about it. randsc is continuous, so a difference of 1e-20 between a long
 * double and a quad build moves the result by 1e-20. randscq is a step
 * function: two builds that place a point on opposite sides of a cell edge
 * return values with nothing to do with each other. The two binaries therefore
 * agree everywhere except on a hairline along the cell edges, and a zoom that
 * reuses a row whose coordinate has drifted within its tolerance can flicker
 * there. That is inherent in asking for hard edges, not a defect to be fixed.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
/* Whether a parsed formula calls randsc or randscq.
 *
 * The engine asks because boundary tracing has to be turned off for such a
 * formula: it walks the edge of a region, finds one colour all the way round,
 * and fills the inside without computing it. That holds for a fractal, whose
 * bands really are solid, and not for a noise field, where the inside is
 * whatever the noise says. Left on, some pixels are filled rather than
 * computed and the noise is simply wrong there.
 *
 * Walking the operation list rather than searching the text: the text would
 * also match a name that merely contains "randsc", and this cannot. */
int sffe_uses_noise(sffe *const parser)
{
    if (parser == NULL)
        return 0;
    for (unsigned int i = 0; i < parser->oprCount; i++)
        if (parser->oprs[i].fnc == sfrandsc ||
            parser->oprs[i].fnc == sfrandscq ||
            parser->oprs[i].fnc == sfrandscp ||
            parser->oprs[i].fnc == sfrandsch ||
            parser->oprs[i].fnc == sfrandsct)
            return 1;
    return 0;
}

sfarg *sfrandscq(sfarg *const p)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p),
                        randsc_unit(randsc_hash(cx, cy, h)), 0);
        return sfaram1(p);
    }

    GSL_SET_COMPLEX(&sfvalue(p), randsc_unit(randsc_hash(cx, cy, h)), 0);
    return sfaram1(p);
}

/* sqrt(3), for the two tilings whose cells are not axis-aligned. Worked out
 * once at the precision in use: a decimal literal would be a double and would
 * cap both grids at sixteen digits in the quad build. */
static const number_t RANDSC_SQRT3 = nsqrt((number_t)3);

/* size means the same thing across the whole family, so that changing one
 * letter of a formula changes the shape of the cells and not their scale.
 *
 * randsc, randscq and randscp all lay one cell over each unit square of the
 * degraded size, so their cells have unit area. A hexagon of circumradius one
 * has an area of 3*sqrt(3)/2, near enough 2.6, and an equilateral triangle of
 * side one has sqrt(3)/4, near enough 0.43; laid out as they come, the two
 * would be six times apart from each other and both wrong against the rest.
 * These factors scale each grid to unit cells: the reciprocal of the square
 * root of the area of one cell of a lattice of pitch one. */
static const number_t RANDSCH_PITCH = nsqrt(3 * nsqrt((number_t)3) / 2);
static const number_t RANDSCT_PITCH = nsqrt(nsqrt((number_t)3) / 4);

/* A salt apiece, so that a cell index landing on the same pair of integers in
 * two tilings does not hand back the same value in both. Without them the
 * hexagons and the triangles agreed with the squares wherever the indices
 * happened to meet, which is often enough near the origin. */
#define RANDSCH_SALT 0xD1B54A32D192ED03ULL
#define RANDSCT_SALT 0x8CB92BA72F3D8DD7ULL
#define RANDSCT_UPPER 0xA5A5A5A5A5A5A5A5ULL

/**
 * @brief The same mosaic on a hexagonal grid: a honeycomb of flat cells.
 * @details randsch takes the arguments randsc takes and means the same by
 * them. It is randscq with the squares replaced by hexagons.
 *
 * Of the three regular polygons that tile the plane, squares are randscq
 * already and triangles alternate in orientation, which reads as a pattern
 * rather than a texture. A honeycomb has no such grain: every cell has the
 * same six neighbours at the same six angles, so it looks less like a grid and
 * more like a material.
 *
 * The point is taken to the axial coordinates of a pointy-topped grid and
 * rounded through cube coordinates -- three numbers summing to zero, of which
 * the one that moved furthest is given whatever the other two leave over.
 * Rounding the two axial numbers on their own would land on the nearest
 * rhombus of the grid, which is not the nearest hexagon.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
sfarg *sfrandsch(sfarg *const p)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p),
                        randsc_unit(randsc_hash(cx, cy, h ^ RANDSCH_SALT)), 0);
        return sfaram1(p);
    }

    /* randsc_setup hands back the square cell and where the point sits inside
     * it; adding them recovers the position in units of the degraded size,
     * which is what the hexagon grid is laid out in. */
    number_t X = ((number_t)cx + u) * RANDSCH_PITCH;
    number_t Y = ((number_t)cy + v) * RANDSCH_PITCH;

    number_t q = RANDSC_SQRT3 / 3 * X - (number_t)1 / 3 * Y;
    number_t r = (number_t)2 / 3 * Y;

    number_t ax = q, az = r, ay = -q - r;
    number_t rx = randsc_round(ax), ry = randsc_round(ay),
             rz = randsc_round(az);
    number_t dx = nfabs(rx - ax), dy = nfabs(ry - ay), dz = nfabs(rz - az);
    if (dx > dy && dx > dz)
        rx = -ry - rz;
    else if (dy > dz)
        ry = -rx - rz;
    else
        rz = -rx - ry;

    GSL_SET_COMPLEX(&sfvalue(p),
                    randsc_unit(randsc_hash((int64_t)rx, (int64_t)rz,
                                            h ^ RANDSCH_SALT)),
                    0);
    return sfaram1(p);
}

/**
 * @brief The same mosaic on a triangular grid: flat equilateral triangles.
 * @details randsct takes the arguments randsc takes and means the same by
 * them. It is randscq with the squares replaced by triangles.
 *
 * The triangular tiling is the rhombic one cut in half. Taking the point to
 * the basis (1, 0) and (1/2, sqrt(3)/2) gives a lattice of rhombi with every
 * side 1; the diagonal joining the two far corners of a rhombus is also 1, so
 * it cuts the rhombus into two equilateral triangles. Which side of that
 * diagonal the point falls on is the sum of its two fractional coordinates
 * against one, and that bit goes into the hash along with the rhombus.
 *
 * Unlike the hexagons the two orientations alternate, which is a property of
 * the tiling and not of the noise: a triangular mosaic has a grain, and that
 * is what one asks for by choosing it.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
sfarg *sfrandsct(sfarg *const p)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p),
                        randsc_unit(randsc_hash(cx, cy, h ^ RANDSCT_SALT)), 0);
        return sfaram1(p);
    }

    number_t X = ((number_t)cx + u) * RANDSCT_PITCH;
    number_t Y = ((number_t)cy + v) * RANDSCT_PITCH;

    /* The inverse of the basis above: the rhombus is the unit square of these
     * two coordinates. */
    number_t a = X - Y / RANDSC_SQRT3;
    number_t b = (number_t)2 * Y / RANDSC_SQRT3;
    int64_t ia, ib;
    number_t fa, fb;
    randsc_cell(a, &ia, &fa);
    randsc_cell(b, &ib, &fb);
    int upper = fa + fb >= 1;

    GSL_SET_COMPLEX(
        &sfvalue(p),
        randsc_unit(randsc_hash(ia, ib,
                                upper ? h ^ RANDSCT_SALT ^ RANDSCT_UPPER
                                      : h ^ RANDSCT_SALT)),
        0);
    return sfaram1(p);
}

sfarg *sfgamma(sfarg *const p)
{
#ifdef USE_FLOAT128
    sfvalue(p) = complex_gamma_stirling(sfvalue(sfaram1(p)));
#else
    sfvalue(p) = complex_gamma_lanczos(sfvalue(sfaram1(p)));
#endif
    return sfaram1(p);
}


/**
 * @brief Principal branch W_0 of the complex Lambert W function.
 * @details W solves w * e^w = z. GSL has no complex version, so this iterates
 * Halley's method, which converges cubically.
 *
 * The starting point is picked by region, because a single formula does not
 * work everywhere: the series around the origin diverges past |z| = 1/e, and
 * the asymptotic expansion needs log z to be well away from zero.
 *
 * This used to run Householder's order-3 method from w = log z. That has a
 * stationary point: at z = 1 the start is w = 0, where the update's numerator
 * 1 + l2*l1/2 vanishes exactly, so the step was zero and the iteration
 * returned the starting guess. lambertw(1) came out as 0 instead of 0.5671.
 * Halley's update has no such cancellation, and costs less per step.
 *
 * @param p The call; its argument is sfaram1(p).
 * @return Pointer to the input argument, per the sffe convention.
 */
sfarg *sflambertw(sfarg *const p)
{
    /* Cubic convergence from these starting points settles in four or five
     * steps; the cap only guards inputs that do not converge at all. */
    const int MAX_ITERATIONS = 20;
    const number_t TOLERANCE = 1e-15;

    gsl_complex z = sfvalue(sfaram1(p));
    gsl_complex w;

    if (GSL_REAL(z) == 0.0 && GSL_IMAG(z) == 0.0) {
        sfvalue(p) = z; /* W(0) = 0, and the iteration cannot start there */
        return sfaram1(p);
    }

    /* The starting point decides which branch Halley converges to, so these
     * three regions are about correctness and not only about speed. */
    gsl_complex ez1 = gsl_complex_add_real(gsl_complex_mul_real(z, N_E), 1.0);
    if (gsl_complex_abs2(ez1) < 1.0) {
        /* Around the branch point z = -1/e the two real branches meet, and an
         * ordinary guess slides onto W_-1 -- for z = -0.3 that answers -1.7813
         * instead of -0.4894. The expansion there is in p = nsqrt(2(e z + 1)),
         * where W_0 takes +p and W_-1 would take -p. */
        gsl_complex q = gsl_complex_sqrt(gsl_complex_mul_real(ez1, 2.0));
        gsl_complex q2 = gsl_complex_mul(q, q);
        w = gsl_complex_add_real(q, -1.0);
        w = gsl_complex_sub(w, gsl_complex_mul_real(q2, 1.0 / 3.0));
        w = gsl_complex_add(
            w, gsl_complex_mul_real(gsl_complex_mul(q2, q), 11.0 / 72.0));
    } else if (gsl_complex_abs2(z) < 0.1296) {
        /* |z| < 0.36, just inside the radius 1/e of the series about the
         * origin, W = z - z^2 + 3z^3/2. Used any further out it lands in the
         * wrong basin: at z = -0.67 - 0.02i it used to answer W_-2. The test
         * comes after the branch point and not before it, because log z is
         * large here and would send the asymptotic form below astray. */
        gsl_complex z2 = gsl_complex_mul(z, z);
        w = gsl_complex_add(gsl_complex_sub(z, z2),
                            gsl_complex_mul_real(gsl_complex_mul(z2, z), 1.5));
    } else if ((gsl_complex_abs2(gsl_complex_log(z)) > 4.0 &&
                gsl_complex_abs2(z) > 0.81) ||
               gsl_complex_abs2(gsl_complex_add_real(z, 1.0)) < 0.25) {
        /* W = L1 - L2 + L2/L1 + ..., L1 = log z, L2 = log L1; keeping the
         * third term is what saves the iterations.
         *
         * Both halves of the test earn their place. Selecting on |log z| alone
         * lets in points of small modulus but large argument, where the
         * expansion does not hold -- z = -0.30 - 0.44i started at -0.76 +
         * 0.25i and converged to a far branch. Selecting on |z| alone lets in
         * z near 1, where log z is nearly zero and L2 diverges: z = 1.05
         * started at -63. The last clause covers z near -1, where the ratio
         * below would divide by nearly nothing, and where this expansion
         * happens to be good anyway. */
        gsl_complex l1 = gsl_complex_log(z);
        gsl_complex l2 = gsl_complex_log(l1);
        w = gsl_complex_add(gsl_complex_sub(l1, l2), gsl_complex_div(l2, l1));
    } else {
        /* A band roughly between |z| = 0.36 and 1, plus the region around
         * z = 1 the asymptotic form cannot serve. W is of order one there and
         * z/(1+z) is the right size throughout. */
        w = gsl_complex_div(z, gsl_complex_add_real(z, 1.0));
    }

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        gsl_complex ew = gsl_complex_exp(w);
        gsl_complex f = gsl_complex_sub(gsl_complex_mul(w, ew), z);

        /* Halley: step = f / (f' - f*f''/(2f')), which for f = w e^w - z is
         * e^w (w+1) - (w+2) f / (2w+2). The exponential cancels out of the
         * second derivative ratio, so no extra exp is needed. */
        gsl_complex wp1 = gsl_complex_add_real(w, 1.0);
        gsl_complex den = gsl_complex_sub(
            gsl_complex_mul(ew, wp1),
            gsl_complex_div(gsl_complex_mul(gsl_complex_add_real(w, 2.0), f),
                            gsl_complex_mul_real(wp1, 2.0)));

        if (gsl_complex_abs2(den) < 1e-300) {
            break;
        }

        gsl_complex step = gsl_complex_div(f, den);
        w = gsl_complex_sub(w, step);

        /* Squared moduli, to keep hypot out of the inner loop. */
        if (gsl_complex_abs2(step) <=
            TOLERANCE * TOLERANCE * (gsl_complex_abs2(w) + 1e-300)) {
            break;
        }
    }

    sfvalue(p) = w;
    return sfaram1(p);
}



// const eval
void sfcPI(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, 4 * natan(1), 0); }

void sfcPI2(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, 2 * natan(1), 0); }

void sfc2PI(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, 8 * natan(1), 0); }

void sfcE(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, nexp(1), 0); }

void sfcI(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, 0, 1); }

void sfcRND(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, rand(), 0); }

#endif
