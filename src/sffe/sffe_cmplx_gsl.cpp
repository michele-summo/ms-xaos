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
#include <cstring>

#include "number_math.h"
#include "randsctile_tables.h"

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
    {sfjulian, SFFE_VARIADIC, "julian\0", NULL, false, 2}, /* |a|^b * e^(i*c*arg(a)) */
    /* inveps(a, b): real(a)/(|a|^2 + real(b)) - i*imag(a)/(|a|^2 + imag(b)),
     * an inverse softened by b so that it stays finite at the origin */
    {sfinveps, SFFE_VARIADIC, "inveps\0", NULL, false, 2},
    /* atan2s(y, x): atan2 of each pair of components. Differs from atan2 only
     * on real arguments, where negating one gives a negative zero and
     * natan2(-0, -0) is -pi rather than 0. */
    {sfatan2s, 2, "atan2s\0"},

    /* ngon(a, b, c, d): folds a about the centre b onto a c-sided polygon,
     * the corner radius raised to the power d */
    {sfngon, SFFE_VARIADIC, "ngon\0", NULL, false, 2},
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
    {sfifiter, SFFE_VARIADIC, "ifiter\0", sfifiter_sel, false, 2},
    {sfifiterl, SFFE_VARIADIC, "ifiterl\0", sfifiterl_sel, false, 2},
    {sfifiterf, 2, "ifiterf\0", sfifiterf_sel},
    /* ifiterr's threshold is the argument the selector reads: the parser
     * evaluates that one first and then chooses between the two before it. */
    {sfifiterr, 3, "ifiterr\0", sfifiterr_sel, true},

    /* Names with no implementation behind them. sffe_parse turns a call to one
     * of these into an unknown-function error rather than jumping through a
     * null pointer, which is what it used to do. */

    {sfrand, 1, "rand\0"}, /* real(a) times a random number in [0, 1) */
    /* Coherent noise over the position, 1 to 3 arguments and so
     * variadic. randsc interpolates and gives blobs; the rest do not and
     * give a mosaic of flat cells, differing only in how they cut the
     * plane up. See each of them for the rest. */
    {sfrandsc, SFFE_VARIADIC, "randsc\0", NULL, false, 2},
    {sfrandscq, SFFE_VARIADIC, "randscq\0", NULL, false, 2},
    {sfrandscp, SFFE_VARIADIC, "randscp\0", NULL, false, 2},
    {sfrandsch, SFFE_VARIADIC, "randsch\0", NULL, false, 2},
    {sfrandsct, SFFE_VARIADIC, "randsct\0", NULL, false, 2},
    /* The same field over a tiling its first argument chooses among
     * forty-five, the seed second; so the first two are needed. */
    {sfrandsctile, SFFE_VARIADIC, "randsctile\0", NULL, false, 3},
    /* fbm(value, seed, ...): octaves of the same noise summed over a point
     * the caller names, rather than over the position. See sffbm. */
    {sffbm, SFFE_VARIADIC, "fbm\0", NULL, false, 3},

    /* Watching the orbit rather than the point: one number about the whole of
     * it, handed back on the last pass. 1 to 4 arguments and so variadic. */
    {sftrap, SFFE_VARIADIC, "trap\0", NULL, false, 2},
    {sfstripe, SFFE_VARIADIC, "stripe\0", NULL, false, 2},

    /* Figures rather than noise: the same field over the position, drawn
     * instead of diced. Every argument has a default, so "snowflake()" is a
     * call and so is "sierpinskyc( ,5)". */
    {sfsierpinskyt, SFFE_VARIADIC, "sierpinskyt\0", NULL, false, 1},
    {sfsierpinskyc, SFFE_VARIADIC, "sierpinskyc\0", NULL, false, 1},
    {sfsnowflake, SFFE_VARIADIC, "snowflake\0", NULL, false, 1},

    /* A polynomial in the first argument, the rest being its coefficients
     * from the highest power down. */
    {sfpoly, SFFE_VARIADIC, "poly\0", NULL, false, 2}};

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
unsigned int sfifiter_sel(unsigned int argc, const sfNumber *)
{ /* ifiter: cycles through its arguments */
    return sffe_iteration % argc;
}

unsigned int sfifiterl_sel(unsigned int argc, const sfNumber *)
{ /* ifiterl: stays on the last argument once past it */
    return sffe_iteration < argc ? sffe_iteration : argc - 1;
}

/* How many passes the picture allows, which is what "the last one" means.
 *
 * A formula cannot know it is on its final pass by escaping -- that is decided
 * by the value it is about to produce -- so the last pass has to mean the last
 * one the iteration limit allows. The engine sets this beside sffe_iteration
 * at the head of every pixel. Zero says nobody has set it, and then ifiterf
 * simply never fires, which is what a test or a bare parser wants. */
thread_local unsigned int sffe_maxiter = 0;

unsigned int sfifiterf_sel(unsigned int argc, const sfNumber *)
{ /* ifiterf: the last argument on the final pass, the first on all the rest */
    return (sffe_maxiter && sffe_iteration + 1 >= sffe_maxiter) ? argc - 1 : 0;
}

unsigned int sfifiterr_sel(unsigned int argc, const sfNumber *probe)
{ /* ifiterr: the second branch once the passes reach the threshold, which the
   * parser has evaluated for us before asking */
    (void)argc;
    return (probe && (number_t)sffe_iteration >= GSL_REAL(*probe)) ? 1 : 0;
}

/* args are held right to left, so source index i sits at args[argc - 1 - i] */
sfarg *sfifiter(sfarg *const p)
{
    sfvalue(p) = sfvalue(p->args[p->argc - 1 - sfifiter_sel(p->argc, NULL)]);
    return p;
}

sfarg *sfifiterl(sfarg *const p)
{
    sfvalue(p) = sfvalue(p->args[p->argc - 1 - sfifiterl_sel(p->argc, NULL)]);
    return p;
}

/**
 * @brief The second formula on the last pass, the first on every other.
 * @details ifiterf(a, b) evaluates a on every pass but the final one, and b on
 * that. The final pass is the last the iteration limit allows, since a formula
 * has no way of knowing which pass will be the one that escapes -- that
 * depends on the value it has not produced yet.
 *
 * Only the chosen one is evaluated, as with ifiter: the selector is consulted
 * before either argument has run.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the first argument written, unused by the evaluator.
 */
sfarg *sfifiterf(sfarg *const p)
{
    sfvalue(p) = sfvalue(p->args[p->argc - 1 - sfifiterf_sel(p->argc, NULL)]);
    return p;
}

/**
 * @brief The second formula once the passes reach a count, the first before.
 * @details ifiterr(a, b, n) evaluates a while the pass number is below n and b
 * from n onwards. n is read as a real number and may be any expression.
 *
 * Only the chosen branch is evaluated, as with ifiter, though the threshold
 * had to be taught to the parser first: the lazy mechanism chose a block
 * before any argument had run, and so could not consult a threshold that is
 * itself an argument. An argument may now be marked as read by the selector
 * rather than chosen by it, in which case it is evaluated first and the
 * branches are the arguments before it.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the first argument written, unused by the evaluator.
 */
sfarg *sfifiterr(sfarg *const p)
{
    /* The same question the selector was asked, and so the same answer. The
     * branch that was not chosen has not been evaluated and holds whatever it
     * last held, so it must not be read. */
    unsigned int k = sfifiterr_sel(p->argc - 1, sfaram1(p)->value);
    sfvalue(p) = sfvalue(p->args[p->argc - 1 - k]);
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

/* The argument in a given place, counting from the first one written, or the
 * value it takes when the call does not give it -- either by stopping short of
 * it, "julian(z)", or by leaving its place empty, "julian(z, ,3)".
 *
 * A call with defaults has to be variadic -- the parser counts what it is
 * given and hands the count over -- and the sfaramN macros count from the
 * other end, which is no use when the end moves. */
static inline cmplx sfarg_or(sfarg *const p, unsigned int place, number_t re,
                             number_t im)
{
    if (place <= p->argc && !p->args[p->argc - place]->omitted)
        return sfvalue(p->args[p->argc - place]);
    cmplx fallback;
    GSL_SET_COMPLEX(&fallback, re, im);
    return fallback;
}

sfarg *sfjulian(sfarg *const p)
{
    /* the modulus raised to the first, the angle multiplied by the first:
     * julian(a) is a itself */
    gsl_complex z = sfarg_or(p, 1, 0, 0);
    gsl_complex m;
    GSL_SET_COMPLEX(&m, gsl_complex_abs(z), 0);
    m = gsl_complex_pow(m, sfarg_or(p, 2, 1, 0));
    gsl_complex b = sfarg_or(p, 3, 1, 0);
    number_t mx = GSL_REAL(m);
    number_t my = GSL_IMAG(m);
    number_t arg = gsl_complex_arg(z);
    number_t byg = nexp(-GSL_IMAG(b)*arg);
    number_t bxg = arg * GSL_REAL(b);
    number_t cosbxg = ncos(bxg);
    number_t sinbxg = nsin(bxg);

    GSL_REAL(sfvalue(p)) = byg*(mx*cosbxg - my*sinbxg);
    GSL_IMAG(sfvalue(p)) = byg*(my*cosbxg + mx*sinbxg);
    return sfaram1(p);
}

sfarg *sfinveps(sfarg *const p)
{ /* cinv */
    cmplx a = sfarg_or(p, 1, 0, 0);
    /* A hundredth, which softens the pole without moving much else. Written
     * as a division: the literal 0.01 is a double, and promoting it gives a
     * different number from the one this build reads out of "0.01". */
    const number_t hundredth = (number_t)1 / 100;
    cmplx eps = sfarg_or(p, 2, hundredth, hundredth);
    number_t x = GSL_REAL(a);
    number_t y = GSL_IMAG(a);
    number_t delta = (x*x + y*y);
    GSL_REAL(sfvalue(p)) = x/(delta + GSL_REAL(eps));
    GSL_IMAG(sfvalue(p)) = -y/(delta + GSL_IMAG(eps));
    return sfaram1(p);
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

    /* three sides about the origin, corners left where they are */
    gsl_complex centre = sfarg_or(p, 2, 0, 0);
    gsl_complex n = sfarg_or(p, 3, 3, 0);
    gsl_complex zc = gsl_complex_sub(sfarg_or(p, 1, 0, 0), centre);
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
    rn = gsl_complex_mul_real(gsl_complex_pow(rn, sfarg_or(p, 4, 1, 0)),
                              gsl_complex_abs(zc));
    gsl_complex argexp = gsl_complex_exp(gsl_complex_mul_real(i, t));
    sfvalue(p) = gsl_complex_add(gsl_complex_mul(rn, argexp), centre);

    return sfaram1(p);
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
/* The value the pass starts from -- what the formula calls z. The engine
 * writes it before every evaluation and reads back what the formula made
 * of it; the figures below follow it from one pass to the next. */
thread_local cmplx sffe_z = {{0, 0}};

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

/* A cell's value, turned by where in the cell the point stands.
 *
 * These fields hand back one number for a whole cell, which is what makes them
 * hard to colour: the engine colours with the two components of the orbit, so
 * every mode draws one tone a cell -- zmag and iter+real come out looking like
 * the value truncated, and imag and angle have nothing at all, the value never
 * leaving the real axis.
 *
 * Varying the value inside a cell cannot be done with one real number: the
 * same number decides whether the point leaves, so a value that moves across a
 * cell takes half the cell out and leaves the other half in, and the cell comes
 * out cut in two. That was tried twice and cut the mosaics both times.
 *
 * The skew does it with the second component, and does it by turning the value
 * rather than by scaling it. How far it turns follows how far out of the middle
 * of its cell the point stands -- t = skew_re*out + skew_im, where out runs from
 * nought at the middle to one at the edge -- and the value is turned by twice
 * the arc tangent of t.
 *
 * Out is measured in the cell's own geometry, so the colour follows the shape
 * the field is cut into rather than cutting across it: squares get square
 * contours, hexagons hexagonal ones, triangles triangular, a Voronoi cell its
 * own polygon shrunk, and the smooth field follows its own blobs. A straight
 * ramp was tried first and drew diagonal, vertical or horizontal bands across
 * every cell alike -- the same lines in the same direction whatever the field
 * was cut into, which is not a picture of anything.
 *
 * The real part of the skew is that gradient and the imaginary part is a turn
 * the whole cell shares, which shifts its colour without drawing anything
 * across it.
 *
 * A turn and not a scaling, which is what this was first written as. Scaling
 * moved the modulus of the value, and a round bailout looks at exactly that: a
 * cell whose level sat near where the bailout falls came out cut in two, and
 * the mosaics lost their shapes. A turn leaves the modulus where it is, so
 * under a circular bailout the skew is free at any strength -- measured over
 * 90000 pixels, at 0.05, 0.3 and 1, over all five fields: not one pixel leaves
 * differently.
 *
 * A bailout polygon reads the components instead of the modulus, and a turn
 * walks the value round a circle that can cross a side, by as much as the
 * shape's corners stand out past them. So there it does move which cells leave:
 * at a skew of 0.3, seven tenths to two and a half per cent of the picture
 * under a hexagon or a square, two and a half to six and a half under a
 * triangle. That cannot be arranged away -- the escape reads the two components
 * and so does the colouring, so anything the colour can read the escape can see
 * as well. A circular bailout is the way to have it for nothing.
 *
 * (1 + it)^2 / (1 + t^2) is the cosine and the sine of that turn without any
 * trigonometry: its modulus is one exactly, and one to the last bit the single
 * division leaves. A skew of nought gives t = 0 and a factor of exactly one,
 * with no arithmetic done at all, so every value is the number it always was.
 *
 * How much colour it gives is what the skew is worth setting by. The engine's
 * index is (iter + quantity) * speed + shift, so a spread of one in the
 * quantity is a band at a speed of one, and inside a cell imag spreads by 0.14
 * at a skew of 0.05, 0.53 at 0.2, and about one at 0.4 to 0.6. Counting
 * distinct values would call 0.01 enough, at nine hundred of them, but a
 * thousand values inside a hundredth of a band are all one colour.
 *
 * What is left flat is zmag, which reads the modulus and so is the one thing a
 * turn cannot touch. Colour with the modes that read the two components
 * apart -- real, imag, angle, real over imag. */
/* Where in its wedge the point stands: nought on either edge of the wedge, one
 * on the line down its middle, and the same in every wedge -- which is what
 * lets the skew vary with it and leave the picture a kaleidoscope.
 *
 * The fold works the angle out anyway, so this costs one division. Written
 * only by randsc_kaleido, so a field that is not folded pays nothing for it.
 * What is left behind from an earlier call is never read: the skew asks for
 * it only when the call names two wedges or more, and then this has just
 * run. */
static thread_local number_t randsc_inwedge = 0;

/* What the skew is asked to do, which is a set rather than a choice: the bits
 * are added together, so a call may have the shape and the rosette at once and
 * get the spiral that is the two of them.
 *
 * A measure is always needed. A call that names neither of the two that give
 * one gets the shape, which is what a call that says nothing has always got.
 */
#define RANDSC_SKEW_SHAPE 1   /* the cell's own outline, shrunk step by step */
#define RANDSC_SKEW_ROSETTE 2 /* the angle round the middle, folded as the
                                 kaleidoscope folds the plane */
#define RANDSC_SKEW_WEDGE 4   /* a turn that differs from wedge to wedge */
#define RANDSC_SKEW_RADIAL 8  /* the modulus as well as the angle */

/* The angle round the middle of a cell, from nought to one, folded into as
 * many turns as the kaleidoscope has wedges -- so each cell carries a rosette
 * with the picture's own symmetry, and a call that folds nothing gets plain
 * spokes.
 *
 * Taken in the cell's own coordinates rather than in the plane. For a square
 * the two are the same; for the hexagon and the triangle they are sheared, so
 * the rosette leans the way the cell leans. That is the principle the shape
 * measure already follows -- a cell is measured in its own geometry -- and it
 * costs one arc tangent instead of a change of basis.
 */
static number_t randsc_rosette(number_t dx, number_t dy, int wedges)
{
    number_t a = natan2(dy, dx) * (number_t)N_1_2PI + (number_t)1 / 2;
    if (wedges > 1)
        a *= (number_t)wedges;
    return a - nfloor(a);
}

/* The turn itself, given the level rather than the hash it comes from: randsc
 * interpolates its level between four corners and has no single hash to hand
 * over. */
static inline void randsc_skew_apply(number_t level, number_t out, number_t ang,
                                     cmplx skew, int mode, int wedges,
                                     number_t *re, number_t *im)
{
    number_t sr = GSL_REAL(skew), si = GSL_IMAG(skew);
    if (sr == 0 && si == 0) {
        /* the way out, and no arithmetic on the way: the same number as ever */
        *re = level;
        *im = 0;
        return;
    }

    number_t m = 0;
    if ((mode & RANDSC_SKEW_SHAPE) ||
        !(mode & (RANDSC_SKEW_SHAPE | RANDSC_SKEW_ROSETTE)))
        m += out;
    if (mode & RANDSC_SKEW_ROSETTE)
        m += ang;

    /* The imaginary part is a turn the whole cell shares. */
    number_t t = sr * m + si;
    number_t q = 1 + t * t;
    number_t vr = level * (1 - t * t) / q;
    number_t vi = level * 2 * t / q;

    /* Asked to vary with the wedge, the value is turned once more, by nothing
     * on the edges of the wedge and by the whole of the turn the imaginary
     * part gives on its own -- twice its arc tangent, 127 degrees for 2 -- on
     * the line down its middle, and by the turn of the imaginary part times
     * how far across between. Every wedge the same, so a kaleidoscope of six
     * is still one: each copy the turn of the one beside it and the mirror of
     * the other, and continuous across both the edges and the middles, which
     * are where the fold mirrors.
     *
     * It was first a turn that differed from one wedge to the next, the
     * imaginary part times the wedge's number plus one. That made the copies
     * different from one another, which is the one thing a kaleidoscope's
     * copies cannot be -- a six-fold picture came out as six unrelated slices
     * with a straight cut at every join. And the multiple ran up against half
     * a turn, so four slices of six sat nearly on the real axis, where the
     * imaginary part of z hovers about nought and real/imag, which divides by
     * it, drew them as snow. Turning continuously round the whole circle
     * instead mended the cuts but still left one mirror of the six.
     *
     * The Cayley pair again, so no trigonometry: no turn at nought, the whole
     * turn at one.
     *
     * Only where there are wedges to tell apart: randsc_kaleido writes where
     * the point stood and nothing clears it afterwards, so a call that folds
     * nothing would otherwise read whatever the last call that did fold left
     * behind -- and answer differently depending on the order the formula
     * happened to be evaluated in. */
    if ((mode & RANDSC_SKEW_WEDGE) && wedges > 1) {
        number_t tp = si * randsc_inwedge;
        number_t inv = 1 / (1 + tp * tp);
        number_t tr = (1 - tp * tp) * inv, ti = 2 * tp * inv;
        number_t x = vr * tr - vi * ti;
        vi = vr * ti + vi * tr;
        vr = x;
    }

    if (mode & RANDSC_SKEW_RADIAL) {
        /* The modulus as well as the angle, and the only one of the four that
         * zmag and the bailout can see -- a turn leaves the modulus where it
         * is, and those two read nothing else. What it costs is that the
         * escape moves with it, so the figure is cut wherever the bailout
         * falls inside a cell. That is the trade, and it is the reason this is
         * a bit one asks for rather than what the skew does.
         *
         * The absolute value keeps the factor from turning the value inside
         * out where t goes past minus one, which would put a seam along that
         * contour. */
        number_t f = 1 + t;
        if (f < 0)
            f = -f;
        vr *= f;
        vi *= f;
    }

    *re = vr;
    *im = vi;
}

static inline void randsc_skewed(uint64_t k, number_t out, number_t ang,
                                 cmplx skew, int mode, int wedges,
                                 number_t *re, number_t *im)
{
    randsc_skew_apply(randsc_unit(k), out, ang, skew, mode, wedges, re, im);
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
/* The point folded into one wedge of a kaleidoscope.
 *
 * The plane is cut into level equal wedges around the origin and every point
 * is brought back into the first of them, so the field is sampled from one
 * wedge and the picture repeats around the origin. What the mode chooses is
 * which mirror does the folding:
 *
 *   0 (and anything else)  the far half of each wedge is a mirror of the near
 *                          half, so every wedge is symmetric about its own
 *                          bisector;
 *   1                      the same the other way about, the near half
 *                          mirroring the far one.
 *
 * Both are continuous across the joins, so the noise stays coherent and the
 * two precisions go on agreeing; a fold that met itself unevenly would show as
 * a seam.
 *
 * This is the only part of the family that costs trigonometry, and it is only
 * reached when a level of two or more is asked for. At a level of one --
 * which is what a call that says nothing gets -- there is a comparison and
 * nothing else.
 */
static void randsc_kaleido(number_t *px, number_t *py, int level, int mode)
{
    number_t x = *px, y = *py;
    number_t radius = nsqrt(x * x + y * y);
    number_t angle = natan2(y, x);
    number_t sector = 2 * N_PI / (number_t)level;
    number_t turns = nfloor(angle / sector);
    number_t s = angle - sector * turns;

    number_t half = sector / 2;

    /* the nearer edge of the wedge, measured in half wedges, before either
     * mirror decides which half the point is brought into */
    randsc_inwedge = (s < half ? s : sector - s) / half;

    if (mode == 1) {
        if (s < half)
            s = sector - s;
    } else {
        if (s > half)
            s = sector - s;
    }

    *px = radius * ncos(s);
    *py = radius * nsin(s);
}

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

/* pass is the iteration to answer for: the one the formula is on, or an
 * earlier one when the self-similar average is catching up (see randsc_sum),
 * and here is the position. Nothing in here reads sffe_iteration or
 * sffe_position itself; randsc_run reads each once, which matters because a
 * thread-local is a call apiece under MinGW -- some twenty nanoseconds, a
 * tenth of what the whole call costs -- and the sum would otherwise read them
 * again at every pass it works out.
 *
 * lead is how many arguments come before the seed: none for the five, one for
 * randsctile, whose first says which tiling. Everything after the seed means
 * the same in all six, one place further along there. */
static RANDSC_INLINE int randsc_setup(sfarg *const p, unsigned int lead,
                                      unsigned int pass, const cmplx *here,
                                      int64_t *cx, int64_t *cy, number_t *u,
                                      number_t *v, uint64_t *hash)
{
    if (p->argc < 1 + lead || p->argc > 8 + lead)
        return RANDSC_STOP;

    /* Seed, cell size, degradation, kaleidoscope level and its mode, in the
     * order they are written; all but the seed have a default, which a call
     * takes by stopping short of them or by leaving their place empty.
     *
     * Degradation halves the cells each pass. One, which was the default,
     * leaves them the size they started and so wastes the argument on a call
     * that says nothing; a half is the shrinking one asks for when one asks. */
    cmplx seed = sfarg_or(p, lead + 1, 0, 0);
    cmplx size = sfarg_or(p, lead + 2, 1, 1);
    cmplx degradation =
        sfarg_or(p, lead + 3, (number_t)1 / 2, (number_t)1 / 2);
    int level = (int)GSL_REAL(sfarg_or(p, lead + 4, 1, 0));
    int mode = (int)GSL_REAL(sfarg_or(p, lead + 5, 0, 0));

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
    if (p->carried == 0 || p->carried > pass) {
        GSL_SET_COMPLEX(&p->carry, 1, 1);
        p->carried = 0;
    }
    if (GSL_REAL(degradation) == 1 && GSL_IMAG(degradation) == 1) {
        /* The default, and multiplying by one leaves the product where it is,
         * so the passes can be counted off without doing any of them. */
        p->carried = pass;
    } else {
        while (p->carried < pass) {
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
    *hash = randsc_hash((int64_t)pass, 0, randsc_seed(seed));
    number_t px = GSL_REAL(*here), py = GSL_IMAG(*here);
    if (level >= 2)
        randsc_kaleido(&px, &py, level, mode);

    if (!randsc_cell(px / wr, cx, u) || !randsc_cell(py / wi, cy, v)) {
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

/* selfsim, the eighth argument: every pass so far instead of this one alone.
 *
 * Each pass of these functions is a field of its own -- the iteration goes
 * into the hash -- with cells the degradation times the size of the pass
 * before. A formula calling one on every pass therefore sees a new and finer
 * field each time, and what a colouring mode reads off the last pass is that
 * last field alone: once its cells are smaller than a pixel, snow.
 *
 * Asked for, the call hands back instead the weighted average of the passes
 * so far, standing to one another as the octaves of a fractional Brownian
 * motion do. The weight of pass n is d^(nH): d is the degradation, taken as
 * the geometric mean of the absolute values of its two components so that it
 * is one number, and H is the argument. One is the plain motion and a half is
 * rougher; either way the fine passes weigh little, so the snow they make lies
 * under the coarse ones rather than over them. The weights are divided by
 * their total, so the answer stays among the values it averages -- in [0, 1]
 * without a skew. With one it is the skewed value that is averaged, so the
 * skew and its modes go on meaning what they meant.
 *
 * A degradation of one gives every pass the same weight and the answer is
 * their plain average, which settles toward a flat middle as the passes add
 * up. That is what the definition says, and it is done as said.
 *
 * H may be complex, and d^(nH) is then complex too: d^(n Hr) as before for
 * its size, and a turn of n Hi ln d. So every pass is turned that much further
 * than the one before and then averaged as a real H averages it, the weights
 * divided by the sum of their sizes -- which keeps the answer within the
 * largest of the values averaged, however the turns fall. Divided by their
 * complex sum instead, the answer would run off wherever that sum came near
 * nought, and it would come near nought at the same pass for every pixel,
 * the weights being the pixel's own business not at all. The turns make the
 * value complex even without a skew, so imag and angle have something to read.
 * A degradation of one has a logarithm of nought and turns nothing.
 *
 * Left out, or its place left empty, it is off and costs one comparison: the
 * answer is then the one pass, to the bit. Nought is not off. d^(n 0) is one
 * for every pass, so nought is the plain average -- and so is every H near it,
 * from either side and complex too. Nought was the switch at first, which
 * made it the one value the curve does not pass through: 0.000001 was the
 * plain average and 0 the last pass alone, a quarter of the range apart. The
 * last pass alone is the far end of the curve instead, H running to minus
 * infinity. And a real H is averaged exactly as it was before H could be
 * complex: no turn is worked out or applied when there is none to apply.
 *
 * The average is kept on the call site, as the degradation's product is, and
 * kept as an average rather than as a sum and a total. Each pass takes a share
 * of it, and the share of the next follows from the share of this one --
 * a' = s a / (1 + s a), s being d^Hr -- so nothing grows with the passes. A
 * weight that would overflow, which a degradation above one or a negative H
 * ask for, has nowhere to do it. d^Hr and the turn a pass takes are worked out
 * once a pixel; the turn of each pass is the one before times that.
 *
 * The passes a call did not see are worked out when it is asked. A call on a
 * branch taken only from the seventh pass on still gets the first seven, at
 * the price of doing them then; one asked for a pass it has already gone past
 * starts again from nought. Either way the answer at a pass is the same
 * whatever route reached it, a pass depending on nothing but the point and
 * the pass. (And on the arguments, which are taken as they stand when asked;
 * for the constants they nearly always are that is no difference at all.)
 *
 * A new pixel is known by its position as well as by its pass. The pass alone,
 * which is what the trap goes by, would let a call that first runs on the
 * seventh pass of one pixel carry on the average the pixel before it left.
 */
typedef sfarg *(*randsc_pass_fn)(sfarg *const, unsigned int, const cmplx *);

static sfarg *randsc_sum(sfarg *const p, unsigned int lead, cmplx h,
                         unsigned int now, const cmplx *here, randsc_pass_fn at)
{
    unsigned int from = p->summed;

    if (from == 0 || from > now ||
        GSL_REAL(p->gathered) != GSL_REAL(*here) ||
        GSL_IMAG(p->gathered) != GSL_IMAG(*here)) {
        /* d^Hr as |dr di|^(Hr/2), which is the same with one root fewer */
        cmplx degradation =
            sfarg_or(p, lead + 3, (number_t)1 / 2, (number_t)1 / 2);
        number_t dd = nfabs(GSL_REAL(degradation) * GSL_IMAG(degradation));
        GSL_SET_COMPLEX(&p->share, 1, npow(dd, GSL_REAL(h) / 2));
        /* and Hi ln d, the turn from one pass to the next. None at all for a
         * real H, nor for a degradation of one; and none for a degradation of
         * nought, whose passes are all nought and whose logarithm is not a
         * number. */
        number_t psi = 0;
        if (GSL_IMAG(h) != 0 && dd > 0)
            psi = GSL_IMAG(h) * nlog(dd) / 2;
        if (psi != 0 && psi - psi == 0) /* neither NaN nor infinite */
            GSL_SET_COMPLEX(&p->turn, ncos(psi), nsin(psi));
        else
            GSL_SET_COMPLEX(&p->turn, 1, 0);
        GSL_SET_COMPLEX(&p->phase, 1, 0);
        from = 0;
    }

    number_t mr = GSL_REAL(p->mean), mi = GSL_IMAG(p->mean);
    number_t a = GSL_REAL(p->share), s = GSL_IMAG(p->share);
    number_t cr = GSL_REAL(p->phase), ci = GSL_IMAG(p->phase);
    number_t ur = GSL_REAL(p->turn), ui = GSL_IMAG(p->turn);
    int turning = ui != 0 || ur != 1;
    sfarg *last = sfaram1(p);
    for (unsigned int k = from; k <= now; k++) {
        last = at(p, k, here);
        number_t vr = GSL_REAL(sfvalue(p)), vi = GSL_IMAG(sfvalue(p));
        if (k == 0) {
            /* the first pass is all there is so far, and is not turned */
            mr = vr;
            mi = vi;
            continue;
        }
        if (turning) {
            /* this pass's turn is the last one's and one step more */
            number_t x = cr * ur - ci * ui;
            ci = cr * ui + ci * ur;
            cr = x;
            x = vr * cr - vi * ci;
            vi = vr * ci + vi * cr;
            vr = x;
        }
        /* written two ways so that neither end divides nought by nought: a
         * share that has run down to nothing stays nothing, and one whose
         * weight has overflowed takes the whole */
        number_t x = s * a;
        a = x < 1 ? x / (1 + x) : 1 / (1 + 1 / x);
        mr += a * (vr - mr);
        mi += a * (vi - mi);
    }

    GSL_SET_COMPLEX(&p->mean, mr, mi);
    GSL_SET_COMPLEX(&p->share, a, s);
    GSL_SET_COMPLEX(&p->phase, cr, ci);
    p->gathered = *here;
    p->summed = now + 1;
    GSL_SET_COMPLEX(&sfvalue(p), mr, mi);
    return last;
}

/* Each of the six as the formula calls it: the pass it is on, alone, or with
 * a selfsim every pass up to it. The one pass is a direct call. Whether there
 * is a selfsim is whether its place holds anything, not what it holds: nought
 * is a value like any other (see randsc_sum). lead as for randsc_setup. */
static RANDSC_INLINE sfarg *randsc_run(sfarg *const p, unsigned int lead,
                                       randsc_pass_fn at)
{
    unsigned int now = sffe_iteration;
    cmplx here = sffe_position;
    if (p->argc < 8 + lead || p->args[p->argc - 8 - lead]->omitted)
        return at(p, now, &here);
    return randsc_sum(p, lead, sfvalue(p->args[p->argc - 8 - lead]), now,
                      &here, at);
}

/**
 * @brief Coherent noise over the position, seeded and reproducible.
 * @details randsc(seed), randsc(seed, size), randsc(seed, size, degradation).
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
static sfarg *randsc_at(sfarg *const p, unsigned int pass,
                        const cmplx *here)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, 0, pass, here, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        /* past the grid: there is no cell to stand in, so no skew either */
        GSL_SET_COMPLEX(&sfvalue(p), randsc_unit(randsc_hash(cx, cy, h)), 0);
        return sfaram1(p);
    }

    /* how far the value is turned by where in the cell the point stands; nought
     * leaves it where it has always been. See randsc_skewed. */
    cmplx skew = sfarg_or(p, 6, 0, 0);
    /* which of the four things the skew does; see randsc_skew_apply */
    int skewmode = (int)GSL_REAL(sfarg_or(p, 7, 0, 0));
    /* the rosette takes the picture's own symmetry, so it wants the count */
    int wedges = (skewmode & (RANDSC_SKEW_ROSETTE | RANDSC_SKEW_WEDGE))
                     ? (int)GSL_REAL(sfarg_or(p, 4, 1, 0))
                     : 1;

    number_t su = u * u * (3 - 2 * u); /* smoothstep: flat at both ends, so */
    number_t sv = v * v * (3 - 2 * v); /* the value meets its neighbour flat */

    number_t a = randsc_unit(randsc_hash(cx, cy, h));
    number_t b = randsc_unit(randsc_hash(cx + 1, cy, h));
    number_t c = randsc_unit(randsc_hash(cx, cy + 1, h));
    number_t d = randsc_unit(randsc_hash(cx + 1, cy + 1, h));
    number_t lo = a + (b - a) * su;
    number_t hi = c + (d - c) * su;
    number_t level = lo + (hi - lo) * sv;

    if (GSL_REAL(skew) == 0 && GSL_IMAG(skew) == 0) {
        GSL_SET_COMPLEX(&sfvalue(p), level, 0);
        return sfaram1(p);
    }

    /* This field has no cell edges to follow, and what it does have is its own
     * blobs: taking the level itself as the measure puts the turn's contours on
     * the field's contours, so the colour follows the shape that is there. It
     * is continuous across a lattice line for the same reason the level is.
     *
     * The rosette has no middle of its own here either, so it turns about the
     * middle of the lattice square, which is the only thing with a middle. */
    number_t ang = (skewmode & RANDSC_SKEW_ROSETTE)
                       ? randsc_rosette(u - (number_t)1 / 2,
                                        v - (number_t)1 / 2, wedges)
                       : 0;
    number_t sre, sim;
    randsc_skew_apply(level, level, ang, skew, skewmode, wedges, &sre, &sim);
    GSL_SET_COMPLEX(&sfvalue(p), sre, sim);
    return sfaram1(p);
}

sfarg *sfrandsc(sfarg *const p)
{
    return randsc_run(p, 0, randsc_at);
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

/* The seed nearest the point among the nine cells round the one it falls in
 * -- the Voronoi cell it stands in, named by that seed's hash -- and, when
 * asked, how far out of that cell's middle it stands and where the seed is
 * seen from it. Shared by randscp and the irregular tiling of randsctile. */
static inline uint64_t randscp_nearest(int64_t cx, int64_t cy, number_t u,
                                       number_t v, uint64_t h, int want,
                                       number_t *pout, number_t *pbx,
                                       number_t *pby)
{
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

    /* How far the point stands from the edge of its cell, which is what puts
     * the contours on the polygon.
     *
     * A cell here is a convex polygon, and the points a given distance inside
     * its edge are that polygon shrunk, so a value that follows the distance
     * to the edge draws the cell's own outline over and over -- polygons, and
     * each cell its own. The distance to the boundary the winning seed shares
     * with another is (|b|^2 - |a|^2) / (2|b - a|), a and b being the two seeds
     * seen from the point, and the edge is the nearest of the eight.
     *
     * The comparison cross-multiplies the squares, so the candidates cost no
     * division and no root and one of each is taken at the end. A distance is
     * never negative here, the winner being the nearest seed there is, so
     * squaring loses nothing. None of it is done at all when there is no skew
     * to feed.
     *
     * Nought at the edge on both sides of it, so the turn meets itself across
     * a boundary as the other four do. */
    *pout = 0;
    *pbx = 0;
    *pby = 0;
    if (want) {
        /* The nine seeds again, kept this time, so that the search above needs
         * to carry nothing for a skew that is usually not there: it runs for
         * every point, and this runs only when there is one. The winner is the
         * one whose hash is the winning hash. */
        number_t dxs[9], dys[9];
        int bidx = 0;
        for (int j = -1, k = 0; j <= 1; j++)
            for (int i = -1; i <= 1; i++, k++) {
                uint64_t hh = randsc_hash(cx + i, cy + j, h);
                number_t fx =
                    RANDSCP_JITTER_LOW +
                    RANDSCP_JITTER_SPAN * randsc_unit32((uint32_t)hh);
                number_t fy =
                    RANDSCP_JITTER_LOW +
                    RANDSCP_JITTER_SPAN * randsc_unit32((uint32_t)(hh >> 32));
                dxs[k] = (number_t)i + fx - u;
                dys[k] = (number_t)j + fy - v;
                if (hh == besth)
                    bidx = k;
            }
        number_t bx = dxs[bidx], by = dys[bidx];
        *pbx = bx;
        *pby = by;
        /* out of reach of any real candidate, so the first one always wins */
        number_t bestn = 256, bests = 1;
        /* and the closest other seed, which is what makes the cell its own
         * measure: the seed's distance to the boundary it shares with another
         * is exactly half their separation, so the smallest of those halves is
         * the distance from the seed to its cell's edge */
        number_t nearsep = 256;
        for (int k = 0; k < 9; k++) {
            if (k == bidx)
                continue;
            number_t dx = dxs[k], dy = dys[k];
            number_t ex = dx - bx, ey = dy - by;
            number_t sep = ex * ex + ey * ey;
            if (sep <= 0)
                continue;
            if (sep < nearsep)
                nearsep = sep;
            number_t num = dx * dx + dy * dy - bestd;
            number_t n2 = num * num;
            if (n2 * bests < bestn * sep) {
                bestn = n2;
                bests = sep;
            }
        }
        /* the distance to the edge over the seed's own distance to it, which
         * is nought at the seed and one along the whole boundary however
         * lopsided the cell. Both halves cancel, leaving one root and one
         * division for the lot. */
        number_t pin = 1 - nsqrt(bestn / (bests * nearsep));
        *pout = pin < 0 ? 0 : (pin > 1 ? 1 : pin);
    }

    return besth;
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
static sfarg *randscp_at(sfarg *const p, unsigned int pass,
                         const cmplx *here)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, 0, pass, here, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p),
                        randsc_unit(randsc_remix(randsc_hash(cx, cy, h))),
                        0);
        return sfaram1(p);
    }

    cmplx skew = sfarg_or(p, 6, 0, 0);
    /* which of the four things the skew does; see randsc_skew_apply */
    int skewmode = (int)GSL_REAL(sfarg_or(p, 7, 0, 0));
    /* the rosette takes the picture's own symmetry, so it wants the count */
    int wedges = (skewmode & (RANDSC_SKEW_ROSETTE | RANDSC_SKEW_WEDGE))
                     ? (int)GSL_REAL(sfarg_or(p, 4, 1, 0))
                     : 1;

    number_t pout, pbx, pby;
    uint64_t besth = randscp_nearest(cx, cy, u, v, h,
                                     GSL_REAL(skew) != 0 || GSL_IMAG(skew) != 0,
                                     &pout, &pbx, &pby);

    number_t pre_, pim_;
    /* round the seed the cell was grown from, which is its middle */
    number_t pang = (skewmode & RANDSC_SKEW_ROSETTE)
                        ? randsc_rosette(-pbx, -pby, wedges)
                        : 0;
    randsc_skewed(randsc_remix(besth), pout, pang, skew, skewmode, wedges, &pre_,
                  &pim_);
    GSL_SET_COMPLEX(&sfvalue(p), pre_, pim_);
    return sfaram1(p);
}

sfarg *sfrandscp(sfarg *const p)
{
    return randsc_run(p, 0, randscp_at);
}


/**
 * @brief A fractional Brownian motion over a point the caller names.
 * @details fbm(value, seed), and up to
 * fbm(value, seed, intensity, frequency, octaves, roughness).
 *
 * The same noise the randsc family is built from, summed in octaves: each at
 * twice the frequency of the one before and keeping a share of its height, so
 * that no octave is large enough to see on its own and none is small enough to
 * disappear. What that draws is wear rather than a pattern -- stains, dents,
 * the surface of something that has been left out.
 *
 * Where randsc reads the position and can only read the position, this reads
 * whatever is written in front of it: fbm(z, 7) moves with the orbit,
 * fbm(x, 7) stands still on the plane, and fbm(z*3+c, 7) is whatever that is.
 * That is the whole reason for it: the colouring modes of the same name apply
 * a motion to a colour, and this puts one where a formula can use it.
 *
 * The two mandatory arguments are the point and the seed. The seed is read as
 * randsc reads one, so the same number means the same field in both.
 *
 * intensity  what the motion is multiplied by. It runs from nought to this and
 *            never below, so that adding it to something cannot pull that
 *            under nought -- the colouring modes learned the same lesson the
 *            long way round.
 * frequency  cells of the lattice to a unit of the plane: the size of the
 *            marks, and what to raise as one zooms in.
 * octaves    how many are summed. Held between one and twenty-four.
 * roughness  what each octave keeps of the height of the one before. A half is
 *            the plain motion; higher is grittier, and only then do the later
 *            octaves carry enough to be worth asking for.
 *
 * No kaleidoscope of its own, where the randsc family has one: parchment and
 * parchmenta already fold the plane into sectors, so fbm(parchmenta(z, 6), 7)
 * says it, and says it where anyone reading the formula can see it.
 *
 * Each octave is given a seed of its own. Sharing one would leave every octave
 * agreeing wherever the lattices agree -- the origin, and every point that
 * lands on a corner -- which shows as a knot in the field.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
sfarg *sffbm(sfarg *const p)
{
    if (p->argc < 2 || p->argc > 6) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }

    cmplx at = sfarg_or(p, 1, 0, 0);
    uint64_t seed = randsc_seed(sfarg_or(p, 2, 0, 0));
    number_t much = GSL_REAL(sfarg_or(p, 3, 4, 0));
    number_t freq = GSL_REAL(sfarg_or(p, 4, 8, 0));
    int octaves = (int)GSL_REAL(sfarg_or(p, 5, 4, 0));
    number_t rough = GSL_REAL(sfarg_or(p, 6, (number_t)1 / 2, 0));

    if (octaves < 1)
        octaves = 1;
    if (octaves > 24)
        octaves = 24;
    if (!(rough > 0))
        rough = (number_t)1 / 2;

    number_t x = GSL_REAL(at) * freq, y = GSL_IMAG(at) * freq;
    number_t sum = 0, amp = 1, norm = 0;

    for (int i = 0; i < octaves; i++) {
        int64_t cx, cy;
        number_t u, v;
        /* Past the resolution of the lattice there is no cell to stand in, and
         * every octave after this one is finer still: what has been summed so
         * far is all there is to have. */
        if (!randsc_cell(x, &cx, &u) || !randsc_cell(y, &cy, &v))
            break;

        number_t su = u * u * (3 - 2 * u); /* smoothstep: flat at both ends, */
        number_t sv = v * v * (3 - 2 * v); /* so a cell meets its neighbour  */

        uint64_t h = randsc_hash((int64_t)i, 0, seed);
        number_t a = randsc_unit(randsc_hash(cx, cy, h));
        number_t b = randsc_unit(randsc_hash(cx + 1, cy, h));
        number_t c = randsc_unit(randsc_hash(cx, cy + 1, h));
        number_t d = randsc_unit(randsc_hash(cx + 1, cy + 1, h));
        number_t lo = a + (b - a) * su;
        number_t hi = c + (d - c) * su;

        sum += amp * (lo + (hi - lo) * sv);
        norm += amp;
        amp *= rough;
        x *= 2;
        y *= 2;
    }

    GSL_SET_COMPLEX(&sfvalue(p), norm > 0 ? sum / norm * much : 0, 0);
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
/**
 * @brief A polynomial in the first argument, the rest being its coefficients.
 * @details poly(z, k1, k2, ..., km) is
 *
 *     k1*z^(m-1) + k2*z^(m-2) + ... + k(m-1)*z + km
 *
 * so the first coefficient written multiplies the highest power and the last
 * stands alone. Writing the coefficients in the order one says them is what
 * makes the call read like the polynomial.
 *
 * Worked out by Horner's rule -- start at k1 and repeatedly multiply by z and
 * add the next -- which is m-1 multiplications rather than the m(m-1)/2 that
 * raising each power separately would take, and is the more accurate of the
 * two into the bargain, every product being formed from one that was already
 * rounded once rather than from a power rounded m times.
 *
 * With no coefficients at all there are no terms, and a sum of no terms is
 * zero; so is a coefficient whose place was left empty, "poly(z, 1, , 1)"
 * being z^2 + 1.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
sfarg *sfpoly(sfarg *const p)
{
    if (p->argc < 2) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        /* the call itself, there being no argument to point at */
        return p;
    }
    /* args[argc - 1] is the first written, so z; the coefficients follow it
     * from args[argc - 2] down to args[0]. */
    cmplx z = sfvalue(p->args[p->argc - 1]);
    cmplx acc = sfvalue(p->args[p->argc - 2]);
    for (int k = p->argc - 3; k >= 0; k--)
        acc = gsl_complex_add(gsl_complex_mul(acc, z), sfvalue(p->args[k]));
    sfvalue(p) = acc;
    return sfaram1(p);
}

/* Orbit traps and stripe averaging.
 *
 * Both watch the orbit go by and keep one number about the whole of it -- the
 * nearest the orbit ever came to a shape, or the average of a wave taken along
 * it. That is the sort of quantity the colouring modes of the engine cannot
 * reach: a calculation loop there is compiled to hand back a colour, and only
 * the last z and the one before it survive to be coloured by. Done here it
 * costs the engine nothing.
 *
 * The running quantity lives on the call site, as the noise degradation does
 * (see sfarg): two traps in one formula keep their own, a thread cannot
 * disturb another, and a pass earlier than the one counted means a new pixel
 * has begun. Unlike the degradation this really is a property of the pixel --
 * it is that pixel's orbit being watched -- so it depends on the engine
 * starting every pixel at pass zero, which INIT does.
 *
 * Both hand back the value they were given until the last pass the iteration
 * limit allows, where they hand back what they have gathered instead. So a
 * whole formula is trap(z^2+c, 3) and nothing else: the fractal iterates as it
 * would, and on the last pass the value becomes the trap, which the inside
 * colouring modes then draw. A point that escapes never reaches that pass and
 * keeps its ordinary outside colour, so what these draw is the inside.
 */

/* Distance from a point to one of the shapes an orbit can be trapped by.
 * Numbered rather than named because the shape is an argument; the numbering
 * is the order the help lists them in. */
static number_t sftrap_distance(number_t re, number_t im, int shape,
                                number_t radius)
{
    number_t ar = nfabs(re), ai = nfabs(im);
    switch (shape) {
        case 1: /* the horizontal line through the centre */
            return ai;
        case 2: /* the vertical one */
            return ar;
        case 3: /* both of them, a cross */
            return ar < ai ? ar : ai;
        case 4: /* a ring of that radius */
            return nfabs(nsqrt(re * re + im * im) - radius);
        case 5: /* the square with that half-side */
            return nfabs((ar > ai ? ar : ai) - radius);
        case 6: /* the diamond with that half-diagonal */
            return nfabs(ar + ai - radius);
        default: /* the centre itself */
            return nsqrt(re * re + im * im);
    }
}

/* Has this pixel only just begun? Clears what the last one left behind. */
static int sftrap_begin(sfarg *const p)
{
    if (p->carried == 0 || p->carried > sffe_iteration) {
        GSL_SET_COMPLEX(&p->carry, 0, 0);
        return 1;
    }
    return 0;
}

/* Is this the last pass the limit allows? Where nobody has said what the
 * limit is, nothing is ever revealed, which is what a bare parser wants. */
static int sftrap_last(void)
{
    return sffe_maxiter && sffe_iteration + 1 >= sffe_maxiter;
}

/**
 * @brief How near the orbit came to a shape.
 * @details trap(a, shape, centre, size) measures the distance from a to the
 * shape, keeps the smallest seen so far, and hands back a unchanged until the
 * last pass, where it hands back that smallest distance instead. shape
 * defaults to 0, centre to the origin, size to 1.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
sfarg *sftrap(sfarg *const p)
{
    if (p->argc < 1 || p->argc > 4) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return p;
    }

    /* the point watched, then the shape, where it sits and how big it is, all
     * three of which have a default */
    cmplx a = sfarg_or(p, 1, 0, 0);
    int shape = (int)GSL_REAL(sfarg_or(p, 2, 0, 0));
    cmplx centre = sfarg_or(p, 3, 0, 0);
    cmplx size = sfarg_or(p, 4, 1, 0);

    number_t d = sftrap_distance(GSL_REAL(a) - GSL_REAL(centre),
                                 GSL_IMAG(a) - GSL_IMAG(centre), shape,
                                 GSL_REAL(size));
    if (sftrap_begin(p) || d < GSL_REAL(p->carry))
        GSL_SET_COMPLEX(&p->carry, d, 0);
    p->carried = sffe_iteration + 1;

    if (sftrap_last())
        GSL_SET_COMPLEX(&sfvalue(p), GSL_REAL(p->carry), 0);
    else
        sfvalue(p) = a;
    return sfaram1(p);
}

/**
 * @brief A wave averaged along the orbit.
 * @details stripe(a, density) averages (sin(density * arg a) + 1) / 2 over the
 * passes so far and hands back a unchanged until the last one, where it hands
 * back that average. density defaults to 4 and is how many stripes go round a
 * turn; a whole number, or the stripes do not meet where the turn closes.
 *
 * An average over the orbit changes smoothly with the point even where the
 * iteration count jumps, which is what draws the fibres the method is used
 * for.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the last argument, per the sffe convention.
 */
sfarg *sfstripe(sfarg *const p)
{
    if (p->argc < 1 || p->argc > 2) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return p;
    }

    cmplx a = sfarg_or(p, 1, 0, 0);
    number_t density = GSL_REAL(sfarg_or(p, 2, 4, 0));

    number_t sample =
        (nsin(density * natan2(GSL_IMAG(a), GSL_REAL(a))) + 1) / 2;
    sftrap_begin(p);
    /* the sum in the real part, how many went into it in the imaginary one */
    GSL_SET_COMPLEX(&p->carry, GSL_REAL(p->carry) + sample,
                    GSL_IMAG(p->carry) + 1);
    p->carried = sffe_iteration + 1;

    if (sftrap_last())
        GSL_SET_COMPLEX(&sfvalue(p), GSL_REAL(p->carry) / GSL_IMAG(p->carry),
                        0);
    else
        sfvalue(p) = a;
    return sfaram1(p);
}


/* --- figures: a field drawn rather than diced -----------------------------
 *
 * These three ask which part of a figure the point is standing in and carry it
 * one step towards the outside of that figure. Nothing about them is random --
 * that is the whole of what they share with the noise functions -- so the same
 * point always takes the same number of steps, and the figure is where it looks
 * like it is. Mandelbrot mode and julia mode draw the same picture, since a
 * figure is a shape in the plane and neither mode changes where a shape is.
 *
 * All three are fractals of their own and not fields to multiply into one.
 * Written alone -- "sierpinskyt()" and nothing else -- each draws its figure,
 * the way the Sierpinski, Sierpinski Carpet and Koch Snowflake under Fractal ->
 * More Formulae draw theirs.
 *
 * All three carry the point to its parent. Every part of one of these figures
 * has one: a hole in a gasket sits inside a bigger hole one level up, a hole in
 * a carpet inside the cell that was cut the same way, a triangle of a snowflake
 * on a side of the hexagon at its middle or on the free edge of another
 * triangle. The step onto the parent is a motion of the plane -- doubling away
 * from the nearest corner for the gasket, blowing a cell up by the number of
 * cells for the carpet, folding a triangle across the hexagon side it stands on
 * or blowing it up by three and turning it to face the way its edge faces for
 * the snowflake. The topmost part has no parent inside the figure, so its step
 * carries it out of the bailout, and a part n levels down takes n steps to get
 * there.
 *
 * The point really travels: the pass it leaves on is the level it stands at,
 * the iteration count is the picture, and z along the way is a point of the
 * plane like any other -- which is what lets the colouring modes that read z,
 * smooth colouring, and writing a figure inside a larger formula all mean
 * something.
 *
 * That holds for the pass it goes on as well, and it has to. Every outside
 * colouring mode but the iteration count reads what came back on that pass, so
 * a number standing in for "gone" would be one number for the whole figure and
 * all of them would draw one flat tone. The topmost part is therefore thrown
 * rather than declared gone: the gasket doubles its middle hole away from a
 * corner like any other and it lands outside the triangle, which is outside the
 * bailout; the carpet scales its middle cell about the cell beside it and it
 * lands a square out; the snowflake throws its middle hexagon out the way the
 * sector faces. All three carry where they came from with them. Only a point
 * that was never in the figure at all is handed the standing-in number.
 *
 * A gasket and a carpet fill their shape, so every point in one has a level and
 * leaves on it. A snowflake does not fill its hexagon: it leaves six corners of
 * ground, which is no part of the figure and has no level. That ground is
 * handed back where it stands and never leaves, so it is drawn in the inside
 * colour rather than taking a band of the outside one off the figure, and the
 * incolouring modes have the point itself to read.
 *
 * Where the levels of a snowflake are counted from is the difference between a
 * picture and a flat tone. Counted from the triangle it grew from, five eighths
 * of it is level one and everything else is a fringe around that. Counted from
 * the regular hexagon at its middle -- which with the six triangles standing on
 * that hexagon's six sides is the same figure, exactly -- the first two levels
 * are five twelfths each and the rest come down evenly, and what is drawn looks
 * like a snowflake instead of like a triangle.
 *
 * They were a field before this: a number between nought and one saying how
 * solidly a point belonged, meant to be multiplied into a formula the way the
 * noise is. Written alone a field draws nothing at all -- a number between
 * nought and one cannot leave a bailout of four -- and what was on the screen
 * was the bailout shape in one flat tone, which is not what anybody asking for
 * a Sierpinski triangle wants.
 *
 * radius means what bailout means and is read the same way: as the square of
 * the distance. What it draws is then the shape a bailout of that number
 * draws, inscribed in it exactly:
 *
 *   sierpinskyt   the triangle of BAILOUT_TRIANGLEM90, corner for corner
 *   sierpinskyc   the square of BAILOUT_SQUARE, side for side
 *   snowflake     the hexagon of BAILOUT_HEXAGON0, its six points on the
 *                 hexagon's six corners
 *
 * So a figure written with the same number as the bailout stands exactly where
 * the bailout stands, and the two lie over one another.
 *
 * The thing to know is that a bailout polygon stands its *sides* the square
 * root of the bailout from the centre -- that is its apothem, not its
 * circumradius -- so its corners are further out than the number says: twice
 * as far for a triangle, and by a seventh for a hexagon. A figure drawn to the
 * number instead of to the shape comes out half the size the shape is, which
 * is what these did until it was measured off a picture of the two together.
 *
 * None of them iterates over the figure. A pass of the gasket is three
 * multiply-adds and a comparison, a pass of the carpet two divisions worth of
 * scaling and two truncations, and a pass of the snowflake four multiplications
 * to name the sector and, for the six sevenths of it that is the hexagon or a
 * triangle on the hexagon, nothing further -- only what is deeper walks the
 * curve, and that walk is one path down and not a search. So what a pass costs
 * does not depend on how much of a figure the picture can see, nor on how deep
 * the level it is looking at.
 */

/* sin 60 degrees, which is half the width of an equilateral triangle of
 * circumradius one and the height of the bump the Koch curve puts on a third of
 * a segment, three times over. Worked out at the precision in use rather than
 * written as a decimal, which would be a double in the quad build. */
static const number_t FIG_SIN60 = nsqrt((number_t)3) / 2;
/* the corners of that triangle, apex up */
static const number_t FIG_VX[3] = {0, -FIG_SIN60, FIG_SIN60};
static const number_t FIG_VY[3] = {1, (number_t)-1 / 2, (number_t)-1 / 2};
/* one over that, which is what turns an apothem into a circumradius */
static const number_t FIG_ACROSS = 1 / (nsqrt((number_t)3) / 2);

/* Each edge of that triangle, as the pair a point is multiplied by to say how
 * far along it and how far out of it the point stands. Both are the edge
 * divided by its own length twice over, and the edge and its length are the
 * same for every picture ever drawn -- so they are worked out here rather than
 * three times a pass for every pixel. At 113 bits of mantissa a division is
 * software and costs about what a hundred additions cost, which is what made
 * this worth doing rather than merely tidy. */
static const number_t FIG_EDGE = 3; /* every edge, squared, of a unit triangle */
static const number_t FIG_AX[3] = {
    (-FIG_SIN60 - 0) / FIG_EDGE, (FIG_SIN60 + FIG_SIN60) / FIG_EDGE,
    (0 - FIG_SIN60) / FIG_EDGE};
static const number_t FIG_AY[3] = {((number_t)-1 / 2 - 1) / FIG_EDGE, 0,
                                   (1 - (number_t)-1 / 2) / FIG_EDGE};
/* The child a triangle in a snowflake grows on each of its edges: a third the
 * size, standing two thirds of the way out along the facing of the corner that
 * edge is opposite, and turned so that its own point faces the same way. The
 * step onto the parent undoes that turn, so what is kept is its cosine and its
 * sine -- half a turn for the child under the base, a sixth either way for the
 * two above it. */
static const int FIG_OPP[3] = {2, 0, 1}; /* the corner an edge is opposite */
static const number_t FIG_CC[3] = {-1, (number_t)1 / 2, (number_t)1 / 2};
static const number_t FIG_CS[3] = {0, -FIG_SIN60, FIG_SIN60};
/* Of a triangle standing on a side of the hexagon at the middle of a snowflake,
 * the base is glued to the hexagon and only these two edges bear children. */
static const int FIG_FREE[2] = {0, 2};

/* A snowflake read from its middle out. The hexagon there has the same side as
 * the six triangles standing on its six sides, and stands its own sides half
 * way out to the points of the figure -- so in units of the figure's reach, its
 * apothem is a half and its corners are at one over the square root of three.
 *
 * The six sectors are told apart by which of six ways the point stands furthest
 * out on, and the one it is in is then turned upright so that one piece of
 * geometry serves all six. Turning by sixty degrees less sixty times the sector
 * does it; the figure has six-fold symmetry from the hexagon out, which is what
 * makes one piece of geometry enough. */
static const number_t FIG_HEX_APOTHEM = (number_t)1 / 2;
static const number_t FIG_SQRT3 = nsqrt((number_t)3);
static const number_t FIG_THIRD = (number_t)1 / 3;
static const number_t FIG_TURN_C[6] = {(number_t)1 / 2, 1,  (number_t)1 / 2,
                                       (number_t)-1 / 2, -1, (number_t)-1 / 2};
static const number_t FIG_TURN_S[6] = {FIG_SIN60,  0, -FIG_SIN60,
                                       -FIG_SIN60, 0, FIG_SIN60};

/* Three times the figure's reach, the way each sector faces.
 *
 * The hexagon at the middle has no parent to walk to and so has to go, and what
 * it hands back on the way out is what every outside colouring mode but the
 * iteration count has to read. A number standing in for "gone" would be the
 * same number for every pixel of it, and every one of those modes would draw
 * one flat tone -- so it is thrown, and what it hands back is a point of the
 * plane that moves as the point moves. The gasket has always done this: a hole
 * with no parent is doubled away from a corner like any other and lands outside
 * the triangle, which is outside the bailout, carrying where it came from with
 * it.
 *
 * Three times the reach clears every bailout shape of the same number. The
 * worst of them is the triangle, whose corners stand at twice its apothem and
 * so at 1.74 of the reach; the nearest this throw can leave a point is three
 * less the hexagon's own corner, 2.42 of the reach, which is past it. Within a
 * sector the throw is a translation, so the point keeps its shape exactly, and
 * the six sectors are the same throw turned.
 *
 * Written in the plane rather than in the figure's units, and so with the
 * turning of an apothem into a circumradius already in it: what is left to do
 * per pass is to multiply by the square root of the radius, which comes back
 * from multiplying the radius into the reciprocal that is already cached. The
 * point itself is then added in the plane, where it already is, so the whole
 * throw is three multiplications and two additions. */
static const number_t FIG_THROW_X[6] = {3, 0, -3, -3, 0, 3};
static const number_t FIG_THROW_Y[6] = {FIG_SQRT3,  2 * FIG_SQRT3, FIG_SQRT3,
                                        -FIG_SQRT3, -2 * FIG_SQRT3, -FIG_SQRT3};


/* The square root of the radius, or its reciprocal, kept on the call site
 * rather than taken again on every pass.
 *
 * At 113 bits of mantissa a square root and a division are both software and
 * either costs about what the whole of a figure costs; taking one of each per
 * pass put these over the noise they are meant to stay under. The radius is
 * almost always a constant, so remembering the one it was last asked for turns
 * that into a comparison and a load. The scratch belongs to the call site (see
 * sfarg), so two figures in one formula keep their own and no thread can
 * disturb another. */
static inline number_t fig_cached(sfarg *const p, number_t radius, int recip)
{
    if (p->carried && GSL_REAL(p->carry) == radius)
        return GSL_IMAG(p->carry);
    number_t v = nsqrt(radius);
    if (recip)
        v = 1 / v;
    GSL_SET_COMPLEX(&p->carry, radius, v);
    p->carried = 1;
    return v;
}

/* What a figure hands back for a point that is no part of it at all: the point
 * turned by the figure's own symmetry.
 *
 * Turned, it is still outside the figure -- a third of a turn leaves a triangle
 * where it was, a quarter leaves a square, half leaves a snowflake -- and still
 * inside the bailout that matches the figure, so it never leaves. That puts it
 * in the inside colour with its own position for the incolouring to read: an
 * empty space stays empty rather than taking the lowest band of the outside
 * colour, which is what filled the margin between a figure and a bailout larger
 * than its radius with one flat tone.
 *
 * Turned rather than handed back where it stood, for the reason the ground of a
 * snowflake is: a fixed point is a stopped orbit and not a bounded one, and a
 * figure written inside a larger formula would hand that formula back the point
 * it was already holding.
 *
 * A turn and not a mirror, though a mirror would keep it outside every one of
 * these and inside every bailout shape as well: a turn is what carries the
 * value round with the point, so that turning the picture turns what is drawn
 * on it, and a mirror turns it the other way. The distance from the centre is
 * left alone either way.
 *
 * The point is handed in rather than read off z, because a kaleidoscope may
 * have folded it: what is turned has to be the point the figure was read at,
 * or the empty space would be turned out of the wedge the rest of the figure
 * lives in. */
#define FIG_TURN(p, cs, sn, zx, zy)                                            \
    GSL_SET_COMPLEX(&sfvalue(p), (cs) * (zx) - (sn) * (zy),                    \
                    (sn) * (zx) + (cs) * (zy))

/* The point a figure reads, folded into one wedge if the call asks for it.
 *
 * The noise fields fold the position, which stands still for the whole orbit.
 * A figure has no use for the position: it reads z, the point it is carrying
 * down towards the part of the figure it stands in, so z is what is folded
 * here. What comes of it is the same -- the figure is drawn in one wedge and
 * repeated round the origin -- and it holds all the way down, every step being
 * taken on the folded point and folded again on the pass after.
 *
 * The two arguments are the ones the noise family takes and mean the same: how
 * many wedges, and which mirror folds them. They come last, so a call written
 * before they existed is the call it was.
 *
 * A level of one -- which is what a call that says nothing gets -- folds
 * nothing and works out no angle, so the figure is the figure it always was,
 * to the bit. */
static inline void fig_point(sfarg *const p, unsigned int place, number_t *zx,
                             number_t *zy)
{
    *zx = GSL_REAL(sffe_z);
    *zy = GSL_IMAG(sffe_z);
    int level = (int)GSL_REAL(sfarg_or(p, place, 1, 0));
    if (level >= 2)
        randsc_kaleido(zx, zy, level,
                       (int)GSL_REAL(sfarg_or(p, place + 1, 0, 0)));
}
/* Twenty-four levels puts the finest bump at three to the minus twenty-four of
 * the radius, which is smaller than anything a picture of a figure this size
 * will ever be asked to show. What the mantissa could carry is beside the
 * point: the levels past this one cost and are not seen. */
#define KOCH_DEPTH 24
/* sin 60 again, in the type the walk is done in */
#define KOCH_SIN60 0.86602540378443864676
/* What the walk says when it ran out of levels without deciding: the point is
 * not known to be ground, so it is carried on as a point of the figure and
 * asked again next pass, one level further down. */
#define KOCH_UNDECIDED (-1)
/**
 * @brief The Sierpinski gasket, as a field over the plane.
 * @details sierpinskyt(radius, kaleidoscope, mode) stands the triangle a
 * triangular bailout of that number draws, point upwards, and cuts the gasket
 * out of it: a point leaves on the pass numbered by the cut that took it, and
 * one on the gasket itself never leaves, so the iteration count is the
 * picture. The last two are the kaleidoscope the noise family takes, one and
 * nought by default; see fig_point.
 *
 * Worked out in barycentric coordinates, where the gasket has a description
 * that costs nothing. Halving the triangle towards each of its corners in turn
 * is halving two of the three weights, so writing those two in binary spells
 * out which corner was taken at every level -- and the gasket is exactly the
 * points where no place holds a one in both, since a one in both is the step
 * that would have to go towards two corners at once. The first such place is
 * the level the point was cut away at.
 *
 * So the figure is two multiplications, two conversions and an AND, at any
 * depth, where subdividing a triangle thirty times would have been thirty
 * rounds of comparisons for every pixel of every pass. The bands it draws
 * therefore follow the halving rather than the three-fold symmetry of the
 * triangle: the figure is symmetric and the shading is symmetric about the one
 * axis, which is a thing to know rather than a thing to fix.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the call, the evaluator having no use for the result.
 */
sfarg *sfsierpinskyt(sfarg *const p)
{
    GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
    if (p->argc > 3)
        return p;
    number_t radius = GSL_REAL(sfarg_or(p, 1, 4, 0));
    if (!(radius > 0))
        return p;

    /* The sides stand the square root of the radius from the centre, where a
     * triangular bailout of that number stands its sides, which puts the
     * corners at twice that: a triangle's apothem is half its circumradius. */
    number_t corner = 2 * fig_cached(p, radius, 0);
    number_t x, y;
    fig_point(p, 2, &x, &y);

    /* The three barycentric weights, each multiplied by three times the corner
     * distance. That is a positive number, so it changes neither which of them
     * is the largest nor which of them is negative, and it takes the divisions
     * out: with the apex at (0,C) and the base along y = -C/2, the weights come
     * out of the point in three multiply-adds. */
    number_t wa = 2 * y + corner;
    number_t across = FIG_SIN60 * 2 * x;
    number_t wb = corner - y - across;
    number_t wc = corner - y + across;
    if (wa < 0 || wb < 0 || wc < 0) {
        /* out of the triangle, so no part of the figure: a third of a turn,
         * which is the triangle's own, and it never leaves */
        FIG_TURN(p, (number_t)-1 / 2, FIG_SIN60, x, y);
        return p;
    }

    /* Double away from the corner whose share of the point is largest, which
     * is the corner of the sub-triangle it stands in. A point in the middle
     * triangle has no share over a half anywhere, so whichever corner is
     * chosen it lands outside -- which is what makes the topmost hole the one
     * that leaves. */
    int k = (wa >= wb && wa >= wc) ? 0 : (wb >= wc ? 1 : 2);
    GSL_SET_COMPLEX(&sfvalue(p), 2 * x - corner * FIG_VX[k],
                    2 * y - corner * FIG_VY[k]);
    return p;
}

/**
 * @brief The Sierpinski carpet, as a field over the plane.
 * @details sierpinskyc(radius, squares, kaleidoscope, mode) fills the square a
 * square bailout of that number draws -- half a side the square root of radius, at the origin --
 * cuts it into squares by squares, keeps the ring of cells along the border,
 * throws away everything the ring encloses, and does the same to each cell it
 * kept. So the picture is one square in the middle and 4*squares-4 around it:
 * eight at three squares to a side, which is the carpet as it is usually
 * drawn, twelve at four, sixteen at five. Two is the one number with no ring
 * to speak of, and there the far corner goes instead, which is a gasket again
 * -- a square cut in four with one corner taken away is what a gasket is.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the call, the evaluator having no use for the result.
 */
sfarg *sfsierpinskyc(sfarg *const p)
{
    GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
    if (p->argc > 4)
        return p;
    number_t radius = GSL_REAL(sfarg_or(p, 1, 4, 0));
    int squares = (int)GSL_REAL(sfarg_or(p, 2, 3, 0));
    if (!(radius > 0) || squares < 2 || squares > 64)
        return p;

    /* How far along the square the point stands, counted in cells rather than
     * in the plane. The reciprocal of the root is what is kept on the call
     * site, a root costing more than anything done with it here; the root
     * itself comes back by multiplying the square into that reciprocal --
     * radius over its own root -- rather than by dividing, a division at 113
     * bits of mantissa being software and costing about what a hundred
     * additions cost. */
    number_t invhalf = fig_cached(p, radius, 1);
    number_t half = radius * invhalf;
    number_t zx, zy;
    fig_point(p, 3, &zx, &zy);
    number_t tu = (zx + half) * invhalf * squares / 2;
    number_t tv = (zy + half) * invhalf * squares / 2;
    if (!(tu >= 0) || tu >= squares || !(tv >= 0) || tv >= squares) {
        /* out of the square, so no part of the figure: a quarter of a turn,
         * which is the square's own, and it never leaves */
        FIG_TURN(p, 0, 1, zx, zy);
        return p;
    }

    /* Which cell the point stands in, and the step onto its parent: scale
     * about that cell so the cell fills the whole square again and the point
     * stands where its parent stood. Counted in cells, that scaling is just
     * subtracting the index -- the multiplication by squares and the division
     * by it cancel, which is worth having where a division is software.
     *
     * The cells the border ring encloses are the ones that were thrown away --
     * one at three squares to a side, four at four, nine at five. They have no
     * parent inside the figure, so a point in one is thrown out of the square
     * instead, which is what makes it leave: one pass for the middle of the
     * square, two for the middles inside each cell of the ring, and so on down.
     *
     * Thrown the way its own quarter of the square faces -- the quarters being
     * cut by the diagonals, so which one a point is in is a comparison of the
     * two distances from the middle and their signs. Four times the half-side
     * clears every bailout shape of the same number, the worst being the
     * triangle, whose corners stand at twice its apothem: the nearest this can
     * leave a point is four less the corner of the square itself, which is
     * 2.59 half-sides against the 2 the triangle reaches.
     *
     * Thrown four ways rather than always the same way, which is what the cell
     * beside it amounted to: the square has four-fold symmetry and so should
     * what is drawn on it, and one direction for the whole middle drew it as a
     * ramp from one side to the other with nothing across it.
     *
     * One cell a pass rather than all of them at once: the pass is doing the
     * counting now, so the loop that used to read thirty digits, and the fixed
     * point it read them in, are both gone. */
    int i = (int)tu;
    int j = (int)tv;
    int cut = squares > 2 ? (i > 0 && i < squares - 1 && j > 0 &&
                             j < squares - 1)
                          : (i == 1 && j == 1);
    if (cut) {
        number_t ax = zx < 0 ? -zx : zx, ay = zy < 0 ? -zy : zy;
        number_t throwx = 0, throwy = 0;
        if (ax >= ay)
            throwx = zx < 0 ? -4 * half : 4 * half;
        else
            throwy = zy < 0 ? -4 * half : 4 * half;
        GSL_SET_COMPLEX(&sfvalue(p), zx + throwx, zy + throwy);
        return p;
    }
    GSL_SET_COMPLEX(&sfvalue(p), half * (2 * (tu - i) - 1),
                    half * (2 * (tv - j) - 1));
    return p;
}

/* Whether a point stands under the Koch curve drawn over the segment from
 * (0,0) to (1,0), and if so how many levels down the answer was found.
 *
 * The curve never doubles back, so how far along the segment the point is picks
 * out exactly one of the four smaller curves the big one is made of, and the
 * question repeats in that one's frame. One path down rather than a search:
 * what it costs is the depth, not four to the depth. Falling below the frame is
 * being inside; running out of levels without falling below it is being
 * outside. */
static inline int koch_under(double x, double y, int depth)
{
    /* Walked in double whatever the picture is drawn in. The figure stands at
     * a fixed size in the plane, so what comes in here is a fraction of one
     * edge and a double holds fifty-two bits of it -- thirty-three levels'
     * worth, where twenty-four are drawn. At 113 bits of mantissa every
     * multiplication below would otherwise be software, and the walk cost more
     * than the noise it sits beside; in double it costs less. */
    const double third = 1.0 / 3;
    const double apex = 0.28867513459481288225; /* sqrt(3)/6, the bump height */
    /* The curve never rises above the triangle on its own chord whose apex is
     * that bump -- that triangle is its convex hull -- so a point above either
     * of the two sides is above the whole curve and is outside now rather than
     * in another twenty levels. Without this the walk ran to the bottom for
     * every point outside the figure, which in a picture of a figure is most
     * of them, and cost twice what the noise costs at 113 bits of mantissa
     * where every multiplication is software. */
    const double slope = apex * 2; /* the hull rises this fast from each end */
    for (int level = 1; level <= depth; level += 1) {
        if (y < 0)
            return level; /* under the curve, and this is how far down */
        if (x < 0 || x > 1)
            return 0; /* out of the frame: ground, and known now */
        if (y > x * slope || y > (1 - x) * slope)
            return 0; /* above the hull: ground, and known now */
        if (x < third) {
            x *= 3;
            y *= 3;
        } else if (x < 0.5) {
            /* the side that climbs away at sixty degrees */
            double ax = x - third, ay = y;
            x = 3 * (ax / 2 + ay * KOCH_SIN60);
            y = 3 * (ay / 2 - ax * KOCH_SIN60);
        } else if (x < 2 * third) {
            /* and the one that comes back down */
            double ax = x - 0.5, ay = y - apex;
            x = 3 * (ax / 2 - ay * KOCH_SIN60);
            y = 3 * (ax * KOCH_SIN60 + ay / 2);
        } else {
            x = 3 * (x - 2 * third);
            y *= 3;
        }
    }
    return KOCH_UNDECIDED; /* out of levels: not yet known either way */
}

/**
 * @brief The Koch snowflake, read from a hexagon out.
 * @details snowflake(radius, kaleidoscope, mode) draws a Koch snowflake with
 * its six points standing on the six corners of the hexagon a hexagonal
 * bailout of that number draws. The last two are the kaleidoscope the noise
 * family takes, one and nought by default; see fig_point.
 *
 * A snowflake comes apart exactly: a regular hexagon at the middle, six
 * triangles of the hexagon's own side standing on its six sides, twelve of a
 * third that on the free edges of those, forty-eight of a ninth on theirs, and
 * so on for ever. Three of the six are the corners of the triangle the figure
 * grew from and three are the first bumps put on its edges, and nothing tells
 * them apart -- from the hexagon out the figure has six-fold symmetry, which is
 * what makes this the natural way to read it and what lets one piece of
 * geometry serve all six sectors.
 *
 * That reading is what the picture wants as well. The bands it gives are even:
 * the hexagon and the six triangles are five twelfths of the figure each and
 * the twelve a further tenth, where reading the figure from its first triangle
 * made one band five eighths of the picture and left the rest to a fringe --
 * one flat triangle, with no snowflake to look at and nothing for the outside
 * colour to say.
 *
 * A triangle standing on a side of the hexagon folds across that side into it,
 * which lands it exactly on one of the six the hexagon is made of. A triangle
 * standing on the free edge of another steps onto that one, blown up by three
 * and turned to face the way it faces. Either way it is one step to the parent,
 * so the pass a point leaves on is the level of the piece it stands in, and the
 * hexagon, having nowhere to go, is thrown out the way its sector faces and
 * leaves on the first -- thrown rather than simply declared gone, so that what
 * it hands back is a point of the plane and the outside colouring modes have
 * something to read there. See FIG_THROW_X.
 *
 * The six corners of ground the figure leaves in its hexagon are no part of it
 * and have no level. A point standing there is turned half about, which lands
 * ground on ground and so never leaves: the space comes out in the inside
 * colour instead of taking a band of the outside one off the figure, and the
 * incolouring modes have a moving point to work with. It costs the membership
 * question once a pass for as long as the pixel is looked at, and that is the
 * price of an empty space that stays empty.
 *
 * Under a circular bailout the six points reach a seventh further out than the
 * number says and are cut off: the honest answer to asking for a figure in a
 * shape it does not fit.
 *
 * @param p The call; the arguments are read right to left, see sfaramN.
 * @return Pointer to the call, the evaluator having no use for the result.
 */
sfarg *sfsnowflake(sfarg *const p)
{
    GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
    if (p->argc > 3)
        return p;
    number_t radius = GSL_REAL(sfarg_or(p, 1, 4, 0));
    if (!(radius > 0))
        return p;

    /* A hexagon stands its corners its apothem over cos thirty degrees from
     * the centre, and a snowflake stands its points as far out as the corners
     * of the triangle it grew from -- so this lands the six points on the six
     * corners of the hexagon a bailout of this number would draw. */
    number_t inv = fig_cached(p, radius, 1);
    number_t scale = FIG_SIN60 * inv;
    number_t zx, zy;
    fig_point(p, 2, &zx, &zy);
    number_t x = zx * scale;
    number_t y = zy * scale;

    /* The six points reach exactly as far as the corners of the triangle the
     * figure grew from, so the whole of it sits inside the circle those corners
     * are on and most of the plane can be turned away in three multiplications.
     */
    if (x * x + y * y > 1) {
        /* past the whole figure, so no part of it: half about, and it never
         * leaves -- see FIG_TURN */
        FIG_TURN(p, -1, 0, zx, zy);
        return p;
    }

    /* How far out the point stands on each of the six ways the middle hexagon
     * faces. Three of them are the negatives of the other three, so four
     * multiplications answer all six, and the largest names the sector. */
    number_t u0 = FIG_SIN60 * x + y / 2;
    number_t u2 = -FIG_SIN60 * x + y / 2;
    number_t far = u0;
    int k = 0;
    if (y > far) {
        far = y;
        k = 1;
    }
    if (u2 > far) {
        far = u2;
        k = 2;
    }
    if (-u0 > far) {
        far = -u0;
        k = 3;
    }
    if (-y > far) {
        far = -y;
        k = 4;
    }
    if (-u2 > far) {
        far = -u2;
        k = 5;
    }
    if (far < FIG_HEX_APOTHEM) {
        /* The hexagon at the middle: no parent to walk to, so it is thrown out
         * the way its sector faces and goes on this pass, carrying where it
         * came from with it for the outside colouring to read. */
        number_t root = radius * inv; /* the square root back, by multiplying */
        GSL_SET_COMPLEX(&sfvalue(p), zx + FIG_THROW_X[k] * root,
                        zy + FIG_THROW_Y[k] * root);
        return p;
    }

    /* that sector turned upright, so that the same geometry serves all six */
    number_t c = FIG_TURN_C[k], sn = FIG_TURN_S[k];
    number_t xr = x * c - y * sn;
    number_t yr = x * sn + y * c;

    number_t nx, ny;
    number_t across = FIG_SQRT3 * (xr < 0 ? -xr : xr);
    if (yr + across <= 1) {
        /* The triangle standing on that side of the hexagon: apex out at the
         * figure's reach, base along the side. Its step is onto the hexagon,
         * and what it lands on is the largest triangle the hexagon holds --
         * corner on corner, three times the area, so the step blows the point
         * up by the square root of three and turns it a twelfth of a turn.
         *
         * Folding it flat across the side it stands on would land it on one of
         * the six the hexagon is made of, which is tidier and was what this
         * did. But a fold is an isometry, and an isometry has no sensitivity to
         * where the point started: written inside a larger formula the figure
         * then handed that formula a rigid picture of itself over the whole of
         * this level, which is five twelfths of it, and the formula drew a flat
         * region there. Every other step these three figures take expands --
         * the gasket doubles, the carpet blows a cell up by the number of
         * cells, the levels below this one blow up by three -- and this one now
         * does too. */
        number_t px = xr, py = yr - (number_t)2 / 3;
        nx = (number_t)3 / 2 * px + FIG_SIN60 * py;
        ny = -FIG_SIN60 * px + (number_t)3 / 2 * py;
    } else {
        /* Deeper, and the walk is the one a triangle has always taken -- in
         * that triangle's own frame, where it has circumradius one and stands
         * at the origin, and where only the two edges that are not glued to the
         * hexagon bear children. */
        number_t px = 3 * xr, py = 3 * yr - 2;
        int edge = -1;
        for (int t = 0; t < 2 && edge < 0; t += 1) {
            int e = FIG_FREE[t];
            number_t rx = px - FIG_VX[e], ry = py - FIG_VY[e];
            /* along the edge, and out from it: the corners are taken in the
             * order that puts the outside of each edge on the positive side */
            number_t along = rx * FIG_AX[e] + ry * FIG_AY[e];
            number_t off = rx * FIG_AY[e] - ry * FIG_AX[e];
            if (off <= 0 || along < 0 || along > 1)
                continue;
            /* One level a pass, and no further.
             *
             * The walk used to run all twenty-four levels in a single
             * evaluation, so a formula written around this figure was handed
             * the whole Koch boundary on its first pass and drew it in full
             * however few passes it was given -- where a gasket, whose answer
             * is one doubling, shows one level more for each pass it is
             * allowed. Tying the depth to the pass puts the two on the same
             * footing: what a picture shows is what its iteration count paid
             * for.
             *
             * Running out of levels is not the same as being outside. A point
             * that is out of the frame or above the hull is ground and is known
             * to be ground now; one that is merely undecided is carried on as
             * part of the figure, climbs a level with the rest of them, and is
             * asked again with one level more to spend. So the figure itself is
             * drawn exactly as it was -- a point in it never needed the deep
             * walk, only the assurance that it was not ground -- and it is the
             * ground that comes in a level at a time. */
            int depth = (int)sffe_iteration + 1;
            if (depth > KOCH_DEPTH)
                depth = KOCH_DEPTH;
            if (koch_under((double)along, (double)off, depth))
                edge = e;
        }
        if (edge < 0) {
            /* The ground, which is no part of the figure: turned half about,
             * which is a motion that keeps it there.
             *
             * It has to be a motion. Handed back where it stood it was a fixed
             * point, and a fixed point is not a bounded orbit but a stopped
             * one: the loop ran every pass to the limit with nothing changing,
             * the incolouring modes that read the pass before this one had
             * nothing to read, and a figure written inside a larger formula
             * gave that formula the point it was already holding -- where the
             * gasket and the carpet move every point they are given.
             *
             * Half about rather than the sixth of a turn the figure's own
             * symmetry would also allow: both land ground on ground, but half
             * about leaves the point inside a square and a circle as well as
             * inside the hexagon, where a sixth of a turn carries it out of the
             * square. It is the motion the other two figures give what is no
             * part of them, each by its own symmetry -- see FIG_TURN. */
            FIG_TURN(p, -1, 0, zx, zy);
            return p;
        }
        /* onto the parent: take the child off the corner its edge faces away
         * from, turn it back the way the parent faces, and blow it up by three
         */
        int j = FIG_OPP[edge];
        number_t qx = px + (number_t)2 / 3 * FIG_VX[j];
        number_t qy = py + (number_t)2 / 3 * FIG_VY[j];
        px = 3 * (qx * FIG_CC[j] + qy * FIG_CS[j]);
        py = 3 * (qy * FIG_CC[j] - qx * FIG_CS[j]);
        nx = px * FIG_THIRD;
        ny = (py + 2) * FIG_THIRD;
    }

    /* the sector turned back, and out of the figure's units into the plane. The
     * root comes back by multiplying the square into its reciprocal, which is
     * two multiplications where dividing by the scale would be a division, and
     * a division at 113 bits of mantissa is software. */
    number_t reach = radius * inv * FIG_ACROSS;
    GSL_SET_COMPLEX(&sfvalue(p), (nx * c + ny * sn) * reach,
                    (ny * c - nx * sn) * reach);
    return p;
}

/* Whether a parsed formula calls one of the functions that watch the orbit.
 *
 * The engine asks because boundary tracing has to be turned off for such a
 * formula: it walks the edge of a region, finds one colour all the way round,
 * and fills the inside without computing it. That holds for a fractal, whose
 * bands really are solid. It does not hold for a noise field, where the
 * inside is whatever the noise says, nor for a trap or a stripe average,
 * where two neighbours that take the same number of passes can still have
 * seen quite different orbits. Left on, some pixels are filled rather than
 * computed and are simply wrong.
 *
 * Walking the operation list rather than searching the text: the text would
 * also match a name that merely contains one of these, and this cannot. */
int sffe_uses_noise(sffe *const parser)
{
    if (parser == NULL)
        return 0;
    for (unsigned int i = 0; i < parser->oprCount; i++)
        if (parser->oprs[i].fnc == sfrandsc ||
            parser->oprs[i].fnc == sfrandscq ||
            parser->oprs[i].fnc == sfrandscp ||
            parser->oprs[i].fnc == sfrandsch ||
            parser->oprs[i].fnc == sfrandsct ||
            parser->oprs[i].fnc == sfrandsctile ||
            parser->oprs[i].fnc == sftrap ||
            parser->oprs[i].fnc == sfstripe)
            return 1;
    return 0;
}

static sfarg *randscq_at(sfarg *const p, unsigned int pass,
                         const cmplx *here)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, 0, pass, here, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        /* past the grid: there is no cell to stand in, so no skew either */
        GSL_SET_COMPLEX(&sfvalue(p), randsc_unit(randsc_hash(cx, cy, h)), 0);
        return sfaram1(p);
    }

    /* how far the value is turned by where in the cell the point stands; nought
     * leaves it where it has always been. See randsc_skewed. */
    cmplx skew = sfarg_or(p, 6, 0, 0);
    /* which of the four things the skew does; see randsc_skew_apply */
    int skewmode = (int)GSL_REAL(sfarg_or(p, 7, 0, 0));
    /* the rosette takes the picture's own symmetry, so it wants the count */
    int wedges = (skewmode & (RANDSC_SKEW_ROSETTE | RANDSC_SKEW_WEDGE))
                     ? (int)GSL_REAL(sfarg_or(p, 4, 1, 0))
                     : 1;
    number_t qre, qim;
    /* the larger of the two distances from the middle, which draws squares */
    number_t qu = nfabs(u - (number_t)1 / 2), qv = nfabs(v - (number_t)1 / 2);
    number_t qang = (skewmode & RANDSC_SKEW_ROSETTE)
                        ? randsc_rosette(u - (number_t)1 / 2,
                                         v - (number_t)1 / 2, wedges)
                        : 0;
    randsc_skewed(randsc_hash(cx, cy, h), 2 * (qu > qv ? qu : qv), qang, skew,
                  skewmode, wedges, &qre, &qim);
    GSL_SET_COMPLEX(&sfvalue(p), qre, qim);
    return sfaram1(p);
}

sfarg *sfrandscq(sfarg *const p)
{
    return randsc_run(p, 0, randscq_at);
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
static sfarg *randsch_at(sfarg *const p, unsigned int pass,
                         const cmplx *here)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, 0, pass, here, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p),
                        randsc_unit(randsc_hash(cx, cy, h ^ RANDSCH_SALT)),
                        0);
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

    /* how far the value is turned by where in the cell the point stands; nought
     * leaves it where it has always been. See randsc_skewed. */
    cmplx skew = sfarg_or(p, 6, 0, 0);
    /* which of the four things the skew does; see randsc_skew_apply */
    int skewmode = (int)GSL_REAL(sfarg_or(p, 7, 0, 0));
    /* the rosette takes the picture's own symmetry, so it wants the count */
    int wedges = (skewmode & (RANDSC_SKEW_ROSETTE | RANDSC_SKEW_WEDGE))
                     ? (int)GSL_REAL(sfarg_or(p, 4, 1, 0))
                     : 1;
    number_t hre, him;
    /* how far out of the hexagon's middle, in the hexagon's own reckoning:
     * half the sum of the three cube coordinates' absolute values, which is
     * one at an edge and draws hexagonal contours */
    number_t ha = q - rx, hb = r - rz;
    number_t hd = (nfabs(ha) + nfabs(hb) + nfabs(ha + hb));
    number_t hang =
        (skewmode & RANDSC_SKEW_ROSETTE) ? randsc_rosette(ha, hb, wedges) : 0;
    randsc_skewed(randsc_hash((int64_t)rx, (int64_t)rz, h ^ RANDSCH_SALT),
                  hd > 1 ? 1 : hd, hang, skew, skewmode, wedges, &hre, &him);
    GSL_SET_COMPLEX(&sfvalue(p), hre, him);
    return sfaram1(p);
}

sfarg *sfrandsch(sfarg *const p)
{
    return randsc_run(p, 0, randsch_at);
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
static sfarg *randsct_at(sfarg *const p, unsigned int pass,
                         const cmplx *here)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int state = randsc_setup(p, 0, pass, here, &cx, &cy, &u, &v, &h);
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p),
                        randsc_unit(randsc_hash(cx, cy, h ^ RANDSCT_SALT)),
                        0);
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

    /* how far the value is turned by where in the cell the point stands; nought
     * leaves it where it has always been. See randsc_skewed. */
    cmplx skew = sfarg_or(p, 6, 0, 0);
    /* which of the four things the skew does; see randsc_skew_apply */
    int skewmode = (int)GSL_REAL(sfarg_or(p, 7, 0, 0));
    /* the rosette takes the picture's own symmetry, so it wants the count */
    int wedges = (skewmode & (RANDSC_SKEW_ROSETTE | RANDSC_SKEW_WEDGE))
                     ? (int)GSL_REAL(sfarg_or(p, 4, 1, 0))
                     : 1;
    number_t tre, tim;
    /* the three shares of the triangle the point stands in: the smallest of
     * them is nought on an edge and a third at the middle, so one less three
     * times it runs from nought in the middle to one at an edge, and its
     * contours are triangles */
    number_t s1 = upper ? 1 - fa : fa;
    number_t s2 = upper ? 1 - fb : fb;
    number_t s3 = 1 - s1 - s2;
    number_t sm = s1 < s2 ? s1 : s2;
    if (s3 < sm)
        sm = s3;
    number_t td = 1 - 3 * sm;
    /* the middle of a triangle is where all three shares are a third, so the
     * deviation of two of them says which way round the middle the point is */
    number_t tang = (skewmode & RANDSC_SKEW_ROSETTE)
                        ? randsc_rosette(s1 - (number_t)1 / 3,
                                         s2 - (number_t)1 / 3, wedges)
                        : 0;
    randsc_skewed(randsc_hash(ia, ib,
                              upper ? h ^ RANDSCT_SALT ^ RANDSCT_UPPER
                                    : h ^ RANDSCT_SALT),
                  td < 0 ? 0 : (td > 1 ? 1 : td), tang, skew, skewmode, wedges, &tre,
                  &tim);
    GSL_SET_COMPLEX(&sfvalue(p), tre, tim);
    return sfaram1(p);
}

sfarg *sfrandsct(sfarg *const p)
{
    return randsc_run(p, 0, randsct_at);
}

/* --- randsctile: the same field over any of forty-five tilings -------------
 *
 * randsctile(tiling, seed, ...) is randscq with the squares replaced by the
 * tiles of the tiling its first argument names, and everything after that
 * argument is what randsc takes and means the same. The numbers:
 *
 *   1-3    the regular tilings: squares, triangles, hexagons;
 *   4-11   the eight Archimedean ones, regular polygons with every corner
 *          alike: 4.8.8, 3.6.3.6, 3.4.6.4, 3.12.12, 4.6.12, 3.3.3.4.4,
 *          3.3.4.3.4, 3.3.3.3.6;
 *   12-19  their duals, which are not regular: tetrakis square, rhombille,
 *          deltoidal trihexagonal, triakis triangular, kisrhombille, and the
 *          prismatic, Cairo and floret pentagons;
 *   20-32  bricks, Flemish bond, herringbone, basketweave, Pythagorean,
 *          chevrons, squares and rhombi, houses, rows of squares and triangles
 *          in two rhythms, hexagons among triangles, Greek crosses, T
 *          tetrominoes;
 *   33-37  Islamic stars: of eight points with crosses, of six with hexagons,
 *          of eight from 4.8.8, of twelve from 3.12.12 and from 4.6.12;
 *   38     irregular: randscp's Voronoi cells;
 *   39-43  quasiperiodic, from de Bruijn's multigrids: Penrose's rhombs,
 *          Penrose's kites and darts, Ammann-Beenker, and rhombs of twelve and
 *          of seven directions;
 *   44-45  by substitution: the pinwheel and the chair.
 *
 * Every tiling is scaled to a tile of unit area on average, which is what the
 * rest of the family lays over each unit of the degraded size, so changing
 * the first argument changes the shape of the cells and not their scale. A
 * number outside 1 to 45 draws nought, as a zero size does.
 *
 * The periodic ones are tables (randsctile_tables.h), built and checked by
 * tools/randsctile-tables.py: the point is taken to the lattice of the
 * tiling, split into a whole number of periods -- exact, as randscq's cell is
 * -- and a place within one, and the tiles that can hold a place of one
 * period are tried in turn. Tiles may be concave, the stars and the crosses,
 * so the test is the crossing number rather than a side at a time.
 *
 * The ones that do not repeat have no table to look in and are located
 * afresh at every point, which makes them the dearest of the family. They
 * are worked in double, in both builds -- in long double they took three
 * times as long, and in quad the forty-eight levels of a Penrose tile are
 * software arithmetic -- and so is the place within a period of the ones
 * that do. What a tile is does not need more: the two builds disagree about
 * it only on a hairline along its edges, as they do about any of the
 * mosaics. Past 2^32 cells from the origin, where double would start to lose
 * the tiles that do not repeat, those go flat, as the rest of the family does
 * past its grid.
 *
 * The skew measures a tile as randscp measures its polygon: how far the point
 * is from the tile's edge against the most any point of the tile is -- the
 * radius of the largest circle it holds -- nought there and one along the
 * whole edge, so the contours are the tile shrunk. Against the distance of
 * the tile's middle instead, as it was first written, a concave tile went
 * flat: the middle of an L or a dart sits close to its notch, and the rest
 * of the tile was further from the edge than that and read nought. The
 * rosette turns about the tile's middle, in the tile's own frame.
 */
#define RANDSCTILE_KINDS 45
#define RANDSCTILE_SALT 0x3C6EF372FE94F82BULL

/* past this many cells from the origin a tiling that does not repeat is flat */
#define RANDSCTILE_FAR 4294967296.0

enum {
    RANDSCTILE_TABLE,
    RANDSCTILE_VORONOI,
    RANDSCTILE_GRID_OF,
    RANDSCTILE_KITES,
    RANDSCTILE_PINWHEEL,
    RANDSCTILE_CHAIR
};

/* How each number is drawn, and which of its kind it is: a row of the
 * periodic table, or a multigrid. */
static const struct {
    unsigned char how, which;
} RANDSCTILE_KIND[RANDSCTILE_KINDS] = {
    {RANDSCTILE_TABLE, 0},    {RANDSCTILE_TABLE, 1},
    {RANDSCTILE_TABLE, 2},    {RANDSCTILE_TABLE, 3},
    {RANDSCTILE_TABLE, 4},    {RANDSCTILE_TABLE, 5},
    {RANDSCTILE_TABLE, 6},    {RANDSCTILE_TABLE, 7},
    {RANDSCTILE_TABLE, 8},    {RANDSCTILE_TABLE, 9},
    {RANDSCTILE_TABLE, 10},   {RANDSCTILE_TABLE, 11},
    {RANDSCTILE_TABLE, 12},   {RANDSCTILE_TABLE, 13},
    {RANDSCTILE_TABLE, 14},   {RANDSCTILE_TABLE, 15},
    {RANDSCTILE_TABLE, 16},   {RANDSCTILE_TABLE, 17},
    {RANDSCTILE_TABLE, 18},   {RANDSCTILE_TABLE, 19},
    {RANDSCTILE_TABLE, 20},   {RANDSCTILE_TABLE, 21},
    {RANDSCTILE_TABLE, 22},   {RANDSCTILE_TABLE, 23},
    {RANDSCTILE_TABLE, 24},   {RANDSCTILE_TABLE, 25},
    {RANDSCTILE_TABLE, 26},   {RANDSCTILE_TABLE, 27},
    {RANDSCTILE_TABLE, 28},   {RANDSCTILE_TABLE, 29},
    {RANDSCTILE_TABLE, 30},   {RANDSCTILE_TABLE, 31},
    {RANDSCTILE_TABLE, 32},   {RANDSCTILE_TABLE, 33},
    {RANDSCTILE_TABLE, 34},   {RANDSCTILE_TABLE, 35},
    {RANDSCTILE_TABLE, 36},   {RANDSCTILE_VORONOI, 0},
    {RANDSCTILE_GRID_OF, 0},  {RANDSCTILE_KITES, 0},
    {RANDSCTILE_GRID_OF, 1},  {RANDSCTILE_GRID_OF, 2},
    {RANDSCTILE_GRID_OF, 3},  {RANDSCTILE_PINWHEEL, 0},
    {RANDSCTILE_CHAIR, 0}};

/* What locating a point in a tiling hands back: the hash that names the tile,
 * and, when the skew is asked for, how far out of the tile's middle the point
 * stands and where it stands seen from that middle. */
struct randsctile_hit {
    uint64_t name;
    double out;
    double rx, ry;
};

/* floor without the library call, for the numbers this works in: well inside
 * what an int64 holds, which RANDSCTILE_FAR sees to */
static inline double randsctile_floor(double x)
{
    double t = (double)(int64_t)x;
    return t > x ? t - 1 : t;
}

/* nought where the tile is deepest, one at its edge: the distance to the edge
 * against the most there is */
static double randsctile_out(double edge, double reach)
{
    double o = 1 - edge / reach;
    return o < 0 ? 0 : (o > 1 ? 1 : o);
}

/* The distance from a point to the nearest side of a polygon, the squares
 * compared and one root taken. */
static double randsctile_edge(const double *x, const double *y, int n,
                              double px, double py)
{
    double best = -1;
    for (int i = 0; i < n; i++) {
        int k = i + 1 == n ? 0 : i + 1;
        double dx = x[k] - x[i], dy = y[k] - y[i];
        double t = ((px - x[i]) * dx + (py - y[i]) * dy) / (dx * dx + dy * dy);
        t = t < 0 ? 0 : (t > 1 ? 1 : t);
        double ex = x[i] + t * dx - px, ey = y[i] + t * dy - py;
        double d = ex * ex + ey * ey;
        if (best < 0 || d < best)
            best = d;
    }
    return sqrt(best);
}

/* Which side of the line through a and b the point is on: positive to the
 * left. */
static inline double randsctile_side(double px, double py, double ax,
                                     double ay, double bx, double by)
{
    return (bx - ax) * (py - ay) - (by - ay) * (px - ax);
}

/* Whether p is on the same side of the line through a and b as q is -- the
 * one test the substitutions below need, a child being cut from the rest of
 * its parent by one line at a time. */
static inline int randsctile_with(double px, double py, double qx, double qy,
                                  double ax, double ay, double bx, double by)
{
    double s = randsctile_side(px, py, ax, ay, bx, by);
    double t = randsctile_side(qx, qy, ax, ay, bx, by);
    return (s >= 0) == (t >= 0);
}

/* The periodic tilings. The point goes to lattice coordinates in number_t and
 * is split there, so the whole number of periods is exact however far out it
 * is; only the place within one period, a number of order one, is worked in
 * double. */
static int randsctile_table(int which, number_t X, number_t Y, uint64_t h,
                            int want, struct randsctile_hit *hit)
{
    const auto &t = RANDSCTILE_PERIODIC[which];
    int64_t ia, ib;
    number_t fa, fb;
    if (!randsc_cell(X * (number_t)t.ia + Y * (number_t)t.ib, &ia, &fa) ||
        !randsc_cell(X * (number_t)t.ja + Y * (number_t)t.jb, &ib, &fb)) {
        /* past what a lattice index holds, a little before randsc_setup's own
         * guard where a period is wider than a cell: flat, as there */
        hit->name = randsc_hash(INT64_MIN, INT64_MIN, h);
        return 1;
    }
    double la = (double)fa, lb = (double)fb;
    double px = la * t.ax + lb * t.bx, py = la * t.ay + lb * t.by;

    /* Two passes. The first asks each tile whether it holds the point. Two
     * tiles that share a side hold it as their own vertices say, and those
     * are the same side only to the last bit or so -- the one tile's corners
     * were written out moved by a period, the other's are moved by it here --
     * so a point exactly on a side can find neither. The second pass, which
     * nothing but such a point ever reaches, gives it to the first tile it is
     * within a billionth of a cell of. A hole wider than that is a hole, and
     * finds nothing. */
    const double wide = 1e-9;
    for (int pass = 0; pass < 2; pass++)
        for (int c = t.first; c < t.first + t.count; c++) {
            const auto &cand = RANDSCTILE_CAND[c];
            const auto &f = RANDSCTILE_FACE[cand.face];
            /* the place, moved back by the periods the tile was moved by */
            double qx = px - cand.i * t.ax - cand.j * t.bx;
            double qy = py - cand.i * t.ay - cand.j * t.by;
            /* the box is only there to save the test below, so it may be a
             * little wide and must not be a little narrow */
            if (qx < f.x0 - wide || qx > f.x1 + wide || qy < f.y0 - wide ||
                qy > f.y1 + wide)
                continue;
            const double(*v)[2] = RANDSCTILE_VERT + f.first;
            int in = 0;
            if (pass == 0) {
                for (int i = 0; i < f.count; i++) {
                    int k = i + 1 == f.count ? 0 : i + 1;
                    if ((v[i][1] > qy) != (v[k][1] > qy) &&
                        qx < v[i][0] + (qy - v[i][1]) * (v[k][0] - v[i][0]) /
                                           (v[k][1] - v[i][1]))
                        in = !in;
                }
            }
            double vx[24], vy[24];
            if (pass == 1 || want)
                for (int i = 0; i < f.count; i++) {
                    vx[i] = v[i][0];
                    vy[i] = v[i][1];
                }
            if (pass == 1)
                in = randsctile_edge(vx, vy, f.count, qx, qy) < wide;
            if (!in)
                continue;
            hit->name = randsc_hash(ia + cand.i, ib + cand.j,
                                    h ^ (uint64_t)(cand.face + 1) *
                                            0xD6E8FEB86659FD93ULL);
            if (want) {
                hit->out = randsctile_out(
                    randsctile_edge(vx, vy, f.count, qx, qy), f.reach);
                hit->rx = qx - f.mx;
                hit->ry = qy - f.my;
            }
            return 1;
        }
    return 0;
}

/* de Bruijn's multigrid. In grid space there are n families of parallel
 * lines, x.e_m + g_m a whole number; where a line a of family j crosses a
 * line b of family k there is a tile, the rhombus with sides e_j and e_k
 * whose first corner is a e_j + b e_k plus, for every other family, the whole
 * number above x.e_m + g_m at the crossing. Its name is the four numbers j,
 * k, a, b -- exact, and the same whoever asks.
 *
 * Going the other way, from a point p of the tiling to its tile, uses that
 * the corners sum n/2 times the grid point, give or take a line of every
 * family: the crossing is one of the lines each side of that guess, in some
 * pair of families. Measured over eighteen thousand points as far out as a
 * million, it always was. Tried in order, that is half of up to eighty-four
 * candidates on average; the families whose lines pass nearest the guess are
 * tried first, and each family's nearer line first, which brings it to four
 * to eight. Should the crossing ever not be among them, one line more each
 * way is tried, in plain order, before giving up. */
static int randsctile_grid_try(const double (*e)[2], const double *g, int n,
                               int j, int k, int64_t a, int64_t b, double px,
                               double py, uint64_t h, int want,
                               struct randsctile_hit *hit)
{
    double ejx = e[j][0], ejy = e[j][1], ekx = e[k][0], eky = e[k][1];
    double det = ejx * eky - ejy * ekx;
    double ra = (double)a - g[j], rb = (double)b - g[k];
    double xs = (ra * eky - rb * ejy) / det, ys = (rb * ejx - ra * ekx) / det;
    double vx = (double)a * ejx + (double)b * ekx;
    double vy = (double)a * ejy + (double)b * eky;
    for (int m = 0; m < n; m++) {
        if (m == j || m == k)
            continue;
        double K = -randsctile_floor(-(xs * e[m][0] + ys * e[m][1] + g[m]));
        vx += K * e[m][0];
        vy += K * e[m][1];
    }
    double dx = px - vx, dy = py - vy;
    double s = (dx * eky - dy * ekx) / det, t = (ejx * dy - ejy * dx) / det;
    if (!(s >= 0 && s < 1 && t >= 0 && t < 1))
        return 0;
    hit->name = randsc_hash(
        a, b, h ^ (uint64_t)(j * RANDSCTILE_GRID_MAX + k + 1) * 0xD6E8FEB86659FD93ULL);
    if (want) {
        /* in the rhombus's own coordinates the edge is the nearest of s,
         * 1 - s, t and 1 - t, and the middle is a half from all four */
        double o = s < 1 - s ? s : 1 - s;
        if (t < o)
            o = t;
        if (1 - t < o)
            o = 1 - t;
        hit->out = randsctile_out(o, 0.5);
        hit->rx = s - 0.5;
        hit->ry = t - 0.5;
    }
    return 1;
}

static int randsctile_grid(int which, double X, double Y, uint64_t h, int want,
                           struct randsctile_hit *hit)
{
    const auto &G = RANDSCTILE_GRID[which];
    int n = G.n;
    double px = X / G.edge, py = Y / G.edge;
    double gx = 0, gy = 0;
    for (int m = 0; m < n; m++) {
        gx += (G.g[m] + 0.5) * G.e[m][0];
        gy += (G.g[m] + 0.5) * G.e[m][1];
    }
    double x0 = (px - gx) * 2 / n, y0 = (py - gy) * 2 / n;
    int64_t fl[RANDSCTILE_GRID_MAX], near[RANDSCTILE_GRID_MAX];
    double d[RANDSCTILE_GRID_MAX];
    int order[RANDSCTILE_GRID_MAX];
    for (int m = 0; m < n; m++) {
        double t = x0 * G.e[m][0] + y0 * G.e[m][1] + G.g[m];
        double f = randsctile_floor(t);
        fl[m] = (int64_t)f;
        near[m] = t - f > 0.5 ? fl[m] + 1 : fl[m];
        d[m] = t - f < 0.5 ? t - f : 1 - (t - f);
        /* by how near a line of the family passes, nearest first */
        int i = m;
        while (i > 0 && d[order[i - 1]] > d[m]) {
            order[i] = order[i - 1];
            i--;
        }
        order[i] = m;
    }
    for (int i1 = 1; i1 < n; i1++)
        for (int i2 = 0; i2 < i1; i2++) {
            int j = order[i1] < order[i2] ? order[i1] : order[i2];
            int k = order[i1] < order[i2] ? order[i2] : order[i1];
            int64_t aa[2] = {near[j], 2 * fl[j] + 1 - near[j]};
            int64_t bb[2] = {near[k], 2 * fl[k] + 1 - near[k]};
            for (int q = 0; q < 4; q++)
                if (randsctile_grid_try(G.e, G.g, n, j, k, aa[q >> 1], bb[q & 1],
                                        px, py, h, want, hit))
                    return 1;
        }
    for (int j = 0; j < n; j++)
        for (int k = j + 1; k < n; k++)
            for (int64_t a = fl[j] - 1; a <= fl[j] + 2; a++)
                for (int64_t b = fl[k] - 1; b <= fl[k] + 2; b++)
                    if (randsctile_grid_try(G.e, G.g, n, j, k, a, b, px, py, h,
                                            want, hit))
                        return 1;
    return 0;
}

/* Penrose's kites and darts, from Robinson's triangles: a wheel of ten round
 * the origin, each triangle cut at every level into two or three a golden
 * ratio smaller, and the point followed into the one it falls in. The half
 * kite has its 36 degree corner second, as the cutting rules want it, and a
 * kite or a dart is two halves mirrored across the side B C. The children of
 * a half kite are cut from one another by the lines Q R and A R, those of a
 * half dart by B P, so a level costs a line or two rather than a triangle
 * apiece.
 *
 * A tile is named by the middle of that side, which both halves reach by
 * different routes and so by different arithmetic: it is rounded to a
 * sixty-fourth of a cell first, where the two routes differ by a few parts
 * in a million at the farthest the tiling reaches. */
#define RANDSCTILE_KITE_LEVELS 48
static const double RANDSCTILE_PHI = (1 + sqrt(5.0)) / 2;
/* the wheel, as far out as the tiling reaches after all its levels */
static const double RANDSCTILE_WHEEL = pow(RANDSCTILE_PHI, RANDSCTILE_KITE_LEVELS);

static int randsctile_kites(double X, double Y, uint64_t h, int want,
                            struct randsctile_hit *hit)
{
    const double IPHI = RANDSCTILE_PHI - 1; /* its reciprocal */
    const double R = RANDSCTILE_WHEEL;
    double px = X / RANDSCTILE_KITES_EDGE, py = Y / RANDSCTILE_KITES_EDGE;

    int i = (int)randsctile_floor((atan2(py, px) / M_PI * 10 + 1) / 2);
    i = ((i % 10) + 10) % 10;
    double ax = R * cos((2 * i - 1) * M_PI / 10), ay = R * sin((2 * i - 1) * M_PI / 10);
    double cx = R * cos((2 * i + 1) * M_PI / 10), cy = R * sin((2 * i + 1) * M_PI / 10);
    if (i % 2 == 0) {
        double t = ax;
        ax = cx;
        cx = t;
        t = ay;
        ay = cy;
        cy = t;
    }
    double bx = 0, by = 0;
    int red = 1;

    for (int level = 0; level < RANDSCTILE_KITE_LEVELS; level++) {
        if (red) {
            double qx = ax + (bx - ax) * IPHI, qy = ay + (by - ay) * IPHI;
            double rx = bx + (cx - bx) * IPHI, ry = by + (cy - by) * IPHI;
            if (randsctile_with(px, py, bx, by, qx, qy, rx, ry)) {
                /* R Q B, a half dart */
                ax = rx;
                ay = ry;
                cx = bx;
                cy = by;
                bx = qx;
                by = qy;
                red = 0;
            } else if (randsctile_with(px, py, qx, qy, ax, ay, rx, ry)) {
                /* Q A R */
                cx = rx;
                cy = ry;
                bx = ax;
                by = ay;
                ax = qx;
                ay = qy;
            } else {
                /* C A R */
                double tx = cx, ty = cy;
                bx = ax;
                by = ay;
                ax = tx;
                ay = ty;
                cx = rx;
                cy = ry;
            }
        } else {
            double qx = cx + (ax - cx) * IPHI, qy = cy + (ay - cy) * IPHI;
            if (randsctile_with(px, py, ax, ay, bx, by, qx, qy)) {
                /* B P A, a half dart */
                double tx = ax, ty = ay;
                ax = bx;
                ay = by;
                bx = qx;
                by = qy;
                cx = tx;
                cy = ty;
            } else {
                /* P C B, a half kite */
                double tx = bx, ty = by;
                ax = qx;
                ay = qy;
                bx = cx;
                by = cy;
                cx = tx;
                cy = ty;
                red = 1;
            }
        }
    }

    double mx = (bx + cx) / 2, my = (by + cy) / 2;
    hit->name = randsc_hash((int64_t)randsctile_floor(mx * 64 + 0.5),
                            (int64_t)randsctile_floor(my * 64 + 0.5),
                            h ^ (red ? 0 : 0xA5A5A5A5A5A5A5A5ULL));
    if (want) {
        /* the other half is this one mirrored across B C */
        double ux = cx - bx, uy = cy - by;
        double t = ((ax - bx) * ux + (ay - by) * uy) / (ux * ux + uy * uy);
        double fx = 2 * (bx + t * ux) - ax, fy = 2 * (by + t * uy) - ay;
        double vx[4] = {ax, bx, fx, cx}, vy[4] = {ay, by, fy, cy};
        double gx = (ax + fx + 2 * bx + 2 * cx) / 6;
        double gy = (ay + fy + 2 * by + 2 * cy) / 6;
        hit->out = randsctile_out(randsctile_edge(vx, vy, 4, px, py),
                                  red ? RANDSCTILE_KITE_REACH
                                      : RANDSCTILE_DART_REACH);
        hit->rx = px - gx;
        hit->ry = py - gy;
    }
    return 1;
}

/* Conway and Radin's pinwheel: the right triangle with legs one and two cut
 * into five like it, root five smaller -- the one the foot of the altitude
 * cuts off, and the other part, which is the triangle half as large again,
 * cut by its midpoints into four. The children turn by an angle that is no
 * fraction of a turn, so the tiles face every way there is. The altitude
 * parts the first child from the rest, and the lines between the midpoints
 * part the three corners from the one in the middle.
 *
 * A triangle is A, B, C with the right angle at B and A at the end of the
 * long leg. Every point of a tile reached it by the one route, so the tile's
 * middle is the same number wherever it is asked from and names it. */
#define RANDSCTILE_PINWHEEL_LEVELS 29
static const double RANDSCTILE_PINWHEEL_SIZE =
    pow(sqrt(5.0), RANDSCTILE_PINWHEEL_LEVELS);

static int randsctile_pinwheel(double X, double Y, uint64_t h, int want,
                               struct randsctile_hit *hit)
{
    const double L = RANDSCTILE_PINWHEEL_SIZE;
    double px = X, py = Y;
    /* the patch: two of them back to back, a rectangle two by one */
    double ax = -L, ay = -L / 2, bx = L, by = -L / 2, cx = L, cy = L / 2;
    if (!randsctile_with(px, py, bx, by, ax, ay, cx, cy)) {
        ax = L;
        ay = L / 2;
        bx = -L;
        by = L / 2;
        cx = -L;
        cy = -L / 2;
    }
    for (int level = 0; level < RANDSCTILE_PINWHEEL_LEVELS; level++) {
        double hx = (ax + 4 * cx) / 5, hy = (ay + 4 * cy) / 5;
        if (randsctile_with(px, py, cx, cy, bx, by, hx, hy)) {
            /* B H C, the part the altitude cuts off */
            ax = bx;
            ay = by;
            bx = hx;
            by = hy;
            continue;
        }
        double abx = (ax + bx) / 2, aby = (ay + by) / 2;
        double bhx = (bx + hx) / 2, bhy = (by + hy) / 2;
        double ahx = (ax + hx) / 2, ahy = (ay + hy) / 2;
        if (randsctile_with(px, py, ax, ay, ahx, ahy, abx, aby)) {
            /* A, the middle of A H, the middle of A B */
            bx = ahx;
            by = ahy;
            cx = abx;
            cy = aby;
        } else if (randsctile_with(px, py, hx, hy, ahx, ahy, bhx, bhy)) {
            ax = ahx;
            ay = ahy;
            bx = hx;
            by = hy;
            cx = bhx;
            cy = bhy;
        } else if (randsctile_with(px, py, bx, by, abx, aby, bhx, bhy)) {
            cx = bx;
            cy = by;
            ax = abx;
            ay = aby;
            bx = bhx;
            by = bhy;
        } else {
            /* the one in the middle */
            ax = bhx;
            ay = bhy;
            bx = abx;
            by = aby;
            cx = ahx;
            cy = ahy;
        }
    }
    double gx = (ax + bx + cx) / 3, gy = (ay + by + cy) / 3;
    hit->name = randsc_hash((int64_t)randsctile_floor(gx * 16),
                            (int64_t)randsctile_floor(gy * 16), h);
    if (want) {
        double vx[3] = {ax, bx, cx}, vy[3] = {ay, by, cy};
        /* the inradius of legs one and two: (1 + 2 - root 5) / 2 */
        hit->out = randsctile_out(randsctile_edge(vx, vy, 3, px, py),
                                  (3 - sqrt(5.0)) / 2);
        hit->rx = px - gx;
        hit->ry = py - gy;
    }
    return 1;
}

/* The chair: an L of three squares cut into four Ls half its size, one in each
 * corner of it and one in the middle, the two at the ends of the arms turned
 * a quarter each way. The point is carried down in the coordinates of the L
 * it is in -- [0,2]^2 without its upper right quarter -- so the numbers stay
 * of order one all the way down, and the route taken names the tile. */
#define RANDSCTILE_CHAIR_LEVELS 34
/* the side of the smallest L's box, three quarters of it squared being the
 * unit area of a tile */
static const double RANDSCTILE_CHAIR_SIDE = 2 / sqrt(3.0);

static int randsctile_chair(double X, double Y, uint64_t h, int want,
                            struct randsctile_hit *hit)
{
    const double S = ldexp(1, RANDSCTILE_CHAIR_LEVELS);
    /* the lower left corner of each child's box, in halves of the parent's
     * coordinates, and the quarter turns it is turned by */
    static const int DX[4] = {0, 1, 2, 0}, DY[4] = {0, 1, 0, 2};
    static const int TURN[4] = {0, 0, 1, 3};
    /* the origin in the middle of the lower left square of the whole L */
    double x = (X / RANDSCTILE_CHAIR_SIDE + S / 4) / (S / 2);
    double y = (Y / RANDSCTILE_CHAIR_SIDE + S / 4) / (S / 2);
    uint64_t lo = 0, hi = 0;
    for (int level = 0; level < RANDSCTILE_CHAIR_LEVELS; level++) {
        /* the last is taken when none of the first three holds the point,
         * which is only ever a matter of the last bit on an edge */
        int pick = 3;
        double nx = x, ny = y;
        for (int c = 0; c < 4; c++) {
            double u = (x - (double)DX[c] / 2) * 2 - 1;
            double v = (y - (double)DY[c] / 2) * 2 - 1;
            for (int q = 0; q < TURN[c]; q++) {
                double w = u;
                u = v;
                v = -w;
            }
            u += 1;
            v += 1;
            if (c == 3 ||
                (u >= 0 && u < 2 && v >= 0 && v < 2 && !(u >= 1 && v >= 1))) {
                pick = c;
                nx = u;
                ny = v;
                break;
            }
        }
        x = nx;
        y = ny;
        if (level < 32)
            lo = lo << 2 | (uint64_t)pick;
        else
            hi = hi << 2 | (uint64_t)pick;
    }
    hit->name = randsc_hash((int64_t)lo, (int64_t)hi, h);
    if (want) {
        static const double LX[6] = {0, 2, 2, 1, 1, 0};
        static const double LY[6] = {0, 0, 1, 1, 2, 2};
        /* the deepest an L of side two goes is a half, in the middle of
         * any of its three squares; its middle is (5/6, 5/6) */
        hit->out = randsctile_out(randsctile_edge(LX, LY, 6, x, y), 0.5);
        hit->rx = x - 5.0 / 6;
        hit->ry = y - 5.0 / 6;
    }
    return 1;
}

static sfarg *randsctile_at(sfarg *const p, unsigned int pass,
                            const cmplx *here)
{
    int64_t cx, cy;
    number_t u, v;
    uint64_t h;

    int kind = p->argc >= 1 ? (int)GSL_REAL(sfarg_or(p, 1, 0, 0)) : 0;
    int state = kind >= 1 && kind <= RANDSCTILE_KINDS
                    ? randsc_setup(p, 1, pass, here, &cx, &cy, &u, &v, &h)
                    : RANDSC_STOP;
    if (state == RANDSC_STOP) {
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }
    /* a salt apiece, so that two tilings -- and this and the five -- do not
     * hand back the same value where their indices happen to meet */
    uint64_t salt = h ^ RANDSCTILE_SALT ^ (uint64_t)kind * 0x9E3779B97F4A7C15ULL;
    if (state == RANDSC_BEYOND) {
        GSL_SET_COMPLEX(&sfvalue(p), randsc_unit(randsc_hash(cx, cy, salt)), 0);
        return sfaram1(p);
    }

    /* everything after the tiling one place further along than in randsc */
    cmplx skew = sfarg_or(p, 7, 0, 0);
    int skewmode = (int)GSL_REAL(sfarg_or(p, 8, 0, 0));
    int wedges = (skewmode & (RANDSC_SKEW_ROSETTE | RANDSC_SKEW_WEDGE))
                     ? (int)GSL_REAL(sfarg_or(p, 5, 1, 0))
                     : 1;
    int want = GSL_REAL(skew) != 0 || GSL_IMAG(skew) != 0;

    /* the position in cells of unit area, as randsch and randsct take it */
    number_t X = (number_t)cx + u, Y = (number_t)cy + v;
    double lx = (double)X, ly = (double)Y;
    int how = RANDSCTILE_KIND[kind - 1].how, which = RANDSCTILE_KIND[kind - 1].which;
    int far = how != RANDSCTILE_TABLE && how != RANDSCTILE_VORONOI &&
              !(fabs(lx) < RANDSCTILE_FAR && fabs(ly) < RANDSCTILE_FAR);
    struct randsctile_hit hit = {0, 0, 0, 0};
    int found = 0;
    if (far) {
        /* past the reach of a tiling that does not repeat: flat, as the rest
         * of the family is past its grid */
        GSL_SET_COMPLEX(&sfvalue(p), randsc_unit(randsc_hash(0, 0, salt)), 0);
        return sfaram1(p);
    }
    switch (how) {
        case RANDSCTILE_TABLE:
            found = randsctile_table(which, X, Y, salt, want, &hit);
            break;
        case RANDSCTILE_VORONOI: {
            number_t pout, pbx, pby;
            hit.name = randsc_remix(
                randscp_nearest(cx, cy, u, v, salt, want, &pout, &pbx, &pby));
            hit.out = (double)pout;
            hit.rx = -(double)pbx;
            hit.ry = -(double)pby;
            found = 1;
            break;
        }
        case RANDSCTILE_GRID_OF:
            found = randsctile_grid(which, lx, ly, salt, want, &hit);
            break;
        case RANDSCTILE_KITES:
            found = randsctile_kites(lx, ly, salt, want, &hit);
            break;
        case RANDSCTILE_PINWHEEL:
            found = randsctile_pinwheel(lx, ly, salt, want, &hit);
            break;
        case RANDSCTILE_CHAIR:
            found = randsctile_chair(lx, ly, salt, want, &hit);
            break;
    }
    if (!found) {
        /* No tile claims the point. The tables are checked never to leave
         * one, and the multigrid never has, so this is a hole in one of them
         * rather than something a picture should ever show: nought, which the
         * tests look for. */
        GSL_SET_COMPLEX(&sfvalue(p), 0, 0);
        return sfaram1(p);
    }

    number_t tre, tim;
    number_t ang = (skewmode & RANDSC_SKEW_ROSETTE)
                       ? randsc_rosette((number_t)hit.rx, (number_t)hit.ry, wedges)
                       : 0;
    randsc_skewed(hit.name, (number_t)hit.out, ang, skew, skewmode, wedges, &tre,
                  &tim);
    GSL_SET_COMPLEX(&sfvalue(p), tre, tim);
    return sfaram1(p);
}

sfarg *sfrandsctile(sfarg *const p)
{
    return randsc_run(p, 1, randsctile_at);
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
