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
 * @details ifiterf(a; b) evaluates a on every pass but the final one, and b on
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
 * @details ifiterr(a; b; n) evaluates a while the pass number is below n and b
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

static RANDSC_INLINE int randsc_setup(sfarg *const p, int64_t *cx, int64_t *cy, number_t *u,
                        number_t *v, uint64_t *hash)
{
    if (p->argc < 1 || p->argc > 5)
        return RANDSC_STOP;

    /* Seed, cell size, degradation, kaleidoscope level and its mode, in the
     * order they are written; all but the seed have a default, which a call
     * takes by stopping short of them or by leaving their place empty.
     *
     * Degradation halves the cells each pass. One, which was the default,
     * leaves them the size they started and so wastes the argument on a call
     * that says nothing; a half is the shrinking one asks for when one asks. */
    cmplx seed = sfarg_or(p, 1, 0, 0);
    cmplx size = sfarg_or(p, 2, 1, 1);
    cmplx degradation = sfarg_or(p, 3, (number_t)1 / 2, (number_t)1 / 2);
    int level = (int)GSL_REAL(sfarg_or(p, 4, 1, 0));
    int mode = (int)GSL_REAL(sfarg_or(p, 5, 0, 0));

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
    number_t px = GSL_REAL(sffe_position), py = GSL_IMAG(sffe_position);
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
/**
 * @brief A polynomial in the first argument, the rest being its coefficients.
 * @details poly(z; k1; k2; ...; km) is
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
 * whole formula is trap(z^2+c; 3) and nothing else: the fractal iterates as it
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
 * @details trap(a; shape; centre; size) measures the distance from a to the
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
 * @details stripe(a; density) averages (sin(density * arg a) + 1) / 2 over the
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
 * left alone either way. */
#define FIG_TURN(p, cs, sn)                                                    \
    GSL_SET_COMPLEX(&sfvalue(p),                                               \
                    (cs) * GSL_REAL(sffe_z) - (sn) * GSL_IMAG(sffe_z),         \
                    (sn) * GSL_REAL(sffe_z) + (cs) * GSL_IMAG(sffe_z))
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
 * @details sierpinskyt(radius) stands the triangle a triangular bailout of
 * that number draws, point upwards, and cuts the gasket out of it: a point
 * leaves on the pass numbered by the cut that took it, and one on the gasket
 * itself never leaves, so the iteration count is the picture.
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
    if (p->argc > 1)
        return p;
    number_t radius = GSL_REAL(sfarg_or(p, 1, 4, 0));
    if (!(radius > 0))
        return p;

    /* The sides stand the square root of the radius from the centre, where a
     * triangular bailout of that number stands its sides, which puts the
     * corners at twice that: a triangle's apothem is half its circumradius. */
    number_t corner = 2 * fig_cached(p, radius, 0);
    number_t x = GSL_REAL(sffe_z);
    number_t y = GSL_IMAG(sffe_z);

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
        FIG_TURN(p, (number_t)-1 / 2, FIG_SIN60);
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
 * @details sierpinskyc(radius, squares) fills the square a square bailout of
 * that number draws -- half a side the square root of radius, at the origin --
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
    if (p->argc > 2)
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
    number_t tu = (GSL_REAL(sffe_z) + half) * invhalf * squares / 2;
    number_t tv = (GSL_IMAG(sffe_z) + half) * invhalf * squares / 2;
    if (!(tu >= 0) || tu >= squares || !(tv >= 0) || tv >= squares) {
        /* out of the square, so no part of the figure: a quarter of a turn,
         * which is the square's own, and it never leaves */
        FIG_TURN(p, 0, 1);
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
        number_t zx = GSL_REAL(sffe_z), zy = GSL_IMAG(sffe_z);
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
 * @details snowflake(radius) draws a Koch snowflake with its six points
 * standing on the six corners of the hexagon a hexagonal bailout of that
 * number draws.
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
    if (p->argc > 1)
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
    number_t x = GSL_REAL(sffe_z) * scale;
    number_t y = GSL_IMAG(sffe_z) * scale;

    /* The six points reach exactly as far as the corners of the triangle the
     * figure grew from, so the whole of it sits inside the circle those corners
     * are on and most of the plane can be turned away in three multiplications.
     */
    if (x * x + y * y > 1) {
        /* past the whole figure, so no part of it: half about, and it never
         * leaves -- see FIG_TURN */
        FIG_TURN(p, -1, 0);
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
        GSL_SET_COMPLEX(&sfvalue(p),
                        GSL_REAL(sffe_z) + FIG_THROW_X[k] * root,
                        GSL_IMAG(sffe_z) + FIG_THROW_Y[k] * root);
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
            FIG_TURN(p, -1, 0);
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
            parser->oprs[i].fnc == sftrap ||
            parser->oprs[i].fnc == sfstripe)
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
