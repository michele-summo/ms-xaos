/* How accurate gamma and lambertw actually are, at whatever precision this is
 * built for.
 *
 * Both were written when number_t was 64 bits of mantissa, and neither carries
 * anything that says so out loud -- a truncated series and a fixed tolerance
 * both keep working at 113 bits, just no better than they did. Only measuring
 * tells the two apart, and the two turned out to differ:
 *
 *   lambertw was already exact to a couple of ulp at quad despite a tolerance
 *   hardcoded at 1e-15, because Halley converges cubically: the step that
 *   brings the error to 1e-15 lands it near 1e-45, and the test that stops the
 *   loop runs after that step has been applied.
 *
 *   gamma was not. The Lanczos coefficients it used are a fixed set chosen for
 *   about fifteen digits, so it held 2.2e-15 at both precisions -- the whole of
 *   a double's accuracy, and a five-thousandth of a quad's. The quad build now
 *   uses Stirling with enough Bernoulli terms instead.
 *
 * The tolerances below therefore differ by implementation and not merely by
 * type: the Lanczos ceiling is a property of those coefficients, so the long
 * double build is held to it rather than to its own epsilon.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#ifdef USE_FLOAT128
#include <quadmath.h>
#endif

#include "config.h"
#include "number_math.h"
#include "sffe.h"
#include "misc-f.h"

const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

#ifdef USE_FLOAT128
/* Stirling, with the shift up to Re(w) >= 50 costing a few dozen roundings in
 * the recurrence that undoes it -- hence 1e-29 rather than the type's 1.9e-34. */
static const number_t GAMMA_TOLERANCE = 1e-29;
static const number_t LAMBERTW_TOLERANCE = 1e-31;
/* erf is a Maclaurin series, whose terms peak near e^(|z|^2) before decaying,
 * so it spends about 0.43|z|^2 digits on cancellation: a few ulp close in and
 * progressively fewer digits further out. Two bands are checked, one inside a
 * bailout of two, where a fractal actually lives, and one twice as far. */
static const number_t ERF_NEAR_TOLERANCE = 1e-32;
static const number_t ERF_FAR_TOLERANCE = 1e-27;
#else
/* The Lanczos ceiling, which is where this sits at any precision. */
static const number_t GAMMA_TOLERANCE = 1e-13;
static const number_t LAMBERTW_TOLERANCE = 1e-17;
static const number_t ERF_NEAR_TOLERANCE = 1e-17;
static const number_t ERF_FAR_TOLERANCE = 1e-12;
#endif

static int failures = 0;

static void check(int ok, const char *what, number_t measured, number_t allowed)
{
    char got[64], limit[64];
    printf("%-6s %-46s %s (limite %s)\n", ok ? "ok" : "FAIL", what,
           xnumtostr(got, sizeof got, measured, 4),
           xnumtostr(limit, sizeof limit, allowed, 2));
    if (!ok)
        failures++;
}

static sffe *compile(const char *expression, sfNumber *z)
{
    sffe *parser = sffe_alloc();
    if (!parser)
        exit(1);
    sffe_regvar(&parser, z, "z");
    if (sffe_parse(&parser, expression)) {
        printf("FAIL   cannot parse %s\n", expression);
        exit(1);
    }
    return parser;
}

static number_t abs2(gsl_complex a)
{
    return GSL_REAL(a) * GSL_REAL(a) + GSL_IMAG(a) * GSL_IMAG(a);
}

/* |a - b| / |b| */
static number_t relative(gsl_complex a, gsl_complex b)
{
    return nsqrt(abs2(gsl_complex_sub(a, b)) / (abs2(b) + 1e-300));
}

int main(void)
{
    sfNumber z;

    /* Gamma(n + 1) = n!, computed alongside as an exact reference: every
     * factorial up to 20! is held exactly by both types. */
    {
        sffe *gamma = compile("gamma(z)", &z);
        number_t worst = 0, factorial = 1;
        for (int n = 1; n <= 20; n++) {
            factorial *= (number_t)n;
            GSL_SET_COMPLEX(&z, (number_t)(n + 1), 0);
            sfNumber r = sffe_eval(gamma);
            number_t got = GSL_REAL(r);
            number_t d = got > factorial ? got - factorial : factorial - got;
            number_t e = d / factorial;
            if (e > worst)
                worst = e;
        }
        check(worst <= GAMMA_TOLERANCE, "gamma against n! for n = 1..20", worst,
              GAMMA_TOLERANCE);
        sffe_free(&gamma);
    }

    /* Over the complex plane there is no closed form to compare against, but
     * two identities hold everywhere gamma is finite. The reflection one is
     * what exercises the Re(z) < 0.5 path, which nothing else here reaches. */
    {
        sffe *gamma = compile("gamma(z)", &z);
        number_t worst_recurrence = 0, worst_reflection = 0;
        for (int i = -30; i <= 40; i++) {
            for (int j = -30; j <= 30; j++) {
                number_t x = (number_t)i / 4, y = (number_t)j / 4;
                /* Poles sit at every non-positive integer, and the reflection
                 * identity brings in 1 - z as well, so on the real axis every
                 * integer has to be skipped. */
                if (y == 0 && x == nfloor(x))
                    continue;

                GSL_SET_COMPLEX(&z, x, y);
                gsl_complex gz = sffe_eval(gamma);
                GSL_SET_COMPLEX(&z, x + 1, y);
                gsl_complex gz1 = sffe_eval(gamma);
                GSL_SET_COMPLEX(&z, 1 - x, -y);
                gsl_complex g1z = sffe_eval(gamma);

                gsl_complex zc;
                GSL_SET_COMPLEX(&zc, x, y);
                number_t e = relative(gsl_complex_mul(zc, gz), gz1);
                if (e > worst_recurrence)
                    worst_recurrence = e;

                gsl_complex pz, pi_c;
                GSL_SET_COMPLEX(&pz, N_PI * x, N_PI * y);
                GSL_SET_COMPLEX(&pi_c, N_PI, 0);
                e = relative(gsl_complex_mul(gz, g1z),
                             gsl_complex_div(pi_c, gsl_complex_sin(pz)));
                if (e > worst_reflection)
                    worst_reflection = e;
            }
        }
        check(worst_recurrence <= GAMMA_TOLERANCE,
              "gamma recurrence G(z+1) = z G(z)", worst_recurrence,
              GAMMA_TOLERANCE);
        check(worst_reflection <= GAMMA_TOLERANCE,
              "gamma reflection G(z)G(1-z) = pi/sin(pi z)", worst_reflection,
              GAMMA_TOLERANCE);
        sffe_free(&gamma);
    }

    /* W(z) e^W(z) = z. This cannot tell one branch from another on its own --
     * every branch satisfies it -- but the branch selection is what the cases
     * in sffe_test.cpp pin down, and what this measures is the accuracy. */
    {
        sffe *lambertw = compile("lambertw(z)", &z);
        number_t worst = 0;
        for (int i = -20; i <= 20; i++) {
            for (int j = -20; j <= 20; j++) {
                if (!i && !j)
                    continue;
                number_t x = (number_t)i / 4, y = (number_t)j / 4;
                GSL_SET_COMPLEX(&z, x, y);
                sfNumber w = sffe_eval(lambertw);
                gsl_complex zc;
                GSL_SET_COMPLEX(&zc, x, y);
                number_t e = relative(
                    gsl_complex_mul(w, gsl_complex_exp(w)), zc);
                if (e > worst)
                    worst = e;
            }
        }
        check(worst <= LAMBERTW_TOLERANCE, "lambertw identity W e^W = z", worst,
              LAMBERTW_TOLERANCE);
        sffe_free(&lambertw);
    }

    /* erf has an exact reference on the real axis in the C library, at every
     * precision this builds at. Off it there is none, but erf(conj z) must be
     * conj(erf z), which a series in z alone satisfies exactly -- anything
     * other than zero would mean the fold of the left half plane is wrong. */
    {
        sffe *erf_fn = compile("erf(z)", &z);
        number_t near_axis = 0, far_axis = 0, symmetry = 0;
        for (int i = -40; i <= 40; i++) {
            number_t x = (number_t)i / 10;
            if (x == 0)
                continue;
            GSL_SET_COMPLEX(&z, x, 0);
            sfNumber r = sffe_eval(erf_fn);
#ifdef USE_FLOAT128
            number_t want = erfq((__float128)x);
#else
            number_t want = (number_t)erfl((long double)x);
#endif
            number_t d = GSL_REAL(r) - want;
            if (d < 0)
                d = -d;
            number_t e = d / (want < 0 ? -want : want);
            if (x >= -2 && x <= 2) {
                if (e > near_axis)
                    near_axis = e;
            } else if (e > far_axis)
                far_axis = e;
        }
        for (int i = -32; i <= 32; i++) {
            for (int j = 1; j <= 32; j++) {
                GSL_SET_COMPLEX(&z, (number_t)i / 8, (number_t)j / 8);
                gsl_complex a = sffe_eval(erf_fn);
                GSL_SET_COMPLEX(&z, (number_t)i / 8, -(number_t)j / 8);
                gsl_complex c = sffe_eval(erf_fn);
                number_t e = relative(gsl_complex_conjugate(a), c);
                if (e > symmetry)
                    symmetry = e;
            }
        }
        check(near_axis <= ERF_NEAR_TOLERANCE, "erf on the real axis, |x| <= 2",
              near_axis, ERF_NEAR_TOLERANCE);
        check(far_axis <= ERF_FAR_TOLERANCE, "erf on the real axis, 2 < |x| <= 4",
              far_axis, ERF_FAR_TOLERANCE);
        check(symmetry == 0, "erf(conj z) = conj(erf z)", symmetry, 0);
        sffe_free(&erf_fn);
    }

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures != 0;
}
