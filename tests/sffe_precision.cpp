/* The formula parser has to evaluate at the precision the build was
 * configured for.
 *
 * It did not: gsl_complex, which carries every value a user-defined formula
 * produces, was a plain double whatever number_t was, and the maths under it
 * called the double entry points. A user formula therefore ran with 53 bits of
 * mantissa in all three builds, and stopped resolving anything at a zoom of
 * about 1e-15 -- while the built-in formulas, written in number_t throughout,
 * kept going as deep as the type allowed. The checksum net in engine_ref only
 * drives the built-in iteration loops, so nothing here would have noticed.
 *
 * Both halves matter and fail differently, so both are checked: that a value
 * survives being carried through the parser, and that the functions applied to
 * it are the ones belonging to the type.
 */

#include <cstdio>
#include <cstring>

#include "config.h"
#include "number_math.h"
#include "sffe.h"
#include "misc-f.h"

const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

/* The step from 1 to the next representable number_t, found rather than looked
 * up so that this says the same thing whatever number_t is today. */
static number_t type_epsilon(void)
{
    number_t e = 1;
    while ((number_t)1 + e / 2 != (number_t)1)
        e /= 2;
    return e;
}

static int failures = 0;

static void check(int ok, const char *what, const char *detail)
{
    printf("%-6s %s%s%s\n", ok ? "ok" : "FAIL", what, detail ? ": " : "",
           detail ? detail : "");
    if (!ok)
        failures++;
}

/* Parses expr with z bound to the given value and returns the real part. */
static number_t eval_at(const char *expr, number_t zre, int *err)
{
    sffe *parser = sffe_alloc();
    sfNumber z;
    *err = 1;
    if (!parser)
        return 0;
    GSL_SET_COMPLEX(&z, zre, 0);
    sffe_regvar(&parser, &z, "z");
    if (sffe_parse(&parser, expr)) {
        sffe_free(&parser);
        return 0;
    }
    sfNumber r = sffe_eval(parser);
    number_t out = GSL_REAL(r);
    sffe_free(&parser);
    *err = 0;
    return out;
}

int main(void)
{
    char b1[80], b2[80];
    const number_t eps = type_epsilon();

    printf("number_t epsilon = %s\n", xnumtostr(b1, sizeof b1, eps, 6));

    /* A value one step above 1 has to come back out of the parser as itself.
     * Rounded to double it would be exactly 1, and the answer would be 0. */
    {
        int err;
        number_t got = eval_at("z-1", (number_t)1 + eps, &err);
        check(!err, "parsed \"z-1\"", NULL);
        if (!err)
            check(got == eps, "one ulp above 1 survives evaluation",
                  xnumtostr(b1, sizeof b1, got, 6));
    }

    /* And the functions have to be the ones belonging to number_t. sin(pi) is
     * zero to whatever precision pi is held and sin is computed at, so the
     * result reports the width of the arithmetic underneath: about 1e-16 if
     * this went through the double entry points, far smaller otherwise. */
    {
        int err;
        number_t got = eval_at("sin(z)", (number_t)N_PI, &err);
        number_t magnitude = got < 0 ? -got : got;
        check(!err, "parsed \"sin(z)\"", NULL);
        if (!err)
            check(magnitude < eps * 16,
                  "sin(pi) lands within a few ulp of zero",
                  xnumtostr(b2, sizeof b2, magnitude, 6));
    }

    /* On a build wider than double, state the part that used to be broken:
     * the extra precision is beyond anything a double could have held. */
#if defined(USE_FLOAT128) || defined(USE_LONG_DOUBLE)
    check(eps < (number_t)2.3e-16,
          "the build resolves more finely than a double",
          xnumtostr(b1, sizeof b1, eps, 6));
#endif

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures != 0;
}
