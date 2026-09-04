/* The history of z a user formula reaches back into.
 *
 * p1 is the value z had on the pass before, p2 the one before that, and the
 * six of them used to be kept by shifting a six-place array along by one every
 * pass. They are kept in a ring now, deep enough for the furthest one the
 * formula names and read only where it names one, which is what lets p9999
 * cost what p6 costs.
 *
 * A ring is not obviously the same thing as a shift, so the shift is written
 * out here and the two are compared value by value: whatever the depth, and
 * from whatever the history starts at, a picture drawn before this change must
 * still draw the same way.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "phist.h"

const char *qt_gettext(const char *, const char *text) { return text; }

static int failures = 0;

static void check(int ok, const char *what)
{
    printf("%-6s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        failures++;
}

/* A parser with the resolver on it, holding one formula, and nothing else --
 * which is all the taps are worked out from. */
static sffe *compile(const char *expr)
{
    sffe *p = sffe_alloc();
    p->resolve = sffe_resolve_p;
    static cmplx z, c;
    GSL_SET_COMPLEX(&z, 0.25, -0.5);
    GSL_SET_COMPLEX(&c, 0.1, 0.2);
    sffe_regvar(&p, &z, "z");
    sffe_regvar(&p, &c, "c");
    if (sffe_parse(&p, expr) != 0) {
        printf("FAIL   will not parse: %s\n", expr);
        failures++;
    }
    return p;
}

/* What the ring says p<back> is, read the way the iteration loops read it:
 * out of the place the parser bound the name to. */
static cmplx tap(sffe *p, const char *name)
{
    sfvariable *v = sffe_var(p, name);
    if (!v) {
        printf("FAIL   %s was never asked for\n", name);
        failures++;
        static cmplx zero;
        GSL_SET_COMPLEX(&zero, 0, 0);
        return zero;
    }
    return *v->value;
}

static int same(cmplx a, cmplx b)
{
    return GSL_REAL(a) == GSL_REAL(b) && GSL_IMAG(a) == GSL_IMAG(b);
}

/* The shift the ring replaced, in full: an array of `depth` places, every one
 * of them moved along by one each pass, the newest written at the front. */
struct shift {
    cmplx place[32];
    unsigned int depth;
};

static void shift_reset(struct shift *s, unsigned int depth, number_t re,
                        number_t im)
{
    s->depth = depth;
    for (unsigned int i = 0; i < depth; i++)
        GSL_SET_COMPLEX(&s->place[i], re, im);
}

static void shift_push(struct shift *s, number_t re, number_t im)
{
    for (unsigned int i = s->depth - 1; i > 0; i--)
        s->place[i] = s->place[i - 1];
    GSL_SET_COMPLEX(&s->place[0], re, im);
}

int main(void)
{
    /* --- what a name means ------------------------------------------------ */
    check(sffe_pindex("p") == 1, "p alone is p1");
    check(sffe_pindex("p1") == 1, "p1 is the pass before this one");
    check(sffe_pindex("p6") == 6, "p6 is six passes back");
    check(sffe_pindex("p9999") == 9999, "p9999 is as far back as one may look");
    check(sffe_pindex("p10000") == 0, "and p10000 is not a variable at all");
    check(sffe_pindex("p01") == 0, "nor is p01, which is not how one writes 1");
    check(sffe_pindex("pi") == 0, "pi is the constant and not a p");
    check(sffe_pindex("parchment") == 0, "nor is a function whose name starts p");
    check(sffe_pindex("z") == 0 && sffe_pindex("") == 0,
          "and neither is anything else");

    /* --- only what the formula names is kept up to date -------------------- */
    {
        sffe *f = compile("z^2+p3+c");
        sffe_ptaps_build(f, NULL);
        check(sffe_ph.depth == 3, "a formula naming p3 keeps three passes");
        check(sffe_ph.tapcount == 1, "and fills one place per pass");
        sffe_free(&f);
    }
    {
        sffe *f = compile("z^2+c");
        sffe_ptaps_build(f, NULL);
        check(sffe_ph.depth == 0 && sffe_ph.tapcount == 0,
              "a formula naming none keeps nothing and fills nothing");
        sffe_free(&f);
    }
    {
        sffe *f = compile("p+p1+p2");
        sffe_ptaps_build(f, NULL);
        check(sffe_ph.tapcount == 2,
              "p and p1 are one place asked for under two names");
        sffe_free(&f);
    }
    {
        /* the initial formula may name them too, and asks for the same ring */
        sffe *f = compile("z^2+p2+c");
        sffe *i = compile("p5");
        sffe_ptaps_build(f, i);
        check(sffe_ph.depth == 5, "the deeper of the two formulas sets the depth");
        check(sffe_ph.tapcount == 2, "and both formulas get their places filled");
        sffe_free(&f);
        sffe_free(&i);
    }
    {
        /* a variable registered by an earlier parse of the same parser must
         * not go on being filled once the formula stops naming it */
        sffe *f = compile("z+p4");
        sffe_ptaps_build(f, NULL);
        check(sffe_ph.depth == 4, "a formula naming p4 keeps four passes");
        sffe_parse(&f, "z+p2");
        sffe_ptaps_build(f, NULL);
        check(sffe_ph.depth == 2 && sffe_ph.tapcount == 1,
              "and stops keeping them when the formula stops naming it");
        sffe_free(&f);
    }

    /* --- the ring against the shift ---------------------------------------
     *
     * Every depth from one to eight, started both ways the menu allows, and
     * driven over enough passes to wrap the ring several times.
     */
    for (unsigned int depth = 1; depth <= 8 && !failures; depth++) {
        for (int startzero = 0; startzero < 2; startzero++) {
            char expr[512];
            /* name every place, so that every one of them is compared */
            snprintf(expr, sizeof(expr), "p1");
            for (unsigned int k = 2; k <= depth; k++) {
                char more[16];
                snprintf(more, sizeof(more), "+p%u", k);
                strcat(expr, more);
            }
            sffe *f = compile(expr);
            sffe_ptaps_build(f, NULL);

            struct shift ref;
            number_t sre = startzero ? 0 : (number_t)-0.75;
            number_t sim = startzero ? 0 : (number_t)0.125;
            shift_reset(&ref, depth, sre, sim);
            sffe_phist_reset(sre, sim);

            int agree = 1;
            for (unsigned int pass = 0; pass < 40; pass++) {
                for (unsigned int k = 1; k <= depth; k++) {
                    char name[16];
                    snprintf(name, sizeof(name), "p%u", k);
                    if (!same(tap(f, name), ref.place[k - 1]))
                        agree = 0;
                }
                number_t re = (number_t)pass + 1;
                number_t im = -(number_t)pass;
                shift_push(&ref, re, im);
                sffe_phist_push(re, im);
            }
            char what[96];
            snprintf(what, sizeof(what),
                     "depth %u from %s: the ring says what the shift said",
                     depth, startzero ? "zero" : "the point");
            check(agree, what);
            sffe_free(&f);
        }
    }

    /* --- and at a depth no shift could have afforded ----------------------- */
    {
        sffe *f = compile("z^2+p9999*0.001+c");
        sffe_ptaps_build(f, NULL);
        check(sffe_ph.depth == 9999, "p9999 keeps nine thousand nine hundred "
                                   "and ninety-nine passes");
        check(sffe_ph.tapcount == 1, "and still fills one place per pass");

        sffe_phist_reset(7, 11);
        cmplx expect;
        GSL_SET_COMPLEX(&expect, 7, 11);
        check(same(tap(f, "p9999"), expect),
              "before any pass it stands at what the history started from");

        for (unsigned int pass = 0; pass < 9998; pass++)
            sffe_phist_push((number_t)pass, 0);
        check(same(tap(f, "p9999"), expect),
              "and goes on standing there until the ring has come round");

        sffe_phist_push(12345, 0);
        GSL_SET_COMPLEX(&expect, 0, 0);
        check(same(tap(f, "p9999"), expect),
              "then hands back the first value pushed, 9999 passes ago");
        sffe_free(&f);
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    else
        printf("\nok     the history reads back the way it always did\n");
    return failures != 0;
}
