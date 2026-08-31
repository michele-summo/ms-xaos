/* randsc: coherent noise over the point being iterated.
 *
 * The properties checked here are the ones the function is for, and each of
 * them was a decision rather than an accident:
 *
 *  - the defaults are 1+i for both size and degradation, so a bare call must
 *    equal the fully written one;
 *  - degradation applies component by component and starts at the power zero,
 *    so the first pass is the size as given;
 *  - a zero in either component of either argument divides by zero once the
 *    degradation reaches it, so the function declines to compute instead;
 *  - and the whole point of coherent noise over a hash: near inputs give near
 *    outputs. That is what keeps the two precisions agreeing and what stops
 *    the picture shimmering when the zoom engine reuses a row.
 */

#include <cstdio>

#include "config.h"
#include "number_math.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "misc-f.h"

const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

static int failures = 0;

static void check(int ok, const char *what)
{
    printf("%-6s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        failures++;
}

static sffe *compile(const char *expression)
{
    sffe *parser = sffe_alloc();
    if (!parser || sffe_parse(&parser, expression)) {
        printf("FAIL   cannot parse %s\n", expression);
        failures++;
        return NULL;
    }
    return parser;
}

static number_t at(sffe *parser, number_t x, number_t y, unsigned int n)
{
    GSL_SET_COMPLEX(&sffe_position, x, y);
    sffe_iteration = n;
    return GSL_REAL(sffe_eval(parser));
}

int main(void)
{
    sffe *bare = compile("randsc({7,0})");
    sffe *full = compile("randsc({7,0};{1,1};{1,1})");
    sffe *fade = compile("randsc({7,0};{1,1};{0.5,0.2})");
    sffe *zsize = compile("randsc({7,0};{0,1})");
    sffe *zfade = compile("randsc({7,0};{1,1};{1,0})");
    if (failures)
        return 1;

    check(at(bare, 0.3, 0.7, 0) == at(full, 0.3, 0.7, 0),
          "the defaults are size 1+i and degradation 1+i");

    /* degradation to the power zero is one, so the first pass uses the size as
     * written -- which is the bare call with no degradation at all. */
    check(at(fade, 0.3, 0.7, 0) == at(bare, 0.3, 0.7, 0),
          "the first pass uses the size as given");
    check(at(fade, 0.3, 0.7, 1) != at(fade, 0.3, 0.7, 0),
          "a later pass uses a smaller size");

    /* With degradation 1+i nothing shrinks, so every pass is alike. */
    check(at(full, 0.3, 0.7, 0) == at(full, 0.3, 0.7, 9),
          "degradation 1+i leaves the size alone");

    check(at(zsize, 0.3, 0.7, 0) == 0, "a zero in size declines to compute");
    check(at(zfade, 0.3, 0.7, 0) == 0,
          "a zero in degradation declines to compute");

    /* Coherence. Over a blob of size 1 a step of 1/10000 may move the value by
     * a small fraction of its range, and never by the whole of it, which is
     * what a hash would do. */
    number_t a = at(bare, 0.3, 0.7, 0);
    number_t b = at(bare, (number_t)0.3 + (number_t)1 / 10000, 0.7, 0);
    number_t step = a > b ? a - b : b - a;
    check(step < (number_t)1 / 100, "near points give near values");

    /* And it is not merely flat: far apart, the values are unrelated. */
    number_t far = at(bare, 5.7, -3.1, 0);
    number_t spread = a > far ? a - far : far - a;
    check(spread > (number_t)1 / 100, "distant points give unrelated values");

    /* Same input, same answer -- no hidden state, unlike rand. */
    check(at(bare, 0.3, 0.7, 0) == a, "the same point always gives the same value");

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures != 0;
}
