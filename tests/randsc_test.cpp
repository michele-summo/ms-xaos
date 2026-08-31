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
#include <cstring>

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
    char what[96];
    sffe *bare = compile("randsc({7,0})");
    sffe *full = compile("randsc({7,0};{1,1};{1,1})");
    sffe *fade = compile("randsc({7,0};{1,1};{0.5,0.2})");
    sffe *sized = compile("randsc({7,0};{0.5,0.2};{1,1})");
    sffe *half = compile("randsc({7,0};{1,1};{0.5,1})");
    sffe *fifth = compile("randsc({7,0};{0.03125,1};{1,1})");
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

    /* The pass after it uses size * degradation, which is that size written
     * out and left alone. Compared at the same pass, because the pass is part
     * of the hash: two fields can only be told apart within one of them. */
    check(at(fade, 0.3, 0.7, 1) == at(sized, 0.3, 0.7, 1),
          "a later pass uses a smaller size");
    check(at(fade, 0.3, 0.7, 1) != at(bare, 0.3, 0.7, 1),
          "and not the size it started with");

    /* One component degraded and the other not. Both are raised to the power
     * together now, so the untouched one is multiplied by a one rather than
     * left alone; that is exact, and this says so. 0.5^5 is 0.03125 with no
     * rounding either, so the two spellings must agree to the bit. */
    check(at(half, 0.3, 0.7, 5) == at(fifth, 0.3, 0.7, 5),
          "a degradation of one on one component leaves it alone");

    /* With degradation 1+i nothing shrinks, so pass nine is still size 1+i. */
    check(at(full, 0.3, 0.7, 9) == at(bare, 0.3, 0.7, 9),
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

    /* A new field every pass, whatever the degradation is doing.
     *
     * This is the point of hashing the iteration and not merely scaling by it.
     * With the size alone carrying the pass, a degradation of one froze the
     * field -- the same point returned the same value for ever -- and one of
     * 0.99 shifted the grid by a percent, which is the same picture again. A
     * formula that calls one of these once per iteration is asking for a new
     * value each time.
     */
    {
        sffe *held = compile("randsc({13,0};{0.3,0.3};{1,1})");
        sffe *slow = compile("randsc({13,0};{0.3,0.3};{0.99,0.99})");
        sffe *heldq = compile("randscq({13,0};{0.3,0.3};{1,1})");
        sffe *heldp = compile("randscp({13,0};{0.3,0.3};{1,1})");
        sffe *heldh = compile("randsch({13,0};{0.3,0.3};{1,1})");
        sffe *heldt = compile("randsct({13,0};{0.3,0.3};{1,1})");
        if (!failures) {
            int mheld = 0, mslow = 0, mq = 0, mp = 0, mh = 0, mt = 0;
            for (unsigned int n = 0; n < 20; n++) {
                if (at(held, 0.3, 0.7, n) != at(held, 0.3, 0.7, n + 1))
                    mheld++;
                if (at(slow, 0.3, 0.7, n) != at(slow, 0.3, 0.7, n + 1))
                    mslow++;
                if (at(heldq, 0.3, 0.7, n) != at(heldq, 0.3, 0.7, n + 1))
                    mq++;
                if (at(heldp, 0.3, 0.7, n) != at(heldp, 0.3, 0.7, n + 1))
                    mp++;
                if (at(heldh, 0.3, 0.7, n) != at(heldh, 0.3, 0.7, n + 1))
                    mh++;
                if (at(heldt, 0.3, 0.7, n) != at(heldt, 0.3, 0.7, n + 1))
                    mt++;
            }
            check(mheld == 20, "randsc moves every pass at degradation 1");
            check(mslow == 20, "randsc moves every pass at degradation 0.99");
            check(mq == 20, "randscq moves every pass at degradation 1");
            check(mp == 20, "randscp moves every pass at degradation 1");
            check(mh == 20, "randsch moves every pass at degradation 1");
            check(mt == 20, "randsct moves every pass at degradation 1");

            /* Moves, and by a real amount: consecutive passes are unrelated
             * fields, not the same one nudged. */
            number_t spread = 0;
            for (unsigned int n = 0; n < 20; n++) {
                number_t a = at(held, 0.3, 0.7, n);
                number_t b = at(held, 0.3, 0.7, n + 1);
                spread += a > b ? a - b : b - a;
            }
            check(spread / 20 > (number_t)1 / 10,
                  "consecutive passes are unrelated, not nudged");
        }
    }

    /* The four mosaics share every argument and all of randsc_setup with
     * randsc; what differs is the last step, so that is what is checked. */
    {
        struct {
            const char *name;
            sffe *cells;
            sffe *zero;
        } mosaic[] = {
            {"randscq", compile("randscq({7,0})"),
             compile("randscq({7,0};{0,1})")},
            {"randscp", compile("randscp({7,0})"),
             compile("randscp({7,0};{1,1};{1,0})")},
            {"randsch", compile("randsch({7,0})"),
             compile("randsch({7,0};{0,1})")},
            {"randsct", compile("randsct({7,0})"),
             compile("randsct({7,0};{1,1};{1,0})")},
        };
        const int nmosaic = (int)(sizeof(mosaic) / sizeof(mosaic[0]));
        sffe *smooth = compile("randsc({7,0})");
        number_t step = (number_t)1 / 1000;

        for (int m = 0; !failures && m < nmosaic; m++) {
            sprintf(what, "%s declines to compute on a zero", mosaic[m].name);
            check(at(mosaic[m].zero, 0.3, 0.7, 0) == 0, what);

            /* Flat cells: within one cell the value does not move, where the
             * interpolated field does. */
            sprintf(what, "%s is flat across a cell", mosaic[m].name);
            check(at(mosaic[m].cells, 0.3, 0.7, 0) ==
                      at(mosaic[m].cells, (number_t)0.3 + step, 0.7, 0),
                  what);

            /* Four ways of cutting the plane, four different fields. */
            for (int n = m + 1; n < nmosaic; n++) {
                int alike = 0;
                for (int i = -8; i <= 8; i++)
                    for (int j = -8; j <= 8; j++) {
                        number_t x = (number_t)i / 3, y = (number_t)j / 3;
                        if (at(mosaic[m].cells, x, y, 0) ==
                            at(mosaic[n].cells, x, y, 0))
                            alike++;
                    }
                sprintf(what, "%s and %s are different mosaics",
                        mosaic[m].name, mosaic[n].name);
                check(alike == 0, what);
            }
        }

        if (!failures)
            check(at(smooth, 0.3, 0.7, 0) !=
                      at(smooth, (number_t)0.3 + step, 0.7, 0),
                  "randsc is not flat across a cell");

        /* One cell to the unit square, whatever its shape.
         *
         * randsc, randscq and randscp lay one cell over each unit square of
         * the degraded size and so have cells of unit area; a hexagon of
         * circumradius one and a triangle of side one do not, and are scaled
         * to match. Without that, changing one letter of a formula changed
         * the scale of the picture by a factor of six from end to end.
         *
         * Counting the flat regions over an area of 144 measures it: the
         * count is the area plus whatever the border cuts through, which is
         * of the order of the perimeter. Well away from a sixfold error.
         */
        for (int m = 0; !failures && m < nmosaic; m++) {
            number_t seen[400];
            int nseen = 0;
            for (int i = 0; i < 300 && nseen < 400; i++)
                for (int j = 0; j < 300 && nseen < 400; j++) {
                    number_t val =
                        at(mosaic[m].cells, (number_t)i / 25, (number_t)j / 25, 0);
                    int k;
                    for (k = 0; k < nseen; k++)
                        if (seen[k] == val)
                            break;
                    if (k == nseen)
                        seen[nseen++] = val;
                }
            sprintf(what, "%s gives one cell per unit square (%d over 144)",
                    mosaic[m].name, nseen);
            check(nseen >= 144 && nseen <= 260, what);
        }
    }

    /* Past the resolution of the grid.
     *
     * A cell is the position divided by the size, and it has to end up in an
     * integer, so the size cannot usefully go below about the position over
     * what an integer holds. A degradation of a half reaches that in some
     * sixty passes. What the field does there is not much of a choice -- one
     * flat value over the whole plane, there being no cell structure left to
     * resolve -- but two things about it matter, and both were wrong once.
     */
    {
        sffe *bh = compile("randsch({13,0};{0.3,0.3};{0.5,0.5})");
        sffe *bt = compile("randsct({14,0};{0.3,0.3};{0.5,0.5})");
        sffe *bq = compile("randscq({13,0};{0.3,0.3};{0.5,0.5})");
        if (!failures) {
            check(at(bh, 0.7, 0.4, 200) == at(bh, -1.3, 0.9, 200),
                  "past the resolution the field is flat");

            /* It must still change from pass to pass. It did not, once: the
             * hash was assigned after the range check and the caller read an
             * uninitialised one, which held the same value for every pass and
             * every seed. A formula then sat on one number to the iteration
             * limit. */
            int moved = 0;
            for (unsigned int n = 200; n < 220; n++)
                if (at(bh, 0.7, 0.4, n) != at(bh, 0.7, 0.4, n + 1))
                    moved++;
            check(moved == 20, "past the resolution it still moves each pass");

            /* And the functions must still differ from one another, or a
             * formula subtracting one from another reaches exactly zero and
             * iterates to the limit for nothing. */
            check(at(bh, 0.7, 0.4, 200) != at(bt, 0.7, 0.4, 200) &&
                      at(bh, 0.7, 0.4, 200) != at(bq, 0.7, 0.4, 200) &&
                      at(bt, 0.7, 0.4, 200) != at(bq, 0.7, 0.4, 200),
                  "past the resolution they still differ from each other");
        }
    }

    /* The values themselves, over a grid.
     *
     * floorl and roundl were replaced by integer conversions, which is where
     * most of the cost of one call was; the substitution is only sound if not
     * one value moved. The behaviour above says what the functions are for;
     * this says they are still computing the same numbers, which is what
     * keeps a saved picture drawing the same way. Per precision, since the
     * mosaics are step functions and the two builds are entitled to disagree
     * along a cell edge.
     */
    {
        static const struct {
            const char *name;
            unsigned long long expected[2]; /* 64 bits, 113 bits */
        } golden[] = {
            {"randsc", {0xe44c35d1b7e00368ULL, 0x599f28f8e0bb8b39ULL}},
            {"randscq", {0xb9b30f3fbedb4800ULL, 0xb9b30f3fbedb4800ULL}},
            {"randscp", {0x6d540398e03f4800ULL, 0x6d540398e03f4800ULL}},
            {"randsch", {0x174f609c68321800ULL, 0x174f609c68321800ULL}},
            {"randsct", {0xa82a262574eb5000ULL, 0xa82a262574eb5000ULL}},
        };
        const int which = NUMBER_MANTISSA_BITS == 113 ? 1 : 0;
        for (int g = 0; g < 5; g++) {
            char expr[64];
            sprintf(expr, "%s({7,0};{0.4,0.6};{0.9,0.95})", golden[g].name);
            sffe *f = compile(expr);
            if (!f)
                break;
            unsigned long long sum = 0;
            for (int i = -20; i < 20; i++)
                for (int j = -20; j < 20; j++)
                    for (unsigned int n = 0; n < 3; n++) {
                        number_t val = at(f, (number_t)i / 7, (number_t)j / 9, n);
                        sum = sum * 0x100000001B3ULL +
                              (unsigned long long)(val * 18446744073709551616.0);
                    }
            sprintf(what, "%s draws what it drew (%#llx)", golden[g].name, sum);
            check(sum == golden[g].expected[which], what);
        }
    }

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures != 0;
}
