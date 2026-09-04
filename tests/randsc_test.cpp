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
#include <thread>
#include <vector>

#include "config.h"
#include "number_math.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "misc-f.h"
#include <cmath>

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


/* --- driving a figure the way the iteration loops drive one ---------------
 *
 * A figure hands back where the point stands next, so it is fed back in, and
 * the pass it leaves the bailout shape on is what the picture is drawn from.
 * The shape is written out here rather than asked of the engine, so that this
 * says what the shape is and not merely that two copies of one mistake agree:
 * a bailout polygon stands its sides the square root of the bailout from the
 * centre -- that is its apothem -- and four is the bailout throughout.
 */
#define NEVER 0xffffffffu

static int fig_inside(int sides, number_t a, number_t x, number_t y)
{
    if (sides == 4) /* the square bailout tests the components apart */
        return x * x < a * a && y * y < a * a;
    number_t turn = sides == 3 ? -(number_t)M_PI / 2 : 0;
    for (int k = 0; k < sides; k++) {
        number_t t = turn + (number_t)k * 2 * (number_t)M_PI / sides;
        if (x * ncos(t) + y * nsin(t) >= a)
            return 0;
    }
    return 1;
}

/* the apothem is the square root of the bailout, and four is the bailout
 * throughout except where a figure of another radius is being compared */
static unsigned int leavesa(sffe *parser, int sides, number_t a, number_t x,
                            number_t y)
{
    if (!fig_inside(sides, a, x, y))
        return 0; /* never even started */
    /* the position is the pixel and stays where it is; z is what travels */
    GSL_SET_COMPLEX(&sffe_position, x, y);
    for (unsigned int n = 0; n < 64; n++) {
        GSL_SET_COMPLEX(&sffe_z, x, y);
        sffe_iteration = n;
        cmplx next = sffe_eval(parser);
        x = GSL_REAL(next);
        y = GSL_IMAG(next);
        if (!fig_inside(sides, a, x, y))
            return n + 1;
    }
    return NEVER;
}

/* Whether a figure puts every point it lets go in a place of its own.
 *
 * What it hands back on the pass a point goes is what every outside colouring
 * mode but the iteration count reads. A number standing in for "gone" would be
 * the same number for the whole figure and all of those modes would draw one
 * flat tone -- so what goes back has to be a point of the plane that moves as
 * the point moves. */
static int lands_apart(sffe *parser, int sides, number_t span)
{
    number_t sx[64], sy[64];
    int n = 0, gone = 0;
    for (int j = 0; j < 8; j++)
        for (int i = 0; i < 8; i++) {
            number_t x = -span + 2 * span * (number_t)i / 7;
            number_t y = -span + 2 * span * (number_t)j / 7;
            if (!fig_inside(sides, 2, x, y))
                continue;
            GSL_SET_COMPLEX(&sffe_position, x, y);
            for (unsigned int t = 0; t < 32; t++) {
                GSL_SET_COMPLEX(&sffe_z, x, y);
                sffe_iteration = t;
                cmplx next = sffe_eval(parser);
                x = GSL_REAL(next);
                y = GSL_IMAG(next);
                if (fig_inside(sides, 2, x, y))
                    continue;
                gone += 1;
                for (int q = 0; q < n; q++)
                    if (sx[q] == x && sy[q] == y)
                        return 0; /* two points, one place */
                if (n < 64) {
                    sx[n] = x;
                    sy[n] = y;
                    n += 1;
                }
                break;
            }
        }
    return gone > 8;
}

static unsigned int leaves(sffe *parser, int sides, number_t x, number_t y)
{
    return leavesa(parser, sides, 2, x, y);
}

int main(void)
{
    char what[96];
    sffe *bare = compile("randsc({7,0})");
    sffe *full = compile("randsc({7,0};{1,1};{1,1})");
    sffe *asdefault = compile("randsc({7,0};{1,1};{0.5,0.5})");
    sffe *fade = compile("randsc({7,0};{1,1};{0.5,0.2})");
    sffe *sized = compile("randsc({7,0};{0.5,0.2};{1,1})");
    sffe *half = compile("randsc({7,0};{1,1};{0.5,1})");
    sffe *fifth = compile("randsc({7,0};{0.03125,1};{1,1})");
    sffe *zsize = compile("randsc({7,0};{0,1})");
    sffe *zfade = compile("randsc({7,0};{1,1};{1,0})");
    if (failures)
        return 1;

    check(at(bare, 0.3, 0.7, 0) == at(full, 0.3, 0.7, 0),
          "the defaults are size 1+i and degradation a half each way");

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

    /* That a degradation of one leaves the size alone was checked here by
     * comparing it against a call that said nothing, which worked only while
     * saying nothing meant one. It says a half now, and the property is kept
     * by "a degradation of one on one component leaves it alone" above --
     * that one exercises the same multiplication by one, and against a
     * component that does move rather than against itself.
     *
     * What is worth checking here is the default: a call that says nothing
     * halves the cells each pass, which is what the two written out say. */
    check(at(bare, 0.3, 0.7, 4) == at(asdefault, 0.3, 0.7, 4),
          "saying nothing is size 1+i halving each pass");

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

    /* The degradation is carried from pass to pass rather than worked out
     * again, so the one thing that has to hold is that the route taken makes
     * no difference: a pixel reaching pass fifty in fifty steps and one
     * reaching it in a single jump must get the same number, or the picture
     * would depend on the order the pixels were computed in and no two
     * redraws would agree.
     *
     * The call site's running product is reset by asking for pass zero, which
     * is what the engine does at the start of every pixel. */
    {
        sffe *walk = compile("randsc({7,0};{0.4,0.6};{0.93,0.87})");
        if (!failures) {
            at(walk, 0.3, 0.7, 0);
            number_t jumped = at(walk, 0.3, 0.7, 50);

            at(walk, 0.3, 0.7, 0);
            for (unsigned int n = 1; n <= 50; n++)
                at(walk, 0.3, 0.7, n);
            number_t stepped = at(walk, 0.3, 0.7, 50);

            at(walk, 0.3, 0.7, 0);
            for (unsigned int n = 7; n <= 49; n += 7)
                at(walk, 0.3, 0.7, n);
            number_t bysevens = at(walk, 0.3, 0.7, 50);

            check(jumped == stepped && jumped == bysevens,
                  "the route to a pass makes no difference to it");

            /* And a pass earlier than the one carried starts again, which is
             * what tells one pixel from the next. */
            at(walk, 0.3, 0.7, 300);
            check(at(walk, 0.3, 0.7, 50) == jumped,
                  "going back to an earlier pass starts the count again");
        }
    }

    /* The carried degradation, computed by many threads at once.
     *
     * Each thread parses its own copy of a formula, so the running product
     * sits on a call site that belongs to one thread, and no two tread on
     * each other. That is the claim; this is the check. Pixels are dealt out
     * to the threads, each walked from pass zero as the engine walks them,
     * and every pixel's own checksum must come back the same whoever computed
     * it and in whatever order.
     *
     * Worth knowing about what is being relied on: the running product is not
     * a property of the pixel. It holds the degradation multiplied by itself
     * as many times as the pass says, and two pixels at the same pass want
     * the same number, so it does not matter which of them left it there. The
     * reset only serves for going backwards, which is what the start of a
     * pixel does.
     */
    {
        const int npixels = 600;
        const char *expr = "randsch({13,0};{0.3,0.3};{0.97,0.97})";
        std::vector<unsigned long long> one(npixels, 0), many(npixels, 0);

        for (int threads = 1; threads <= 4 && !failures; threads += 3) {
            std::vector<unsigned long long> &out = threads == 1 ? one : many;
            std::vector<std::thread> pool;
            for (int t = 0; t < threads; t++)
                pool.emplace_back([&, t, threads]() {
                    sffe *f = sffe_alloc();
                    if (!f || sffe_parse(&f, expr))
                        return;
                    /* dealt out unevenly and out of order, as a renderer
                     * hands out rows */
                    for (int k = t; k < npixels; k += threads) {
                        int i = (k * 373) % npixels;
                        unsigned int passes = 1 + (unsigned int)((i * 37) % 300);
                        unsigned long long sum = 1469598103934665603ULL;
                        for (unsigned int n = 0; n < passes; n++) {
                            number_t v = at(f, (number_t)(i % 71) / 23,
                                            (number_t)(i % 59) / 17, n);
                            sum ^= (unsigned long long)((double)v *
                                                        18446744073709551616.0);
                            sum *= 1099511628211ULL;
                        }
                        out[i] = sum;
                    }
                });
            for (auto &th : pool)
                th.join();
        }
        check(one == many && one[0] != 0,
              "four threads draw what one thread draws");
    }

    /* The kaleidoscope.
     *
     * The plane is folded into level wedges around the origin, so the field is
     * sampled from one of them and the picture repeats. The properties below
     * are what "folded" means, and the first of them is what says the two new
     * arguments cost nothing to a call that does not use them.
     */
    {
        sffe *bare = compile("randsc({7,0};{0.5,0.5};{1,1})");
        sffe *one = compile("randsc({7,0};{0.5,0.5};{1,1};{1,0};{0,0})");
        sffe *left = compile("randsc({7,0};{0.5,0.5};{1,1};{4,0};{0,0})");
        sffe *right = compile("randsc({7,0};{0.5,0.5};{1,1};{4,0};{1,0})");
        /* anything that is not one folds the way zero does, which is what
             * "and anything else" in the description has to mean */
        sffe *other = compile("randsc({7,0};{0.5,0.5};{1,1};{4,0};{7,0})");
        if (!failures) {
            /* A level of one is what a call with three arguments already
             * does, to the bit. */
            int same = 1;
            for (int i = -6; i <= 6; i++)
                for (int j = -6; j <= 6; j++)
                    if (at(bare, (number_t)i / 5, (number_t)j / 5, 3) !=
                        at(one, (number_t)i / 5, (number_t)j / 5, 3))
                        same = 0;
            check(same, "a kaleidoscope of one level changes nothing");

            /* Four wedges, so turning by a quarter of the plane lands on the
             * same value. Not to the bit: the fold goes through an angle and
             * back, which rounds. */
            const number_t quarter = N_PI / 2;
            const number_t radius = (number_t)17 / 10;
            const number_t angle = (number_t)4 / 10;
            number_t tol = (number_t)1 / 1000000000000ULL;
            int turns_ok = 1, mirror_ok = 1;
            for (int t = 1; t < 4; t++) {
                number_t a = angle + quarter * t;
                if (nfabs(at(left, radius * ncos(angle), radius * nsin(angle), 3) -
                          at(left, radius * ncos(a), radius * nsin(a), 3)) > tol)
                    turns_ok = 0;
            }
            check(turns_ok, "four wedges: a quarter turn lands on itself");

            /* Inside a wedge the far half mirrors the near one. */
            number_t m = quarter - angle;
            if (nfabs(at(left, radius * ncos(angle), radius * nsin(angle), 3) -
                      at(left, radius * ncos(m), radius * nsin(m), 3)) > tol)
                mirror_ok = 0;
            check(mirror_ok, "the wedge is a mirror of itself");

            /* The two modes are two different fields, and a mode that is
             * neither is the first of them. Not point by point -- the two
             * agree wherever neither happens to fold -- so over a grid. */
            int lr = 0, lo = 0;
            for (int i = -8; i <= 8; i++)
                for (int j = -8; j <= 8; j++) {
                    number_t x = (number_t)i / 5, y = (number_t)j / 5;
                    number_t l = at(left, x, y, 3), r = at(right, x, y, 3);
                    lr += l != r;
                    lo += l != at(other, x, y, 3);
                }
            check(lr > 40, "the two kaleidoscope modes are two fields");
            check(lo == 0, "a mode that is neither folds the way zero does");
        }
    }

    /* Orbit traps and stripe averaging.
     *
     * Both keep one number about a whole orbit, which is what the colouring
     * modes of the engine cannot do -- a calculation loop there hands back a
     * colour and only the last two values of z survive. Both hand back what
     * they were given until the last pass the limit allows, and what they
     * gathered on it.
     */
    {
        /* With a limit of one pass, pass zero is the last one, so a trap
         * reveals what it measured straight away: the distances below are the
         * shapes themselves, at the point 3+4i. */
        struct {
            const char *expr;
            double want;
            const char *what;
        } shapes[] = {
            {"trap({3,4};0)", 5, "shape 0 is the distance to the centre"},
            {"trap({3,4};1)", 4, "shape 1 is a horizontal line"},
            {"trap({3,4};2)", 3, "shape 2 is a vertical one"},
            {"trap({3,4};3)", 3, "shape 3 is the nearer of the two, a cross"},
            {"trap({3,4};4;{0,0};{2,0})", 3, "shape 4 is a ring"},
            {"trap({3,4};5;{0,0};{2,0})", 2, "shape 5 is a square"},
            {"trap({3,4};6;{0,0};{2,0})", 5, "shape 6 is a diamond"},
            {"trap({3,4};0;{3,0};{1,0})", 4, "the centre moves the shape"},
        };
        for (int i = 0; i < 8; i++) {
            sffe *f = compile(shapes[i].expr);
            if (!f)
                break;
            sffe_maxiter = 1;
            number_t got = at(f, 0.3, 0.7, 0);
            check(nfabs(got - (number_t)shapes[i].want) <
                      (number_t)1 / 1000000000,
                  shapes[i].what);
        }

        /* Handed straight back while the passes have not run out, and it is
         * the argument itself, to the bit. */
        sffe *pass = compile("trap({3,4};0)");
        if (!failures) {
            sffe_maxiter = 10;
            check(at(pass, 0.3, 0.7, 0) == 3 && at(pass, 0.3, 0.7, 5) == 3,
                  "a trap hands its argument back until the last pass");
        }

        /* The smallest of the whole orbit, not the last of it. ifiter walks
         * the argument round a cycle, so the orbit here is 3+4i, 0.5, 2, and
         * the smallest distance to the centre is a half. */
        sffe *walk = compile("trap(ifiter({3,4};{0.5,0};{2,0});0)");
        if (!failures) {
            sffe_maxiter = 6;
            for (unsigned int n = 0; n < 5; n++)
                at(walk, 0.3, 0.7, n);
            check(nfabs(at(walk, 0.3, 0.7, 5) - (number_t)0.5) <
                      (number_t)1 / 1000000000,
                  "a trap keeps the smallest of the orbit");

            /* And forgets it when the next pixel begins, or every pixel
             * after a near miss would report that near miss. The next pixel
             * here runs one pass only, so it must report the first point of
             * the orbit, five, and not the half the last pixel found. */
            sffe_maxiter = 1;
            check(nfabs(at(walk, 0.3, 0.7, 0) - 5) < (number_t)1 / 1000000000,
                  "and starts again with the next pixel");
        }

        /* The stripe average of a constant argument is the one sample, since
         * every pass measures the same angle. arg(1+i) is a quarter turn, so
         * a density of four asks for sin(pi) -- a half after the shift. */
        sffe *flat = compile("stripe({1,1};4)");
        if (!failures) {
            sffe_maxiter = 1;
            check(nfabs(at(flat, 0.3, 0.7, 0) - (number_t)0.5) <
                      (number_t)1 / 1000000,
                  "a stripe of one angle is that angle's sample");

            sffe_maxiter = 20;
            for (unsigned int n = 0; n < 19; n++)
                at(flat, 0.3, 0.7, n);
            check(nfabs(at(flat, 0.3, 0.7, 19) - (number_t)0.5) <
                      (number_t)1 / 1000000,
                  "and averaging it twenty times does not move it");
        }

        /* Whatever the orbit, an average of samples between zero and one is
         * between zero and one. */
        sffe *mixed = compile("stripe(ifiter({3,4};{0.5,-2};{2,1});5)");
        if (!failures) {
            sffe_maxiter = 12;
            for (unsigned int n = 0; n < 11; n++)
                at(mixed, 0.3, 0.7, n);
            number_t v = at(mixed, 0.3, 0.7, 11);
            check(v >= 0 && v <= 1, "a stripe average stays between 0 and 1");
        }
        sffe_maxiter = 0;
    }

    /* poly: a polynomial in its first argument, the rest its coefficients
     * from the highest power down. */
    {
        sffe *quad = compile("poly({2,0};{3,0};{4,0};{5,0})");
        sffe *linear = compile("poly({2,0};{3,0};{4,0})");
        sffe *constant = compile("poly({2,0};{7,0})");
        sffe *empty = compile("poly({2,0})");
        sffe *cube = compile("poly({2,0};{1,0};{0,0};{0,0};{0,0})");
        if (!failures) {
            /* 3*4 + 4*2 + 5 = 25 */
            check(at(quad, 0.3, 0.7, 0) == 25,
                  "poly with three coefficients is a quadratic");
            check(at(linear, 0.3, 0.7, 0) == 10, "with two, a line");
            check(at(constant, 0.3, 0.7, 0) == 7, "with one, a constant");
            check(at(empty, 0.3, 0.7, 0) == 0,
                  "with none, a sum of no terms");
            check(at(cube, 0.3, 0.7, 0) == 8, "the first coefficient takes "
                                              "the highest power");

            /* The imaginary part goes through the same arithmetic, so a
             * complex argument is worth one check of its own: (1+i)^2 = 2i,
             * so poly(1+i; 1; 0; 0) is 2i. */
            sffe *comp = compile("poly({1,1};{1,0};{0,0};{0,0})");
            if (!failures) {
                GSL_SET_COMPLEX(&sffe_position, 0.3, 0.7);
                sffe_iteration = 0;
                sfNumber v = sffe_eval(comp);
                check(GSL_REAL(v) == 0 && GSL_IMAG(v) == 2,
                      "and the arithmetic is complex throughout");
            }
        }
    }


    /* --- the figures: fractals of their own --------------------------------
     *
     * sierpinskyt and sierpinskyc carry the point to its parent: every hole in
     * a gasket sits inside a bigger one a level up, and doubling away from the
     * nearest corner is the step between them. The topmost hole has no parent
     * inside the figure, so that step carries it out, and a hole n levels down
     * takes n steps to get there -- the pass a point leaves on is the level it
     * was cut away at, and the iteration count is the picture.
     *
     * The snowflake has no parent to walk to, its levels adding to its edge
     * rather than cutting into its middle, so it is drawn whole: a point in it
     * is handed back where it stands and never leaves, one outside is thrown
     * past any bailout at once.
     */
    {
        sffe *tri = compile("sierpinskyt()");
        sffe *tri4 = compile("sierpinskyt(4)");
        sffe *tri16 = compile("sierpinskyt(16)");
        sffe *carpet = compile("sierpinskyc()");
        sffe *carpet33 = compile("sierpinskyc(4;3)");
        sffe *carpet5 = compile("sierpinskyc(4;5)");
        sffe *carpet5e = compile("sierpinskyc( ;5)");
        sffe *flake = compile("snowflake()");
        sffe *flake4 = compile("snowflake(4)");
        if (!failures) {
            /* every argument has a default, so the bare call is a call */
            check(leaves(tri, 3, (number_t)1 / 3, (number_t)1 / 7) ==
                      leaves(tri4, 3, (number_t)1 / 3, (number_t)1 / 7),
                  "sierpinskyt defaults to a radius of four");
            check(leaves(carpet, 4, (number_t)1 / 3, (number_t)1 / 7) ==
                      leaves(carpet33, 4, (number_t)1 / 3, (number_t)1 / 7),
                  "sierpinskyc to four and to three squares");
            check(leaves(flake, 6, (number_t)1 / 3, (number_t)1 / 7) ==
                      leaves(flake4, 6, (number_t)1 / 3, (number_t)1 / 7),
                  "and snowflake to four");
            check(leaves(carpet5e, 4, (number_t)1 / 3, (number_t)1 / 7) ==
                      leaves(carpet5, 4, (number_t)1 / 3, (number_t)1 / 7),
                  "and the radius may be left empty and the squares given");

            /* The figure is inscribed in the shape a bailout of the same number
             * draws, and for a triangle the corners stand at twice the number:
             * a triangle's apothem is half its circumradius. */
            /* A corner is a fixed point of the map and lies on the shape it
             * fills, so exactly on it nothing starts and a hair inside it the
             * point creeps away for many passes before it leaves. Which of the
             * two the last bit of a coordinate lands on is not the figure
             * talking, so the check is a hair inside. */
            check(leaves(tri4, 3, 0, (number_t)3999 / 1000) > 8,
                  "the apex stands on the corner of a triangular bailout of 4");
            check(leaves(tri4, 3, 0, (number_t)401 / 100) == 0,
                  "and just past it the point never started");
            check(leaves(tri4, 3, 0, (number_t)-199 / 100) != 0 &&
                      leaves(tri4, 3, 0, (number_t)-201 / 100) == 0,
                  "and its base lies along a side of that bailout");

            /* the middle of the triangle is the hole with no parent, so it is
             * carried out on the very first step -- which is the picture the
             * biggest triangle of all draws */
            check(leaves(tri4, 3, 0, 0) == 1,
                  "the middle of the gasket is the hole that leaves first");
            /* and the three below it take one step onto that one and a second
             * out: the middle of the top sub-triangle sits at half the height */
            check(leaves(tri4, 3, 0, 2) == 2,
                  "and a hole one level down takes one step more");

            /* the same shape twice the size, read at twice the distance */
            int scaled = 1;
            for (int i = 1; i < 40; i++) {
                number_t x = (number_t)(i % 7) / 3 - 1;
                number_t y = (number_t)(i % 11) / 5 - 1;
                if (leavesa(tri16, 3, 4, 2 * x, 2 * y) !=
                    leavesa(tri4, 3, 2, x, y))
                    scaled = 0;
            }
            check(scaled,
                  "and four times the radius is twice the figure, as bailout");

            /* the carpet: the middle square goes first, a corner never goes */
            check(leaves(carpet33, 4, 0, 0) == 1,
                  "the middle of the carpet is what leaves first");
            check(leaves(carpet33, 4, (number_t)4 / 3, (number_t)4 / 3) == 2,
                  "and a middle one level down takes one step more");
            check(leaves(carpet33, 4, (number_t)-1999 / 1000,
                         (number_t)-1999 / 1000) > 6,
                  "and the corner survives cut after cut");
            check(leaves(carpet33, 4, 5, 0) == 0,
                  "outside the square the point never started");
            /* What the cut leaves, counted rather than described: a cell in
             * the middle block goes on the first pass, and one in the border
             * ring is carried onto the middle block and goes on the second.
             * So the picture is the one big square the block draws with the
             * ring's 4n-4 around it, whatever n is. */
            struct {
                const char *expr;
                int n;
            } cuts[3] = {{"sierpinskyc(4;3)", 3},
                         {"sierpinskyc(4;4)", 4},
                         {"sierpinskyc(4;5)", 5}};
            for (int c = 0; c < 3; c++) {
                sffe *f = compile(cuts[c].expr);
                if (failures)
                    break;
                int n = cuts[c].n, ring = 0, hole = 0;
                for (int j = 0; j < n; j++) {
                    for (int i = 0; i < n; i++) {
                        /* the middle of that cell, in the plane */
                        number_t cx = (2 * ((number_t)i + (number_t)1 / 2) / n
                                       - 1) * 2;
                        number_t cy = (2 * ((number_t)j + (number_t)1 / 2) / n
                                       - 1) * 2;
                        unsigned int t = leaves(f, 4, cx, cy);
                        if (t == 1)
                            hole += 1;
                        else if (t == 2)
                            ring += 1;
                    }
                }
                check(hole == (n - 2) * (n - 2) && ring == 4 * n - 4,
                      cuts[c].expr);
                sffe_free(&f);
            }

            /* two squares to a side has no ring to speak of, so the far
             * quarter goes instead, which is a gasket again */
            sffe *carpet2 = compile("sierpinskyc(4;2)");
            if (!failures) {
                check(leaves(carpet2, 4, (number_t)3 / 2, (number_t)3 / 2) == 1,
                      "cut in two, the quarter that goes is the far corner");
                check(leaves(carpet2, 4, (number_t)-1999 / 1000,
                             (number_t)-1999 / 1000) > 6,
                      "and the near one stays");
                sffe_free(&carpet2);
            }

            /* The snowflake is read from its middle out: a hexagon there, six
             * triangles of the hexagon's own side standing on its six sides,
             * twelve of a third that on their free edges, and so on. The
             * hexagon goes on the first pass and a piece k levels out on the
             * kth, as the other two figures count.
             *
             * At a radius of four the figure reaches 2.3094 -- two over cos
             * thirty degrees -- so the hexagon stands its sides half of that
             * out, 1.1547, and its corners at 1.3333. Its corners point along
             * the axes and its sides face the ways between. */
            check(leaves(flake4, 6, 0, (number_t)115 / 100) == 1 &&
                      leaves(flake4, 6, 0, (number_t)116 / 100) == 2,
                  "a side of the middle hexagon stands half way out");
            check(leaves(flake4, 6, (number_t)132 / 100, 0) == 1,
                  "and a corner of it at that over cos thirty");
            check(leaves(flake4, 6, (number_t)134 / 100, 0) == NEVER,
                  "where two of the six triangles meet at a point and the "
                  "ground reaches in between them");

            /* the six triangles, each with its point on a corner of the
             * hexagonal bailout, all one level out from the middle */
            int six = 1;
            for (int d = 0; d < 6; d++) {
                double t = 3.14159265358979323846 / 6 +
                           d * 3.14159265358979323846 / 3;
                number_t r = (number_t)2309 / 1000;
                if (leaves(flake4, 6, (number_t)(r * ncos((number_t)t)),
                           (number_t)(r * nsin((number_t)t))) != 2)
                    six = 0;
            }
            check(six, "six triangles stand on it, and their points reach the "
                       "corners of a hexagonal bailout of 4");
            check(leaves(flake4, 6, (number_t)4444 / 10000,
                         (number_t)17963 / 10000) == 3,
                  "and twelve more on their free edges, one level further");

            /* From the hexagon out a snowflake is six-fold symmetric -- three
             * of the six triangles are the corners of the triangle it grew
             * from and three are the first bumps on its edges, and nothing
             * tells them apart. Reading it from the middle is what makes that
             * visible; the pass a point leaves on has to show it. */
            int sixfold = 1;
            for (int i = 1; i < 30; i++) {
                number_t px = (number_t)(i % 7) / 2 - (number_t)3 / 2;
                number_t py = (number_t)(i % 11) / 3 - (number_t)3 / 2;
                /* the same point turned by sixty degrees */
                number_t qx = px / 2 - py * (number_t)8660254 / 10000000;
                number_t qy = px * (number_t)8660254 / 10000000 + py / 2;
                if (leaves(flake4, 6, px, py) != leaves(flake4, 6, qx, qy))
                    sixfold = 0;
            }
            check(sixfold, "and a sixth of a turn leaves the picture alone");

            check(leaves(flake4, 6, 1, (number_t)-15 / 10) == NEVER,
                  "the ground beside the snowflake never leaves");
            check(leaves(flake4, 6, 0, (number_t)-2313 / 1000) == 0,
                  "and past the points nothing ever started");

            /* The hexagon has no parent, so it is thrown out rather than
             * declared gone: two points of it land in two places, and the
             * outside colouring modes have something to read. The gasket has
             * always done this -- its topmost hole is doubled away from a
             * corner like any other and lands outside the triangle, which is
             * outside the bailout, carrying where it came from with it. */
            check(lands_apart(flake4, 6, (number_t)12 / 10),
                  "a snowflake puts every point it lets go in its own place");
            check(lands_apart(tri4, 3, (number_t)18 / 10),
                  "and so does a gasket, as it always did");
            /* the ground is handed back untouched, which is what leaves the
             * incolouring something to say about it */
            {
                GSL_SET_COMPLEX(&sffe_z, 1, (number_t)-15 / 10);
                cmplx same = sffe_eval(flake4);
                check(GSL_REAL(same) == 1 &&
                          GSL_IMAG(same) == (number_t)-15 / 10,
                      "and stands exactly where it stood");
            }
            /* A triangle on a free edge is its parent cut down by three and
             * turned to face the way that edge faces, so a step off it lands on
             * the parent and it goes one pass later. Walking out along one of
             * the twelve, from the parent's edge to its own point. */
            int child = 1;
            for (int i = 0; i < 6; i++) {
                /* out along the way the bump on the upper right edge faces */
                number_t g = (number_t)(i + 1) / 20;
                number_t bx = (number_t)3333 / 10000 + g * (number_t)866 / 1000;
                number_t by = (number_t)17321 / 10000 + g / 2;
                if (leaves(flake4, 6, bx, by) != 3)
                    child = 0;
            }
            check(child, "a child is its parent cut down by three and turned "
                         "to face the way its edge faces");

            /* the whole silhouette, in a hundred and twenty directions: just
             * outside the shape nothing ever started, and every corner of the
             * shape is reached */
            struct {
                sffe *f;
                int sides;
                double turn;
                const char *what;
            } shape[2] = {{tri4, 3, -3.14159265358979323846 / 2,
                           "the gasket fills a triangular bailout of four"},
                          {flake4, 6, 0,
                           "and the snowflake reaches a hexagonal one"}};
            for (int k = 0; k < 2; k++) {
                double step = 2 * 3.14159265358979323846 / shape[k].sides;
                int agree = 1;
                for (int d = 0; d < 120; d++) {
                    double t = d * 2 * 3.14159265358979323846 / 120;
                    double off = t - shape[k].turn;
                    off -= step * (double)nfloor((number_t)(off / step) +
                                                 (number_t)0.5);
                    double reach = 2 / (double)ncos((number_t)off);
                    double cs = (double)ncos((number_t)t);
                    double sn = (double)nsin((number_t)t);
                    if (leaves(shape[k].f, shape[k].sides,
                               (number_t)(reach * 1.02 * cs),
                               (number_t)(reach * 1.02 * sn)) != 0)
                        agree = 0;
                }
                check(agree, shape[k].what);
            }
            /* And each reaches its corners. A hair inside them rather than on
             * them: a corner is one point of the plane, and whether a cosine
             * lands on the inside or the outside of it is decided by the last
             * bit rather than by the figure. */
            int corners = 1;
            for (int k = 0; k < 3; k++) {
                double t = 3.14159265358979323846 / 2 +
                           k * 2 * 3.14159265358979323846 / 3;
                number_t r = (number_t)4 * (number_t)999 / 1000;
                if (leaves(tri4, 3, (number_t)(r * ncos((number_t)t)),
                           (number_t)(r * nsin((number_t)t))) <= 3)
                    corners = 0;
            }
            for (int k = 0; k < 6; k++) {
                double t = 3.14159265358979323846 / 6 +
                           k * 3.14159265358979323846 / 3;
                number_t r = (number_t)23090 / 10000;
                /* in the figure rather than in the ground beside it: all
                 * six points belong to the six triangles standing on the
                 * middle hexagon, so all six go on the second pass */
                unsigned int at = leaves(flake4, 6,
                                         (number_t)(r * ncos((number_t)t)),
                                         (number_t)(r * nsin((number_t)t)));
                if (at == 0 || at == NEVER)
                    corners = 0;
            }
            check(corners, "and both reach every corner of the shape they fill");

            /* nothing in them is drawn from a seed */
            check(leaves(tri4, 3, (number_t)1 / 3, (number_t)1 / 7) ==
                      leaves(tri4, 3, (number_t)1 / 3, (number_t)1 / 7),
                  "and the same point always leaves on the same pass");
        }
        sffe_free(&tri);
        sffe_free(&tri4);
        sffe_free(&tri16);
        sffe_free(&carpet);
        sffe_free(&carpet33);
        sffe_free(&carpet5);
        sffe_free(&carpet5e);
        sffe_free(&flake);
        sffe_free(&flake4);
    }

    /* Defaults, where a call that stops short of an argument, or leaves its
     * place empty, takes them.
     *
     * Each is checked against the same call with the value written out: the
     * two must agree to the bit, or the default is not the value it is said
     * to be. The empty places matter as much as the short calls: a function
     * that reads its arguments by counting what it was given cannot tell an
     * empty place from a value, and would take the nothing for a zero. */
    {
        struct {
            const char *shorthand;
            const char *written;
            const char *what;
        } defaults[] = {
            {"julian({0.4,0.7})", "julian({0.4,0.7};{1,0};{1,0})",
             "julian defaults to the first power and the first turn"},
            {"julian({0.4,0.7};{2,0})", "julian({0.4,0.7};{2,0};{1,0})",
             "and to the first turn when only the power is given"},
            {"inveps({0.4,0.7})", "inveps({0.4,0.7};{0.01,0.01})",
             "inveps softens by a hundredth each way"},
            {"ngon({0.4,0.7})", "ngon({0.4,0.7};{0,0};{3,0};{1,0})",
             "ngon defaults to a triangle about the origin"},
            {"ngon({0.4,0.7};{0.1,0.2})",
             "ngon({0.4,0.7};{0.1,0.2};{3,0};{1,0})",
             "and to a triangle when only the centre is given"},
            {"randsc({7,0};{1,1})", "randsc({7,0};{1,1};{0.5,0.5})",
             "randsc halves its cells each pass"},
            {"randscq({7,0};{1,1})", "randscq({7,0};{1,1};{0.5,0.5};{1,0})",
             "and so does the mosaic, unfolded"},
            {"trap({0.4,0.7})", "trap({0.4,0.7};{0,0};{0,0};{1,0})",
             "trap measures to the centre by default"},
            {"stripe({0.4,0.7})", "stripe({0.4,0.7};{4,0})",
             "stripe lays four to a turn"},

            {"julian({0.4,0.7}; ;{2,0})", "julian({0.4,0.7};{1,0};{2,0})",
             "an empty place takes the default julian declares"},
            {"inveps({0.4,0.7}; )", "inveps({0.4,0.7};{0.01,0.01})",
             "and so does one left empty at the end"},
            {"ngon({0.4,0.7}; ;{5,0})", "ngon({0.4,0.7};{0,0};{5,0};{1,0})",
             "ngon about the origin with only the sides given"},
            {"randsc({7,0}; ; )", "randsc({7,0};{1,1};{0.5,0.5})",
             "randsc fills in size and degradation alike"},
            {"randscq({7,0}; ;{0.5,0.5}; )",
             "randscq({7,0};{1,1};{0.5,0.5};{1,0})",
             "and the mosaic its size and its kaleidoscope level"},
            {"trap({0.4,0.7}; ; ;{2,0})", "trap({0.4,0.7};{0,0};{0,0};{2,0})",
             "trap keeps its shape and its centre"},
            {"stripe({0.4,0.7}; )", "stripe({0.4,0.7};{4,0})",
             "stripe keeps its four to a turn"},
        };
        for (int i = 0; i < (int)(sizeof(defaults) / sizeof(defaults[0])) &&
                        !failures;
             i++) {
            sffe *shorthand = compile(defaults[i].shorthand);
            sffe *written = compile(defaults[i].written);
            if (failures)
                break;
            /* over several passes, since a default may only show itself once
             * the degradation has had something to bite on */
            int same = 1;
            for (unsigned int n = 0; n < 4; n++)
                if (at(shorthand, 0.3, 0.7, n) != at(written, 0.3, 0.7, n))
                    same = 0;
            check(same, defaults[i].what);
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
            {"randsc", {0x82c777cee11cf0acULL, 0x256878588c8eb491ULL}},
            {"randscq", {0x7b41886762e6b000ULL, 0x7b41886762e6b000ULL}},
            {"randscp", {0x7f716aced74a4000ULL, 0x7f716aced74a4000ULL}},
            {"randsch", {0x46b88bf1d206f000ULL, 0x46b88bf1d206f000ULL}},
            {"randsct", {0x224342c97c8b0800ULL, 0x224342c97c8b0800ULL}},
        };
        const int which = NUMBER_MANTISSA_BITS == 113 ? 1 : 0;
        for (int g = 0; g < 5; g++) {
            char expr[64];
            sprintf(expr, "%s({7,0};{0.4,0.6};{0.9,0.95})", golden[g].name);
            sffe *f = compile(expr);
            if (!f)
                break;
            unsigned long long sum = 0;
            static const unsigned int passes[] = {0, 1, 2, 7, 33, 100};
            for (int i = -20; i < 20; i++)
                for (int j = -20; j < 20; j++)
                    for (int k = 0; k < 6; k++) {
                        number_t val = at(f, (number_t)i / 7, (number_t)j / 9,
                                          passes[k]);
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
