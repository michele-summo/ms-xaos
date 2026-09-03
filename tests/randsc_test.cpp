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
