/* The View dialog's tools in viewtools.cpp.
 *
 * Turning by slices has to come round to a whole turn exactly, however many
 * slices, or the angle field fills with the error of every step. Centering on
 * an axis has to leave the origin in the middle of the picture as the engine
 * turns it, keep the part it does not move exactly as it was where the angle
 * allows, and give nought where nought is meant -- not 1e-17, nor -0.
 */

#include <cmath>
#include <cstdio>

#include "config.h"
#include "number_math.h"
#include "viewtools.h"

static int failures = 0;

static void check(int ok, const char *what)
{
    printf("%-6s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        failures++;
}

/* Nought, and not the negative one. */
static int nought(number_t x) { return x == 0 && ncopysign(1, x) > 0; }

/* The center along the picture's axes, worked out the way the engine does:
 * precalculate_rotation and my_rotate. */
static void as_drawn(number_t angle, number_t cr, number_t ci, number_t *x,
                     number_t *y)
{
    number_t s, c;
    nsincosdeg(angle, &s, &c);
    *x = cr * c - ci * s;
    *y = cr * s + ci * c;
}

static number_t size(number_t cr, number_t ci) { return nfabs(cr) + nfabs(ci); }

/* A few units in the last place of a number as large as `of`. */
static int tiny(number_t x, number_t of)
{
    return nfabs(x) <= 8 * nldexp((number_t)1, 1 - NUMBER_MANTISSA_BITS) * of;
}

int main()
{
    /* --- turning by slices --------------------------------------------- */
    check(view_turn(30, 1, 0) == 30 && view_turn(30, -3, 0) == 30,
          "no slices, no turning");
    check(view_turn(0, 1, 1) == 360 && view_turn(0, -1, 1) == -360,
          "one slice is a whole turn either way");
    int whole = 1;
    for (int n = 1; n <= 360; n++)
        whole &= view_turn(0, n, n) == 360 && view_turn(0, -n, n) == -360 &&
                 view_turn(10, n, n) == 370;
    check(whole, "n steps of n slices are a whole turn exactly, n up to 360");
    check(view_turn(0, 1, 6) == 60 && view_turn(0, 4, 6) == 240 &&
              view_turn(15, -1, 4) == -75,
          "six slices are 60 degrees, four 90");
    number_t seventh = 360;
    seventh /= 7;
    check(view_turn(0, 1, 7) == seventh, "a seventh of a turn, rounded once");

    /* --- the turn itself ---------------------------------------------------- */
    /* A turn and the turn back have to leave a point where it was, to the
     * precision of the build. The sine and cosine were a double pair in the
     * quad build, a turn only to sixteen digits, and a deep view turned by
     * them showed another place. */
    int unit = 1, there = 1;
    for (int k = -720; k <= 720; k += 7) {
        const number_t a = k + (number_t)1 / 3;
        number_t s, c;
        nsincosdeg(a, &s, &c);
        unit &= tiny(s * s + c * c - 1, 1);
        const number_t x0 = -0.743643887037158704752191506114774,
                       y0 = 0.131825904205311970493132056385139;
        const number_t x1 = x0 * c - y0 * s, y1 = x0 * s + y0 * c; /* my_rotate */
        const number_t x2 = x1 * c + y1 * s, y2 = -x1 * s + y1 * c; /* rotateback */
        there &= tiny(x2 - x0, 1) && tiny(y2 - y0, 1);
    }
    check(unit, "the sine and cosine of every angle make a turn to the last digits");
    check(there, "turned and turned back, a point is where it was");

    /* --- centering, square to the axes ------------------------------------- */
    number_t cr, ci;
    cr = -1.25, ci = 0.55;
    view_centeronaxis(0, VIEW_AXIS_R, &cr, &ci);
    check(nought(cr) && ci == (number_t)0.55, "angle 0: R puts the real part to nought");
    cr = -1.25, ci = 0.55;
    view_centeronaxis(0, VIEW_AXIS_I, &cr, &ci);
    check(cr == (number_t)-1.25 && nought(ci),
          "angle 0: I puts the imaginary part to nought");

    cr = -1.25, ci = 0.55;
    view_centeronaxis(90, VIEW_AXIS_R, &cr, &ci);
    check(cr == (number_t)-1.25 && nought(ci),
          "angle 90: R, across the picture, is the imaginary part");
    cr = -1.25, ci = 0.55;
    view_centeronaxis(90, VIEW_AXIS_I, &cr, &ci);
    check(nought(cr) && ci == (number_t)0.55,
          "angle 90: I, up the picture, is the real part");

    int square = 1;
    const number_t quarters[] = {-450, -270, -180, -90, 180, 270, 360, 450, 720};
    for (number_t a : quarters) {
        number_t r1 = 0.5, i1 = 0.3, r2 = 0.5, i2 = 0.3;
        view_centeronaxis(a, VIEW_AXIS_R, &r1, &i1);
        view_centeronaxis(a, VIEW_AXIS_I, &r2, &i2);
        const int odd = (int)nfmod(nfabs(a) / 90, 2);
        if (odd)
            square &= r1 == (number_t)0.5 && nought(i1) && nought(r2) &&
                      i2 == (number_t)0.3;
        else
            square &= nought(r1) && i1 == (number_t)0.3 &&
                      r2 == (number_t)0.5 && nought(i2);
    }
    check(square, "every quarter turn, either way round and past a whole "
                  "turn: one part nought, the other untouched, no -0");

    /* --- centering, turned ------------------------------------------------- */
    int turned = 1, both = 1, again = 1;
    const number_t angles[] = {30, 45, 60, 120, -30, 51.4285714285714285714,
                               1e-9, 389.5};
    for (number_t a : angles) {
        for (int axis = VIEW_AXIS_R; axis <= VIEW_AXIS_I; axis++) {
            const number_t r0 = -1.25, i0 = 0.55;
            number_t x0, y0, x, y;
            as_drawn(a, r0, i0, &x0, &y0);
            cr = r0, ci = i0;
            view_centeronaxis(a, axis, &cr, &ci);
            as_drawn(a, cr, ci, &x, &y);
            /* the origin in the middle that way, the other way unmoved */
            if (axis == VIEW_AXIS_R)
                turned &= tiny(x, size(r0, i0)) && tiny(y - y0, size(r0, i0));
            else
                turned &= tiny(y, size(r0, i0)) && tiny(x - x0, size(r0, i0));

            const number_t r1 = cr, i1 = ci;
            view_centeronaxis(a, axis, &cr, &ci);
            again &= cr == r1 && ci == i1;

            view_centeronaxis(a, axis == VIEW_AXIS_R ? VIEW_AXIS_I : VIEW_AXIS_R,
                              &cr, &ci);
            both &= nought(cr) && nought(ci);
        }
    }
    check(turned, "turned: the origin in the middle across or up the picture "
                  "as the engine draws it, the other way unmoved");
    check(again, "pressing the same axis twice changes nothing");
    check(both, "R then I, or I then R, is the origin exactly");

    cr = 0, ci = 0.75;
    view_centeronaxis(0, VIEW_AXIS_R, &cr, &ci);
    check(nought(cr) && ci == (number_t)0.75, "a center already on the axis stays");
    cr = 0, ci = 0;
    view_centeronaxis(30, VIEW_AXIS_I, &cr, &ci);
    check(nought(cr) && nought(ci), "the origin stays the origin");

    /* At a depth where the picture is far narrower than the rounding of the
     * turn: the center goes where it must to the precision of the build. */
    cr = -0.743643887037158704752191506114774, ci = 0.131825904205311970493132056385139;
    {
        number_t x, y;
        view_centeronaxis(30, VIEW_AXIS_R, &cr, &ci);
        as_drawn(30, cr, ci, &x, &y);
        check(tiny(x, 1), "deep in: R still leaves the origin in the middle");
    }

    printf("%d failures\n", failures);
    return failures != 0;
}
