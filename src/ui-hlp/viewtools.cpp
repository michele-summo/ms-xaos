#include "config.h"
#include "number_math.h"
#include "viewtools.h"

number_t view_turn(number_t from, int steps, int slices)
{
    if (slices <= 0)
        return from;
    /* steps * 360 is whole and exact, and dividing it rounds once: n steps of
     * n slices are 360 to the last digit, where adding 360 / n n times would
     * leave the error of every addition in. */
    return from + (number_t)steps * 360 / slices;
}

void view_centeronaxis(number_t angle, int axis, number_t *cr, number_t *ci)
{
    /* The turn the picture is drawn at, as precalculate_rotation in
     * fractal.cpp works it out: the middle the center is moved to is then
     * the middle of the picture as it is drawn. */
    number_t sine, cosine;
    nsincosdeg(angle, &sine, &cosine);

    /* The center along the picture's axes, as my_rotate turns it. */
    number_t x = *cr * cosine - *ci * sine;
    number_t y = *cr * sine + *ci * cosine;

    /* A part that is nought comes out of the turning as the rounding of the
     * other, a few units in its last place: the center R left it in, read
     * back for I. Taken for nought, so that R and then I give the origin
     * exactly, and pressing either twice changes nothing. */
    const number_t slack = 4 * nldexp((number_t)1, 1 - NUMBER_MANTISSA_BITS) *
                           (nfabs(*cr) + nfabs(*ci));
    if (nfabs(x) <= slack)
        x = 0;
    if (nfabs(y) <= slack)
        y = 0;

    if (axis == VIEW_AXIS_R) {
        if (x == 0)
            return;
        x = 0;
    } else {
        if (y == 0)
            return;
        y = 0;
    }

    /* And back, as rotateback does. */
    *cr = x * cosine + y * sine;
    *ci = -x * sine + y * cosine;
    /* Nought times a negative cosine is a negative nought, which the dialog
     * would print as -0. */
    if (*cr == 0)
        *cr = 0;
    if (*ci == 0)
        *ci = 0;
}
