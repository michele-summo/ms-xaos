#ifndef VIEWTOOLS_H
#define VIEWTOOLS_H

/* The View dialog's tools, in viewtools.cpp: turning the angle by slices of a
 * full turn, and centering the view on an axis of the picture as it is turned.
 * Neither is a parameter of the view: they only fill in the dialog's angle and
 * center, which OK then applies. */

#include "config.h"

/* The angle `steps` slices on from `from`, a full turn being cut into `slices`
 * equal ones: steps * 360 / slices degrees, the other way when steps is
 * negative. Counted from where the turning started rather than added a slice
 * at a time, so that as many steps as there are slices come round to a whole
 * turn exactly. `from` itself when there are no slices to turn by. */
number_t view_turn(number_t from, int steps, int slices);

/* The axes a center goes on, as they lie in the picture: R across it and I up
 * it, which are the real and imaginary axes when the angle is nought. */
#define VIEW_AXIS_R 0
#define VIEW_AXIS_I 1

/* Moves the center (*cr, *ci) of a picture turned by `angle` degrees straight
 * across the picture (VIEW_AXIS_R) or straight up it (VIEW_AXIS_I) until the
 * origin is in the middle that way -- until that part of the center, measured
 * along the picture's own axes, is nought. At angle nought this is the real or
 * imaginary part becoming nought and the other left as it was; R and then I
 * bring the origin to the middle. A center already there is left alone. */
void view_centeronaxis(number_t angle, int axis, number_t *cr, number_t *ci);

#endif
