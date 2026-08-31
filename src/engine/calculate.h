
#ifdef USE_SFFE
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#endif

static pixel32_t inline calculate(number_t x, number_t y, int periodicity);
static pixel32_t inline calculate(number_t x, number_t y, int periodicity)
{
    pixel32_t i;

    rotateback(cfractalc, x, y);
    if (cfractalc.plane) {
        recalculate(cfractalc.plane, &x, &y);
    }
    STAT(ncalculated2++);
#ifdef USE_SFFE
    /* What randsc and randscq hash. Here and not in the INIT macro, because
     * there the pixel is c in mandelbrot mode and z in julia mode -- binding
     * to either makes the noise change with the mode, or stand still. At this
     * point x and y are the pixel in both, already turned by the rotation and
     * mapped through the plane, which is the position the user sees. */
    cmplxset(sffe_position, x, y);

#endif
    if (cfractalc.mandelbrot) {
        if (cformula.flags & STARTZERO)
            i = cfractalc.calculate[periodicity](cfractalc.bre, cfractalc.bim,
                                                 x, y);
        else
            i = cfractalc.calculate[periodicity](x + cfractalc.bre,
                                                 y + cfractalc.bim, x, y);
    } else
        i = cfractalc.calculate[periodicity](x, y, cfractalc.pre,
                                             cfractalc.pim);
    return (i);
}
