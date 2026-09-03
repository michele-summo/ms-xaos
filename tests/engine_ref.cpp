/* Regression net for the fractal iteration loops.
 *
 * The formula table exposes every loop as an iterationfunc, so the loops can be
 * driven directly over a grid without going through the renderer. Each formula
 * gets a checksum of the iteration counts it produces; any change to the
 * arithmetic shows up as a different checksum.
 *
 * Usage:  engref            dump one line per formula
 *         engref --list     names only
 */
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "config.h"
#include "filter.h"
#include "fractal.h"
#include "xthread.h"

/* --- what formulas.cpp needs from the rest of XaoS ------------------------ */
struct fractal_context cfractalc;
struct palette cpalette;
struct taskinfo definfo;
int ethreads = 0;
int nthreads = 1;

void pth_function(void (*fn)(void *, struct taskinfo *, int, int), void *data,
                  int range)
{
    fn(data, &definfo, 0, range);
}
void pth_synchronize(void) {}
const char *qt_gettext(const char *, const char *t) { return t; }
void init_julia(struct image *, number_t, number_t, number_t, number_t) {}

/* -------------------------------------------------------------------------- */

static unsigned long long checksum(iterationfunc fn)
{
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = 0; iy < 48; iy++) {
        for (int ix = 0; ix < 48; ix++) {
            number_t pre = (number_t)-2.2 + (number_t)ix * (number_t)0.0687;
            number_t pim = (number_t)-1.3 + (number_t)iy * (number_t)0.0541;
            unsigned int it = fn(0, 0, pre, pim);
            h ^= it;
            h *= 1099511628211ULL;
        }
    }
    return h;
}

int main(int argc, char **argv)
{
    int names_only = (argc > 1 && !strcmp(argv[1], "--list"));

    /* The loops end in color_output, which indexes the palette, so a small
     * deterministic one has to exist. Its contents are arbitrary: the point is
     * that the same input gives the same output before and after a change. */
    static unsigned int pix[257];
    for (int k = 0; k < 257; k++)
        pix[k] = (unsigned int)(k * 2654435761u);
    cpalette.pixels = pix;
    cpalette.size = 256;
    cpalette.maxentries = 256;
    cpalette.start = 0;
    cpalette.end = 255;
    cpalette.type = 0;

    cfractalc.maxiter = 200;
    cfractalc.bailout = 4;
    /* What make_fractalc gives a real context. Left at zero the colouring
     * speed multiplies every mode's value by nothing, so every pixel comes
     * out the same colour and the checksums stop telling one mode from
     * another -- which is what hid smooth colouring skipping the speed and
     * the shift entirely. */
    cfractalc.incolorspeed = 1.0f;
    cfractalc.outcolorspeed = 1.0f;
    cfractalc.newtonconvergence = 1E-6;
    cfractalc.periodicity_limit = 1e-8;
    cfractalc.periodicity = 1;
    cfractalc.range = 2;
    cfractalc.incoloringmode = 0;
    cfractalc.pre = 0;
    cfractalc.pim = 0;
    cfractalc.bre = 0;
    cfractalc.bim = 0;

    int n = 0;
    for (const struct formula *f = formulas; f->magic == FORMULAMAGIC; f++, n++) {
        if (names_only) {
            printf("%s\n", f->name[0]);
            continue;
        }
        printf("%-24s", f->name[0]);
        printf(" calc=%016llx", f->calculate ? checksum(f->calculate) : 0ULL);
        printf(" peri=%016llx",
               f->calculate_periodicity ? checksum(f->calculate_periodicity)
                                        : 0ULL);
        printf(" scalc=%016llx",
               f->smooth_calculate ? checksum(f->smooth_calculate) : 0ULL);
        printf(" speri=%016llx",
               f->smooth_calculate_periodicity
                   ? checksum(f->smooth_calculate_periodicity)
                   : 0ULL);
        /* The same smooth loop with the colouring shifted, which has to come
         * out differently from the column before it. Smooth colouring used to
         * return a pixel of its own and never meet the shift or the speed at
         * all, so both controls did nothing whenever it was chosen; were that
         * to come back, these two columns would agree. */
        cfractalc.outcolorshift = 40;
        printf(" sshift=%016llx\n",
               f->smooth_calculate ? checksum(f->smooth_calculate) : 0ULL);
        cfractalc.outcolorshift = 0;
    }
    if (!names_only)
        printf("formule: %d\n", n);
    return 0;
}
