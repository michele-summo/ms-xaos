/* How the outside's fbm marks are sized, and the inside's tessellation.
 *
 * Fixed, where it starts, the three fbm modes must draw what they drew before
 * there was a choice, to the bit: their checksums are pinned here, taken from
 * the code before the change, the same at both precisions. A factor of
 * nought is fixed whatever the choice. Away from it the frequency doubles
 * factor times every ten iterations -- exactly, powers of two scaling
 * exactly -- the smooth count is the whole one where it is whole, and
 * nothing but the three fbm modes outside ever moves.
 *
 * The tessellation is the angle of z turned by a tile of its own: a tile is
 * found for every point, the turn is the same over one tile, and the tilings
 * draw apart from one another.
 *
 * formulas.cpp is included rather than linked, since color_output,
 * incolor_output and fbm_at_pixel are its own and static.
 */
#include <cstdio>
#include <cstring>

#include "formulas.cpp"

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

static int failures = 0;
static char what[512];

static void check(int ok, const char *text)
{
    printf("%s   %s\n", ok ? "ok  " : "FAIL", text);
    if (!ok)
        failures++;
}

/* An outside mode over escape points, counts and places in the plane: the
 * place is what the noise is read at, so it moves with the point. */
static unsigned long long outside(OutColormodeClass::OutColormode m)
{
    cfractalc.coloringmode = m;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = -40; iy < 40; iy++)
        for (int ix = -40; ix < 40; ix++) {
            number_t zre = ((number_t)ix + (number_t)0.5) / 9;
            number_t zim = ((number_t)iy + (number_t)0.5) / 11;
            color_px = (number_t)-0.75 + (number_t)ix / 64;
            color_py = (number_t)0.1 + (number_t)iy / 64;
            unsigned int iter = 1 + (unsigned int)((ix + 40) * 7 + iy + 40) % 97;
            h ^= color_output(zre, zim, iter);
            h *= 1099511628211ULL;
        }
    return h;
}

/* And an inside one over orbits. */
static unsigned long long inside(int m)
{
    cfractalc.incoloringmode = m;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = -40; iy < 40; iy++)
        for (int ix = -40; ix < 40; ix++) {
            number_t zre = ((number_t)ix + (number_t)0.5) / 23;
            number_t zim = ((number_t)iy + (number_t)0.5) / 29;
            color_px = zre;
            color_py = zim;
            h ^= incolor_output(zre, zim, (number_t)0.3, (number_t)-0.4, 100);
            h *= 1099511628211ULL;
        }
    return h;
}

static void scale(int scaling, number_t factor)
{
    cfractalc.outfbmscaling = scaling;
    cfractalc.outfbmfactor = factor;
}

static const OutColormodeClass::OutColormode fbm_iter =
    OutColormodeType::ColOut_fbm_iter;
static const OutColormodeClass::OutColormode fbm_alone =
    OutColormodeType::ColOut_fbm;

int main(void)
{
    static unsigned int pix[257];
    for (int k = 0; k < 257; k++)
        pix[k] = k * 2654435761u;
    cpalette.pixels = pix;
    cpalette.size = 200;
    cpalette.maxentries = 256;
    cpalette.type = C256;
    cfractalc.maxiter = 1000;
    /* as make_fractalc starts a context */
    cfractalc.incolorspeed = cfractalc.outcolorspeed = 1;
    cfractalc.incolorfuninitial = cfractalc.incolorfuninitialspeed = 1;
    cfractalc.incolorfunspeed = cfractalc.incolorfun2speed = 1;
    cfractalc.outcolorfuninitial = cfractalc.outcolorfuninitialspeed = 1;
    cfractalc.outcolorfunspeed = cfractalc.outcolorfun2speed = 1;
    cfractalc.incolorfrequency = 1;
    cfractalc.infbmintensity = cfractalc.outfbmintensity = 4;
    cfractalc.infbmfrequency = cfractalc.outfbmfrequency = 8;
    cfractalc.infbmoctaves = cfractalc.outfbmoctaves = 4;
    cfractalc.infbmroughness = cfractalc.outfbmroughness = (number_t)0.5;
    cfractalc.infbmseed = cfractalc.outfbmseed = 1;
    cfractalc.indecomprepeats = cfractalc.outdecomprepeats = 1;
    cfractalc.indecompspeed = cfractalc.outdecompspeed = 1;
    cfractalc.outfbmscaling = FBM_SCALE_FIXED;
    cfractalc.outfbmfactor = 1;
    cfractalc.intessellation = 3;
    color_change_prepare(&cfractalc);

    /* --- fixed, what there was ------------------------------------------ */
    static const unsigned long long was[3] = {
        0x48a6dc546e35c4ffULL, /* fbm + iter */
        0xfc38d2310e06af9dULL, /* fbm */
        0x4e6ca2eb201ff101ULL, /* fbm, inside */
    };
    const unsigned long long fixed[3] = {outside(fbm_iter), outside(fbm_alone),
                                         inside(INCOLORING_FBM)};
    static const char *const name[3] = {"fbm + iter", "fbm", "fbm inside"};
    for (int k = 0; k < 3; k++) {
        sprintf(what, "%s with the marks fixed draws as it did (%#llx)",
                name[k], fixed[k]);
        check(fixed[k] == was[k], what);
    }

    /* --- a factor of nought is fixed ------------------------------------ */
    for (int s = FBM_SCALE_ITER; s <= FBM_SCALE_SMOOTH; s++) {
        scale(s, 0);
        sprintf(what, "%s at a factor of nought is fixed",
                fbm_scalings[s]);
        check(outside(fbm_iter) == fixed[0] && outside(fbm_alone) == fixed[1],
              what);
    }

    /* --- the frequency doubles factor times every ten iterations --------- */
    {
        /* at a factor of ten, once an iteration: the motion at a count of n
         * is the motion at 2^n times the frequency, fixed, to the bit */
        int same = 1;
        for (int n = 0; n < 24 && same; n++)
            for (int k = 0; k < 50; k++) {
                color_px = (number_t)-0.75 + (number_t)k / 37;
                color_py = (number_t)0.1 + (number_t)k / 53;
                scale(FBM_SCALE_ITER, 10);
                number_t a = fbm_at_pixel(0, (number_t)n);
                scale(FBM_SCALE_FIXED, 1);
                number_t f = cfractalc.outfbmfrequency;
                cfractalc.outfbmfrequency = f * (number_t)((uint64_t)1 << n);
                number_t b = fbm_at_pixel(0);
                cfractalc.outfbmfrequency = f;
                if (a != b) {
                    printf("       count %d, point %d\n", n, k);
                    same = 0;
                    break;
                }
            }
        check(same, "at a factor of ten a count of n is the motion at 2^n "
                    "times the frequency");

        /* and at a factor of one, ten iterations to a doubling: where the
         * count is a whole number of tens it is one field alone, and between
         * them a blend held within nought to the intensity */
        int whole = 1, held = 1;
        for (int k = 0; k < 200; k++) {
            color_px = (number_t)-0.75 + (number_t)k / 91;
            color_py = (number_t)0.1 + (number_t)k / 67;
            scale(FBM_SCALE_ITER, 1);
            number_t a = fbm_at_pixel(0, 30);
            scale(FBM_SCALE_ITER, 10);
            whole &= a == fbm_at_pixel(0, 3);
            scale(FBM_SCALE_ITER, 1);
            for (int n = 31; n < 40; n++) {
                number_t v = fbm_at_pixel(0, (number_t)n);
                held &= v >= 0 && v <= cfractalc.outfbmintensity;
            }
        }
        check(whole, "at a factor of one, thirty iterations are three "
                     "doublings");
        check(held, "and between them the blend stays within the intensity");
        scale(FBM_SCALE_FIXED, 1);
    }

    /* --- the smooth count is the whole one where it is whole ------------- */
    {
        int same = 1;
        for (int n = 1; n < 60 && same; n++) {
            color_px = (number_t)n / 13;
            color_py = (number_t)n / 17;
            scale(FBM_SCALE_ITER, (number_t)1.5);
            cfractalc.coloringmode = fbm_iter;
            unsigned int a = color_output(0, 0, (unsigned int)n);
            scale(FBM_SCALE_SMOOTH, (number_t)1.5);
            color_smooth = (number_t)n * 256;
            color_smooth_held = 1;
            unsigned int b = color_output(0, 0, (unsigned int)n);
            same = a == b && !color_smooth_held;
        }
        check(same, "a smooth count that is whole sizes the marks as the "
                    "whole count does, and is taken down after");
        /* and without one handed on, the whole count stands in for it */
        scale(FBM_SCALE_ITER, (number_t)1.5);
        unsigned long long it = outside(fbm_iter);
        scale(FBM_SCALE_SMOOTH, (number_t)1.5);
        check(outside(fbm_iter) == it,
              "with no smooth count worked out the whole one is used");
        check(it != fixed[0], "and the marks move");
        scale(FBM_SCALE_FIXED, 1);
    }

    /* --- it moves the outside's fbm modes and nothing else --------------- */
    {
        int untouched = 1;
        unsigned long long before[OUTCOLORING], after[OUTCOLORING];
        for (int m = 0; m < OUTCOLORING; m++)
            before[m] = outside((OutColormodeClass::OutColormode)m);
        unsigned long long in[INCOLORING], in2[INCOLORING];
        for (int m = 0; m < INCOLORING; m++)
            in[m] = inside(m);
        scale(FBM_SCALE_ITER, 2);
        for (int m = 0; m < OUTCOLORING; m++) {
            after[m] = outside((OutColormodeClass::OutColormode)m);
            int fbm = m == OutColormodeType::ColOut_fbm_smooth ||
                      m == OutColormodeType::ColOut_fbm_iter ||
                      m == OutColormodeType::ColOut_fbm;
            /* fbm + smooth is drawn by the smooth loop, which this does not
             * drive; the other two move here */
            if (fbm && m != OutColormodeType::ColOut_fbm_smooth &&
                after[m] == before[m]) {
                printf("       %s did not move\n", outcolorname[m]);
                untouched = 0;
            }
            if (!fbm && after[m] != before[m]) {
                printf("       %s moved\n", outcolorname[m]);
                untouched = 0;
            }
        }
        for (int m = 0; m < INCOLORING; m++) {
            in2[m] = inside(m);
            if (in2[m] != in[m]) {
                printf("       %s inside moved\n", incolorname[m]);
                untouched = 0;
            }
        }
        check(untouched, "the sizing moves the outside's fbm modes and "
                         "nothing else, inside or out");
        scale(FBM_SCALE_FIXED, 1);
    }

    /* --- the tessellation ------------------------------------------------ */
    {
        int everywhere = 1, apart = 1;
        unsigned long long seen[RANDSCTILE_TILINGS + 1];
        for (int t = 1; t <= RANDSCTILE_TILINGS; t++) {
            uint64_t name;
            for (int k = 0; k < 400 && everywhere; k++) {
                number_t x = ((number_t)(k % 20) - 10 + (number_t)0.37) * 3;
                number_t y = ((number_t)(k / 20) - 10 + (number_t)0.61) * 3;
                if (!randsctile_name(t, x, y, &name)) {
                    printf("       tiling %d has no tile at (%g, %g)\n", t,
                           (double)x, (double)y);
                    everywhere = 0;
                }
            }
            cfractalc.intessellation = t;
            seen[t] = inside(INCOLORING_TESSELLATION);
            for (int u = 1; u < t; u++)
                if (seen[u] == seen[t]) {
                    printf("       tilings %d and %d draw the same\n", u, t);
                    apart = 0;
                }
        }
        check(everywhere, "every tiling has a tile for every point");
        check(apart, "and every one draws a picture of its own");

        /* a tile names itself and no neighbour: squares of unit area, the
         * middles of a row of them */
        uint64_t a, b, c;
        int one = 1, two = 1;
        for (int k = -20; k < 20; k++) {
            randsctile_name(1, (number_t)k + (number_t)0.5, (number_t)0.5, &a);
            randsctile_name(1, (number_t)k + (number_t)0.45, (number_t)0.55,
                            &b);
            randsctile_name(1, (number_t)k + (number_t)1.5, (number_t)0.5, &c);
            one &= a == b;
            two &= a != c;
        }
        check(one, "two points of one square are one tile");
        check(two, "and the next square is another");

        check(!randsctile_name(0, 0, 0, &a) && !randsctile_name(46, 0, 0, &a),
              "a tiling outside 1 to 45 has none");
        cfractalc.intessellation = 3;
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    return failures != 0;
}
