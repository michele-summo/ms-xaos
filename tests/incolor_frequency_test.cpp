/* The frequency of the inside modes that repeat.
 *
 * At one, the six modes that repeated before it must give what they gave
 * before there was such a thing, to the bit, since every position saved
 * until now draws at one. Away from one it must do what it says -- make what
 * repeats repeat that much faster -- which for the modes that read z alone is
 * the same as drawing them at z times the frequency, exactly; the five written
 * since are held to that too, where they read z alone. And the modes that do
 * not repeat must never read it.
 *
 * The same for the speed function they bend the argument of every periodic
 * function by: ident, where it starts, changes nothing; another bends the
 * argument before the frequency multiplies it, k*f(a), in all eleven -- which
 * for sin(real)*sin(imag) is the same as drawing it at z bent -- and the
 * modes that do not repeat never read it.
 *
 * formulas.cpp is included rather than linked, since incolor_output is its own
 * and static.
 */
#include <cstdio>
#include <cstring>
#include <initializer_list>

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

/* The six that repeated before the frequency, by their numbers in
 * incolorname, and with the five written since, all that repeat. */
static const int periodic[6] = {5, 6, 7, 9, 18, 20};
static const int repeating[11] = {5,
                                  6,
                                  7,
                                  9,
                                  18,
                                  20,
                                  INCOLORING_COS_ABS,
                                  INCOLORING_HEXAGONS,
                                  INCOLORING_SPIRAL,
                                  INCOLORING_COS_ABS_ZC,
                                  INCOLORING_TESSELLATION};

/* Neither repeats nor is drawn by incolor_output's palette: true colour and
 * the fbm modes, whose noise is not a periodic function. */
static int set_apart(int m)
{
    return m == INCOLORING_TRUECOLOR ||
           (m >= INCOLORING_FBM && m <= INCOLORING_FBM_DECOMP);
}

/* What mode m gives at frequency f and periodic function g over a grid of
 * orbits, z scaled by s and then, if there is one, put through map. */
static unsigned long long sweep(int m, number_t f, number_t s, int g = 0,
                                number_t (*map)(number_t) = NULL)
{
    cfractalc.incoloringmode = m;
    cfractalc.incolorfrequency = f;
    cfractalc.incolorperiodicfun = g;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = -60; iy <= 60; iy++)
        for (int ix = -60; ix <= 60; ix++) {
            number_t zre = (number_t)ix / 23, zim = (number_t)iy / 29;
            number_t pre = (number_t)0.3 + zre / 7;
            number_t pim = (number_t)-0.4 + zim / 5;
            number_t wre = zre * s, wim = zim * s;
            if (map != NULL) {
                wre = map(wre);
                wim = map(wim);
            }
            unsigned int r = incolor_output(wre, wim, pre, pim, 100);
            h ^= r;
            h *= 1099511628211ULL;
        }
    return h;
}

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
    /* as make_fractalc starts a context: the speed and every part of the
     * speed change at one, which leaves the value alone */
    cfractalc.incolorspeed = 1;
    cfractalc.incolorfuninitial = cfractalc.incolorfuninitialspeed = 1;
    cfractalc.incolorfunspeed = cfractalc.incolorfun2speed = 1;
    cfractalc.indecomprepeats = cfractalc.outdecomprepeats = 1;
    cfractalc.indecompspeed = cfractalc.outdecompspeed = 1;
    cfractalc.intessellation = 3;
    color_change_prepare(&cfractalc);

    /* --- at one, what there was ---------------------------------------
     *
     * Checksums of the six over the grid, taken before the frequency was
     * added; the same at both precisions. */
    {
        static const unsigned long long was[6] = {
            0x769b49a3d3a78e88ULL, /* cos(mag) */
            0xc76df1b9d86ac38eULL, /* mag*cos(real^2) */
            0xa0c905c68fd924d6ULL, /* sin(real^2-imag^2) */
            0xc0ea6727c519c9f3ULL, /* squares */
            0x5b53e0d76e58c628ULL, /* sin(real)*sin(imag) */
            0xe90ffcb7543e56caULL, /* frac(mag) */
        };
        for (int k = 0; k < 6; k++) {
            sprintf(what, "%s at a frequency of one draws as it did",
                    incolorname[periodic[k]]);
            check(sweep(periodic[k], 1, 1) == was[k], what);
        }
    }

    /* --- away from one, what it says -----------------------------------
     *
     * The modes that read z alone: at frequency f, z, the picture they draw
     * at frequency one and z scaled by f -- exactly, since the frequencies
     * here are powers of two and the scaling is exact in binary. For
     * sin(real^2-imag^2) the argument is z squared, so a frequency of four is
     * z twice as far. cos(|z|) and hexagons read the modulus and the
     * projections of z, which scale with it, and tessellation its parts,
     * as squares does. */
    {
        struct {
            int mode;
            number_t f, s;
        } same[7] = {{7, 4, 2},
                     {9, 2, 2},
                     {18, 2, 2},
                     {20, 2, 2},
                     {INCOLORING_COS_ABS, 2, 2},
                     {INCOLORING_HEXAGONS, 2, 2},
                     {INCOLORING_TESSELLATION, 2, 2}};
        for (int k = 0; k < 7; k++) {
            sprintf(what,
                    "%s at a frequency of %g is the picture at one with z "
                    "%g times as far",
                    incolorname[same[k].mode], (double)same[k].f,
                    (double)same[k].s);
            check(sweep(same[k].mode, same[k].f, 1) ==
                      sweep(same[k].mode, 1, same[k].s),
                  what);
        }
        /* and the ones that read c as well, or the logarithm of the modulus:
         * something else than at one */
        for (int m : {5, 6, INCOLORING_SPIRAL, INCOLORING_COS_ABS_ZC}) {
            sprintf(what, "%s at a frequency of two is not what it is at one",
                    incolorname[m]);
            check(sweep(m, 2, 1) != sweep(m, 1, 1), what);
        }
    }

    /* --- the rest never read it ----------------------------------------- */
    {
        int untouched = 1;
        for (int m = 0; m < INCOLORING; m++) {
            int repeats = 0;
            for (int k = 0; k < 11; k++)
                if (repeating[k] == m)
                    repeats = 1;
            if (repeats || set_apart(m))
                continue;
            if (sweep(m, 3, 1) != sweep(m, 1, 1)) {
                printf("       %s moved\n", incolorname[m]);
                untouched = 0;
            }
        }
        check(untouched, "the modes that do not repeat draw the same at any "
                         "frequency");
    }

    /* --- the function they bend the argument by ------------------------- */
    {
        /* It reaches every mode that repeats -- squares and the bands of
         * cos(mag) as well as the sines and frac(mag). */
        int moved = 1;
        for (int k = 0; k < 11; k++)
            for (int g = 1; g < COLORFUN - 1; g++)
                if (sweep(repeating[k], 1, 1, g) ==
                    sweep(repeating[k], 1, 1)) {
                    printf("       %s under %s did not move\n",
                           incolorname[repeating[k]], colorfun[g]);
                    moved = 0;
                }
        check(moved, "every function but ident moves each of the modes that "
                     "repeat");

        /* sin(real)*sin(imag) reads z alone, a part at a time: bent by a
         * function, it is the picture with each part of z bent by it,
         * exactly, as the test does the arithmetic the function does. */
        check(sweep(18, 1, 1, 6) ==
                  sweep(18, 1, 1, 0, [](number_t x) { return x * x * x; }),
              "sin(real)*sin(imag) under cube is the picture at z cubed");
        check(sweep(18, 1, 1, 5) ==
                  sweep(18, 1, 1, 0,
                        [](number_t x) { return x > 0 ? x * x : -(x * x); }),
              "sin(real)*sin(imag) under square is the picture at z squared "
              "with its sign");
        /* and the frequency multiplies what the function gave */
        check(sweep(18, 2, 1, 6) ==
                  sweep(18, 1, 1, 0,
                        [](number_t x) { return 2 * (x * x * x); }),
              "at a frequency of two, twice the cube");

        int untouched = 1;
        for (int m = 0; m < INCOLORING; m++) {
            int reached = 0;
            for (int k = 0; k < 11; k++)
                if (repeating[k] == m)
                    reached = 1;
            if (reached || set_apart(m))
                continue;
            for (int g = 1; g < COLORFUN - 1; g++)
                if (sweep(m, 1, 1, g) != sweep(m, 1, 1)) {
                    printf("       %s moved under %s\n", incolorname[m],
                           colorfun[g]);
                    untouched = 0;
                }
        }
        check(untouched, "the modes that do not repeat draw the same under any "
                         "function");
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    return failures != 0;
}
