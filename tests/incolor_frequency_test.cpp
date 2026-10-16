/* The frequency of the inside modes that repeat.
 *
 * At one, the six modes it reaches must give what they gave before there was
 * such a thing, to the bit, since every position saved until now draws at
 * one. Away from one it must do what it says -- make what repeats repeat that
 * much faster -- which for the modes that read z alone is the same as drawing
 * them at z times the frequency, exactly. And the modes that do not repeat
 * must never read it.
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

/* The six that repeat, by their numbers in incolorname. */
static const int periodic[6] = {5, 6, 7, 9, 18, 20};

/* What mode m gives at frequency f over a grid of orbits, z scaled by s. */
static unsigned long long sweep(int m, number_t f, number_t s)
{
    cfractalc.incoloringmode = m;
    cfractalc.incolorfrequency = f;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = -60; iy <= 60; iy++)
        for (int ix = -60; ix <= 60; ix++) {
            number_t zre = (number_t)ix / 23, zim = (number_t)iy / 29;
            number_t pre = (number_t)0.3 + zre / 7;
            number_t pim = (number_t)-0.4 + zim / 5;
            unsigned int r = incolor_output(zre * s, zim * s, pre, pim, 100);
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
     * z twice as far. */
    {
        struct {
            int mode;
            number_t f, s;
        } same[4] = {{7, 4, 2}, {9, 2, 2}, {18, 2, 2}, {20, 2, 2}};
        for (int k = 0; k < 4; k++) {
            sprintf(what,
                    "%s at a frequency of %g is the picture at one with z "
                    "%g times as far",
                    incolorname[same[k].mode], (double)same[k].f,
                    (double)same[k].s);
            check(sweep(same[k].mode, same[k].f, 1) ==
                      sweep(same[k].mode, 1, same[k].s),
                  what);
        }
        /* and the two that read c as well: something else than at one */
        for (int m : {5, 6}) {
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
            for (int k = 0; k < 6; k++)
                if (periodic[k] == m)
                    repeats = 1;
            if (repeats || m == INCOLORING_TRUECOLOR || m >= INCOLORING_FBM)
                continue;
            if (sweep(m, 3, 1) != sweep(m, 1, 1)) {
                printf("       %s moved\n", incolorname[m]);
                untouched = 0;
            }
        }
        check(untouched, "the modes that do not repeat draw the same at any "
                         "frequency");
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    return failures != 0;
}
