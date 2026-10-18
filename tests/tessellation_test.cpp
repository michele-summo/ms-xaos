/* The inside's tessellation mode: what it colours by and how its tiles stand
 * apart.
 *
 * Where it starts -- the angle of z, turned, with the decomposition settings
 * where they start -- it must draw what it drew before it had a choice, to
 * the bit: its checksums are pinned here, taken from the code before the
 * change, the same at both precisions. randsctile_place names every tile as
 * randsctile_name does and says where in it the point stands, which on
 * squares is plain arithmetic. Each setting moves what it belongs to and
 * nothing else: the decomposition settings reach the angles and not zmag or
 * the tiles alone, the angle of z-c reaches the angle of z and not the one
 * about the tile's middle, a flat tile is the same whatever it would have
 * been coloured by, and the frequency of the periodic modes is the tiles at
 * z scaled.
 *
 * formulas.cpp is included rather than linked, since incolor_output is its
 * own and static.
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

static void check(int ok, const char *text)
{
    printf("%s   %s\n", ok ? "ok  " : "FAIL", text);
    if (!ok)
        failures++;
}

/* The inside over orbits, z scaled by s; the grid is the one the fbm test
 * draws the tessellation over, so that the pins below are the same picture. */
static unsigned long long inside(number_t s = 1)
{
    cfractalc.incoloringmode = INCOLORING_TESSELLATION;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = -40; iy < 40; iy++)
        for (int ix = -40; ix < 40; ix++) {
            number_t zre = ((number_t)ix + (number_t)0.5) / 23 * s;
            number_t zim = ((number_t)iy + (number_t)0.5) / 29 * s;
            color_px = zre;
            color_py = zim;
            h ^= incolor_output(zre, zim, (number_t)0.3, (number_t)-0.4, 100);
            h *= 1099511628211ULL;
        }
    return h;
}

static void style(int colour, int fill)
{
    cfractalc.intessellcolor = colour;
    cfractalc.intessellfill = fill;
}

/* The inside's decomposition settings, as uih_setdecomposition takes them. */
static void decomp(int sectors, number_t rotation, int of, int radial,
                   number_t repeats, int way, number_t speed)
{
    cfractalc.indecompsectors = sectors;
    cfractalc.indecomprotation = rotation;
    cfractalc.indecompof = of;
    cfractalc.indecompradial = radial;
    cfractalc.indecomprepeats = repeats;
    cfractalc.indecompway = way;
    cfractalc.indecompspeed = speed;
}
static void plain(void) { decomp(0, 0, 0, 0, 1, 0, 1); }

/* Each decomposition setting moved on its own, the rest where they start. */
static void moved(int which)
{
    plain();
    switch (which) {
        case 0: cfractalc.indecompsectors = 6; break;
        case 1: cfractalc.indecomprotation = 30; break;
        case 2: cfractalc.indecompof = 1; break;
        case 3: cfractalc.indecompradial = 3; break;
        case 4: cfractalc.indecomprepeats = 2; break;
        case 5: cfractalc.indecompway = 2; break;
        case 6: cfractalc.indecompspeed = 2; break;
    }
}
static const char *const setting[] = {"sectors", "rotation", "angle of",
                                      "rings", "repetitions", "direction",
                                      "speed"};

int main()
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
    cfractalc.indecomprepeats = cfractalc.outdecomprepeats = 1;
    cfractalc.indecompspeed = cfractalc.outdecompspeed = 1;
    cfractalc.intessellation = 3;
    style(TESSELL_ANGLE, TESSELL_TURNED);
    color_change_prepare(&cfractalc);

    /* --- where it starts, what there was -------------------------------- */
    {
        static const struct {
            int tiling;
            unsigned long long sum;
        } was[] = {{1, 0x70a351df50f5f5fdULL},  {3, 0x2058fe816d9e4bfeULL},
                   {18, 0x5c9a2116d7ecbaa6ULL}, {38, 0xc651db659223f3dcULL},
                   {39, 0x2f13603e2dfffa57ULL}, {44, 0x856e53e54aeb74cdULL},
                   {45, 0xc9cd6499a83fb957ULL}};
        int same = 1;
        for (const auto &w : was) {
            cfractalc.intessellation = w.tiling;
            unsigned long long now = inside();
            if (now != w.sum) {
                printf("       tiling %d: %016llx, was %016llx\n", w.tiling,
                       now, w.sum);
                same = 0;
            }
        }
        check(same, "the angle of z, turned, draws what it drew before, to "
                    "the bit");
        cfractalc.intessellation = 18;
    }

    /* --- where in its tile a point stands -------------------------------- */
    {
        int names = 1, inrange = 1;
        for (int t = 1; t <= RANDSCTILE_TILINGS; t++)
            for (int k = 0; k < 400; k++) {
                number_t x = ((number_t)(k % 20) - 10 + (number_t)0.37) * 3;
                number_t y = ((number_t)(k / 20) - 10 + (number_t)0.61) * 3;
                uint64_t a = 0, b = 1;
                double out = -1, rx, ry;
                int found = randsctile_name(t, x, y, &a);
                if (found != randsctile_place(t, x, y, &b, &out, &rx, &ry) ||
                    a != b) {
                    if (names)
                        printf("       tiling %d at (%g, %g) named apart\n",
                               t, (double)x, (double)y);
                    names = 0;
                }
                if (!(out >= 0 && out <= 1))
                    inrange = 0;
            }
        check(names, "randsctile_place finds the tile randsctile_name does");
        check(inrange, "and the point stands between its middle and its "
                       "edge");

        /* squares of unit area: from the middle (k + 1/2, 1/2) the point is
         * where it was put, and as far out as the larger of the two, against
         * the half side */
        int square = 1;
        for (int k = -6; k < 6; k++)
            for (int j = 0; j < 9; j++) {
                const double dx = (j % 3 - 1) * 0.25, dy = (j / 3 - 1) * 0.125;
                uint64_t name;
                double out, rx, ry;
                if (!randsctile_place(1, (number_t)k + (number_t)0.5 + dx,
                                      (number_t)0.5 + dy, &name, &out, &rx,
                                      &ry) ||
                    fabs(rx - dx) > 1e-12 || fabs(ry - dy) > 1e-12 ||
                    fabs(out - 2 * (fabs(dx) > fabs(dy) ? fabs(dx) : fabs(dy))) >
                        1e-12)
                    square = 0;
            }
        check(square, "on squares that is the point less the square's "
                      "middle");
    }

    /* --- each setting moves what it belongs to --------------------------- */
    {
        unsigned long long base[4][4];
        int apart = 1;
        for (int c = 0; c < 4; c++)
            for (int f = 0; f < 4; f++) {
                style(c, f);
                base[c][f] = inside();
                for (int d = 0; d < c * 4 + f; d++)
                    if (base[d / 4][d % 4] == base[c][f] &&
                        f != TESSELL_FLAT && d % 4 != TESSELL_FLAT &&
                        !(c == TESSELL_NONE && f == TESSELL_TURNED)) {
                        printf("       %s, %s draws as %s, %s\n",
                               tessell_colors[c], tessell_fills[f],
                               tessell_colors[d / 4], tessell_fills[d % 4]);
                        apart = 0;
                    }
            }
        check(apart, "every colouring and every fill draws a picture of its "
                     "own");
        check(base[0][3] == base[1][3] && base[0][3] == base[2][3] &&
                  base[0][3] == base[3][3],
              "but a flat tile is the same whatever it would be coloured by");
        check(base[TESSELL_NONE][TESSELL_TURNED] == base[0][TESSELL_FLAT],
              "and a tile turned from nothing is a flat one");

        int reach = 1, stays = 1;
        for (int c = 0; c < 4; c++)
            for (int f = 0; f < 4; f++)
                for (int s = 0; s < 7; s++) {
                    moved(s);
                    style(c, f);
                    unsigned long long now = inside();
                    plain();
                    /* the angles read every setting, the one about the
                     * tile's middle all but whose angle it is; zmag, the
                     * tiles alone and a flat tile none */
                    int should = f != TESSELL_FLAT &&
                                 (c == TESSELL_ANGLE ||
                                  (c == TESSELL_MIDDLE && s != 2));
                    if (should && now == base[c][f]) {
                        printf("       %s does not reach %s, %s\n", setting[s],
                               tessell_colors[c], tessell_fills[f]);
                        reach = 0;
                    }
                    if (!should && now != base[c][f]) {
                        printf("       %s moves %s, %s\n", setting[s],
                               tessell_colors[c], tessell_fills[f]);
                        stays = 0;
                    }
                }
        check(reach, "every decomposition setting reaches the angles");
        check(stays, "and nothing else");

        /* the frequency of the periodic modes is the tiles at z scaled,
         * exactly, two scaling exactly: for the tiles alone, which read
         * nothing but the tile */
        int scaled = 1;
        for (int f = 0; f < 4; f++) {
            style(TESSELL_NONE, f);
            cfractalc.incolorfrequency = 2;
            unsigned long long twice = inside();
            cfractalc.incolorfrequency = 1;
            if (twice != inside(2))
                scaled = 0;
        }
        check(scaled, "frequency 2 is the tiles at z scaled by 2");
        style(TESSELL_ANGLE, TESSELL_TURNED);
    }

    /* --- the dialog's names ---------------------------------------------- */
    {
        int n = 0, numbered = 1;
        for (; tessell_tilings[n] != NULL; n++) {
            char at[8];
            snprintf(at, sizeof at, "%d: ", n + 1);
            if (strncmp(tessell_tilings[n], at, strlen(at)))
                numbered = 0;
        }
        check(n == RANDSCTILE_TILINGS && numbered,
              "the dialog names every tiling, after its number");
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    return failures != 0;
}
