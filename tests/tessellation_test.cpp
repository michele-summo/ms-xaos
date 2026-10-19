/* The tessellation modes: what they colour by, how their tiles stand apart,
 * and how far apart they are set.
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
 * The intensity at one is what there was, and at nought leaves only what the
 * colour comes from: zmag is zmag and the angle the decomposition, to the
 * bit. The outside's three modes, at nought with nothing to colour by, are
 * the count alone; sized by the count, a whole number of doublings is the
 * tiles at a fixed frequency that many times doubled, and in between a
 * palette of 256 entries takes one of the two sizes.
 *
 * formulas.cpp is included rather than linked, since incolor_output is its
 * own and static.
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

/* The outside over the same grid, the pixel's place the point itself and a
 * count of its own at each, in the mode given. */
static unsigned int outside_at(int ix, int iy)
{
    number_t x = ((number_t)ix + (number_t)0.5) / 23;
    number_t y = ((number_t)iy + (number_t)0.5) / 29;
    color_px = x;
    color_py = y;
    return color_output(x * 3 + 1, y * 3 - 2,
                        (unsigned int)((ix * 7 + iy * 3 + 400) % 23));
}
static unsigned long long outside(OutColormodeType::OutColormode mode)
{
    cfractalc.coloringmode = mode;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = -40; iy < 40; iy++)
        for (int ix = -40; ix < 40; ix++) {
            h ^= outside_at(ix, iy);
            h *= 1099511628211ULL;
        }
    return h;
}
static void outstyle(int tiling, int colour, int fill, number_t intensity,
                     number_t frequency, int scaling, number_t factor)
{
    cfractalc.outtessellation = tiling;
    cfractalc.outtessellcolor = colour;
    cfractalc.outtessellfill = fill;
    cfractalc.outtessellintensity = intensity;
    cfractalc.outtessellfrequency = frequency;
    cfractalc.outtessellscaling = scaling;
    cfractalc.outtessellfactor = factor;
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
    cfractalc.intessellintensity = 1;
    cfractalc.outtessellation = 3;
    cfractalc.outtessellintensity = 1;
    cfractalc.outtessellfrequency = 8;
    cfractalc.outtessellfactor = 1;
    style(TESSELL_ANGLE, TESSELL_TURNED);
    color_change_prepare(&cfractalc);

    /* --- where it starts, what there was -------------------------------- */
    {
        /* the same at 64 and 113 bits; with a long double that is double, as
         * on Apple Silicon, tiling 18 puts a few points in other tiles, and
         * its sum there was taken from the old code built with double */
        static const struct {
            int tiling;
            unsigned long long sum, sum53;
        } was[] = {{1, 0x70a351df50f5f5fdULL, 0x70a351df50f5f5fdULL},
                   {3, 0x2058fe816d9e4bfeULL, 0x2058fe816d9e4bfeULL},
                   {18, 0x5c9a2116d7ecbaa6ULL, 0x14dd02c8e49acd02ULL},
                   {38, 0xc651db659223f3dcULL, 0xc651db659223f3dcULL},
                   {39, 0x2f13603e2dfffa57ULL, 0x2f13603e2dfffa57ULL},
                   {44, 0x856e53e54aeb74cdULL, 0x856e53e54aeb74cdULL},
                   {45, 0xc9cd6499a83fb957ULL, 0xc9cd6499a83fb957ULL}};
        int same = 1;
        for (const auto &w : was) {
            cfractalc.intessellation = w.tiling;
            unsigned long long now = inside();
            const unsigned long long sum =
                NUMBER_MANTISSA_BITS == 53 ? w.sum53 : w.sum;
            if (now != sum) {
                printf("       tiling %d: %016llx, was %016llx\n", w.tiling,
                       now, sum);
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

    /* --- the inside's intensity ------------------------------------------ */
    {
        /* at nought, what the colour comes from alone */
        int zmag = 1, angle = 1, notile = 1;
        cfractalc.intessellintensity = 0;
        for (int t = 1; t <= RANDSCTILE_TILINGS; t++) {
            cfractalc.intessellation = t;
            style(TESSELL_ZMAG, TESSELL_TURNED);
            unsigned long long tz = inside();
            cfractalc.incoloringmode = 1;
            unsigned long long hz = 1469598103934665603ULL;
            for (int iy = -40; iy < 40; iy++)
                for (int ix = -40; ix < 40; ix++) {
                    number_t zre = ((number_t)ix + (number_t)0.5) / 23;
                    number_t zim = ((number_t)iy + (number_t)0.5) / 29;
                    hz ^= incolor_output(zre, zim, (number_t)0.3,
                                         (number_t)-0.4, 100);
                    hz *= 1099511628211ULL;
                }
            if (tz != hz)
                zmag = 0;
            /* the angle with its settings moved, which is when
             * Decomposition-like is worked as the tessellation's angle is */
            decomp(6, 15, 1, 3, 2, 2, (number_t)1.5);
            style(TESSELL_ANGLE, TESSELL_TURNED);
            unsigned long long ta = inside();
            cfractalc.incoloringmode = 2;
            unsigned long long ha = 1469598103934665603ULL;
            for (int iy = -40; iy < 40; iy++)
                for (int ix = -40; ix < 40; ix++) {
                    number_t zre = ((number_t)ix + (number_t)0.5) / 23;
                    number_t zim = ((number_t)iy + (number_t)0.5) / 29;
                    ha ^= incolor_output(zre, zim, (number_t)0.3,
                                         (number_t)-0.4, 100);
                    ha *= 1099511628211ULL;
                }
            plain();
            if (ta != ha)
                angle = 0;
        }
        check(zmag, "at intensity nought zmag is zmag, on every tiling");
        check(angle, "and the angle the decomposition, on every tiling");
        /* nor does the tiling show, but in the angle about its middle */
        for (int c = 0; c < 4; c++)
            for (int f = 0; f < 4; f++) {
                if (c == TESSELL_MIDDLE && f != TESSELL_FLAT)
                    continue;
                style(c, f);
                cfractalc.intessellation = 3;
                unsigned long long a = inside();
                cfractalc.intessellation = 18;
                if (a != inside()) {
                    printf("       %s, %s shows its tiles at nought\n",
                           tessell_colors[c], tessell_fills[f]);
                    notile = 0;
                }
            }
        check(notile, "nor do the tiles, but about their middles");

        /* away from one every colouring and fill moves */
        int moves = 1;
        for (int c = 0; c < 4; c++)
            for (int f = 0; f < 4; f++) {
                style(c, f);
                cfractalc.intessellintensity = 1;
                unsigned long long one = inside();
                cfractalc.intessellintensity = (number_t)0.5;
                if (one == inside()) {
                    printf("       %s, %s does not move\n", tessell_colors[c],
                           tessell_fills[f]);
                    moves = 0;
                }
            }
        check(moves, "and away from one every colouring and fill moves");
        cfractalc.intessellintensity = 1;
        style(TESSELL_ANGLE, TESSELL_TURNED);
    }

    /* --- the outside ----------------------------------------------------- */
    {
        typedef OutColormodeType M;
        cfractalc.outcolorspeed = 1;
        outstyle(3, TESSELL_NONE, TESSELL_TURNED, 0, 8, FBM_SCALE_FIXED, 1);
        const unsigned long long iter = outside(M::ColOut_iter);
        int count = 1;
        for (int t = 1; t <= RANDSCTILE_TILINGS; t += 11)
            for (int s = 0; s < 3; s++) {
                outstyle(t, TESSELL_NONE, TESSELL_TURNED, 0, 8, s, 1);
                if (outside(M::ColOut_tessell_iter) != iter)
                    count = 0;
            }
        check(count, "at intensity nought, by nothing, tessellation + iter is "
                     "iter");
        outstyle(18, TESSELL_NONE, TESSELL_TURNED, 0, 8, FBM_SCALE_FIXED, 1);
        cfractalc.coloringmode = M::ColOut_tessell;
        unsigned int first = outside_at(0, 0);
        int flat = 1;
        for (int k = -40; k < 40; k++)
            if (outside_at(k, k / 2) != first)
                flat = 0;
        check(flat, "and tessellation alone one colour");

        /* each mode, colouring and fill a picture of its own at one */
        int apart = 1;
        static const M::OutColormode modes[3] = {
            M::ColOut_tessell_smooth, M::ColOut_tessell_iter, M::ColOut_tessell};
        for (int c = 0; c < 4; c++)
            for (int f = 0; f < 3; f++) {
                outstyle(3, c, f, 1, 8, FBM_SCALE_FIXED, 1);
                unsigned long long a = outside(M::ColOut_tessell);
                outstyle(3, c, f, (number_t)0.5, 8, FBM_SCALE_FIXED, 1);
                if (a == outside(M::ColOut_tessell))
                    apart = 0;
            }
        outstyle(3, TESSELL_ANGLE, TESSELL_TURNED, 1, 8, FBM_SCALE_FIXED, 1);
        if (outside(modes[1]) == outside(modes[2]) ||
            outside(modes[1]) == iter)
            apart = 0;
        check(apart, "at one the tiles show, and the intensity moves them");

        /* the decomposition settings reach the outside's angles alone */
        int reach = 1;
        for (int c = 0; c < 4; c++) {
            outstyle(3, c, TESSELL_TURNED, 1, 8, FBM_SCALE_FIXED, 1);
            unsigned long long a = outside(M::ColOut_tessell);
            cfractalc.outdecompsectors = 6;
            unsigned long long b = outside(M::ColOut_tessell);
            cfractalc.outdecompsectors = 0;
            if ((c == TESSELL_ANGLE || c == TESSELL_MIDDLE) != (a != b))
                reach = 0;
        }
        check(reach, "the outside's decomposition settings reach its angles "
                     "and nothing else");

        /* the inside's intensity is not the outside's, nor the other way */
        outstyle(3, TESSELL_ANGLE, TESSELL_SHADED, 1, 8, FBM_SCALE_FIXED, 1);
        style(TESSELL_ANGLE, TESSELL_SHADED);
        unsigned long long o = outside(M::ColOut_tessell),
                           i = inside();
        cfractalc.intessellintensity = 3;
        int own = outside(M::ColOut_tessell) == o;
        cfractalc.intessellintensity = 1;
        cfractalc.outtessellintensity = 3;
        own = own && inside() == i;
        check(own, "each side has an intensity of its own");
        style(TESSELL_ANGLE, TESSELL_TURNED);

        /* sized by the count: fixed is a factor of nought, the smooth count
         * where it is whole the whole one, and a whole number of doublings
         * the fixed frequency that many times doubled */
        outstyle(3, TESSELL_ANGLE, TESSELL_TURNED, 1, 8, FBM_SCALE_FIXED, 1);
        const unsigned long long fixed = outside(M::ColOut_tessell);
        outstyle(3, TESSELL_ANGLE, TESSELL_TURNED, 1, 8, FBM_SCALE_ITER, 0);
        int nought = outside(M::ColOut_tessell) == fixed;
        outstyle(3, TESSELL_ANGLE, TESSELL_TURNED, 1, 8, FBM_SCALE_SMOOTH, 0);
        nought = nought && outside(M::ColOut_tessell) == fixed;
        check(nought, "at a factor of nought the count changes nothing");
        outstyle(3, TESSELL_ANGLE, TESSELL_TURNED, 1, 8, FBM_SCALE_ITER, 3);
        unsigned long long whole = outside(M::ColOut_tessell);
        cfractalc.outtessellscaling = FBM_SCALE_SMOOTH;
        check(outside(M::ColOut_tessell) == whole && whole != fixed,
              "the smooth count where it is whole sizes as the whole one");

        int doubled = 1, between = 1;
        cpalette.type = C256;
        cfractalc.coloringmode = M::ColOut_tessell;
        for (int iy = -40; iy < 40; iy += 3)
            for (int ix = -40; ix < 40; ix += 3) {
                const int n = (ix * 7 + iy * 3 + 400) % 23;
                outstyle(5, TESSELL_ANGLE, TESSELL_TURNED, 1, 2,
                         FBM_SCALE_ITER, 10);
                unsigned int scaled = outside_at(ix, iy);
                outstyle(5, TESSELL_ANGLE, TESSELL_TURNED, 1,
                         nldexp((number_t)2, n), FBM_SCALE_FIXED, 1);
                if (scaled != outside_at(ix, iy))
                    doubled = 0;
                /* three doublings in ten iterations: between two sizes,
                 * one or the other */
                outstyle(5, TESSELL_ANGLE, TESSELL_TURNED, 1, 2,
                         FBM_SCALE_ITER, 3);
                unsigned int mid = outside_at(ix, iy);
                const int w = (int)nfloor((number_t)3 * n / 10);
                outstyle(5, TESSELL_ANGLE, TESSELL_TURNED, 1,
                         nldexp((number_t)2, w), FBM_SCALE_FIXED, 1);
                unsigned int lo = outside_at(ix, iy);
                cfractalc.outtessellfrequency = nldexp((number_t)2, w + 1);
                unsigned int hi = outside_at(ix, iy);
                if (mid != lo && mid != hi)
                    between = 0;
            }
        cpalette.type = 0;
        check(doubled, "a whole number of doublings is the tiles at that "
                       "frequency");
        check(between, "and between two, with 256 colours, one of the two");
        outstyle(3, TESSELL_ANGLE, TESSELL_TURNED, 1, 8, FBM_SCALE_FIXED, 1);
    }

    /* --- the angle and zmag across an edge of the iterations ------------
     *
     * Read at the escape point, which one pass more takes to its square, the
     * angle and zmag jumped at every edge of the iterations. With a smooth
     * count, the colour where the point stands just past the line the passes
     * stop at (tau nought) is the colour of the point one pass on, and on the
     * other side of the edge (tau one) the colour of that side's own point,
     * which is that same point: the two must agree. tau one is the colour as
     * it was, and the colourings that do not read z never notice tau. */
    {
        cfractalc.coloringmode = OutColormodeType::ColOut_tessell;
        cpalette.type = C256;
        auto at = [](number_t zre, number_t zim, number_t tau, unsigned int n) {
            color_smooth = (number_t)n * 256;
            color_smooth_held = 1;
            smooth_tau = tau;
            smooth_degree = 2;
            return color_output(zre, zim, n);
        };
        int agree = 1, as_was = 1, apart = 1, moved = 0;
        for (int c : {TESSELL_ANGLE, TESSELL_ZMAG, TESSELL_MIDDLE,
                      TESSELL_NONE}) {
            outstyle(3, c, TESSELL_SHADED, 1, 8, FBM_SCALE_FIXED, 1);
            const int reads_z = c == TESSELL_ANGLE || c == TESSELL_ZMAG;
            for (int iy = -6; iy < 6; iy++)
                for (int ix = -6; ix < 6; ix++) {
                    color_px = ((number_t)ix + (number_t)0.5) / 7;
                    color_py = ((number_t)iy + (number_t)0.5) / 5;
                    number_t zre = ((number_t)ix + (number_t)0.5) * 300;
                    number_t zim = ((number_t)iy + (number_t)0.5) * 170;
                    unsigned int plain = color_output(zre, zim, 9);
                    unsigned int near = at(zre, zim, 0, 9);
                    unsigned int far = at(zre * zre - zim * zim,
                                          2 * zre * zim, 1, 10);
                    if (reads_z && near != far)
                        agree = 0;
                    if (at(zre, zim, 1, 9) != plain)
                        as_was = 0;
                    if (!reads_z && near != plain)
                        apart = 0;
                    if (reads_z && near != plain)
                        moved = 1;
                }
        }
        smooth_tau = -1;
        cpalette.type = 0;
        check(agree, "the angle and zmag agree across an edge of the "
                     "iterations");
        check(as_was && moved, "the point's own colour at one, another "
                               "below");
        check(apart, "tile-angle and none never read where between the "
                     "passes the point stands");
        outstyle(3, TESSELL_ANGLE, TESSELL_TURNED, 1, 8, FBM_SCALE_FIXED, 1);
    }

    /* --- an intensity between nought and one, dissolved ------------------
     *
     * With a palette of colours between its entries, an intensity between
     * nought and one is the picture at nought dissolved into the picture at
     * one by that much: each channel of the pixel lies between the two, and
     * the tiles are not moved along the palette, which took a tile to
     * another colour already at a tenth. Inside and out. */
    {
        const int was = cpalette.type, mode = cfractalc.incoloringmode;
        const struct truec keep = cpalette.info.truec;
        cfractalc.incoloringmode = INCOLORING_TESSELLATION;
        cpalette.type = TRUECOLOR;
        cpalette.info.truec.rmask = 0xff0000;
        cpalette.info.truec.gmask = 0x00ff00;
        cpalette.info.truec.bmask = 0x0000ff;
        cpalette.info.truec.allmask = 0xffffff;
        auto within = [](unsigned int v, unsigned int a, unsigned int b) {
            for (int shift = 0; shift < 24; shift += 8) {
                int x = (v >> shift) & 255, p = (a >> shift) & 255,
                    q = (b >> shift) & 255;
                int lo = p < q ? p : q, hi = p < q ? q : p;
                if (x < lo - 1 || x > hi + 1)
                    return 0;
            }
            return 1;
        };
        int between_in = 1, moved_in = 0, between_out = 1, moved_out = 0;
        for (int c : {TESSELL_ANGLE, TESSELL_ZMAG})
            for (int f : {TESSELL_TURNED, TESSELL_SHADED, TESSELL_EDGED}) {
                style(c, f);
                outstyle(3, c, f, 1, 8, FBM_SCALE_FIXED, 1);
                cfractalc.coloringmode = OutColormodeType::ColOut_tessell;
                for (int iy = -8; iy < 8; iy++)
                    for (int ix = -8; ix < 8; ix++) {
                        number_t x = ((number_t)ix + (number_t)0.5) / 23;
                        number_t y = ((number_t)iy + (number_t)0.5) / 29;
                        color_px = x;
                        color_py = y;
                        unsigned int in[3], out[3];
                        const number_t much[3] = {0, 1, (number_t)0.25};
                        for (int k = 0; k < 3; k++) {
                            cfractalc.intessellintensity = much[k];
                            cfractalc.outtessellintensity = much[k];
                            in[k] = incolor_output(x, y, (number_t)0.3,
                                                   (number_t)-0.4, 100);
                            out[k] = color_output(x * 3 + 1, y * 3 - 2, 7);
                        }
                        if (!within(in[2], in[0], in[1]))
                            between_in = 0;
                        if (in[2] != in[0] && in[2] != in[1])
                            moved_in = 1;
                        if (!within(out[2], out[0], out[1]))
                            between_out = 0;
                        if (out[2] != out[0] && out[2] != out[1])
                            moved_out = 1;
                    }
            }
        cfractalc.intessellintensity = 1;
        cfractalc.incoloringmode = mode;
        style(TESSELL_ANGLE, TESSELL_TURNED);
        outstyle(3, TESSELL_ANGLE, TESSELL_TURNED, 1, 8, FBM_SCALE_FIXED, 1);
        cpalette.type = was;
        cpalette.info.truec = keep;
        check(between_in && moved_in, "inside, an intensity of a quarter is "
                                      "the two pictures dissolved");
        check(between_out && moved_out, "and outside");
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
