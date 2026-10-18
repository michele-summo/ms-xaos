/* The settings of the decompositions.
 *
 * Binary decomposition and color decomposition outside, sign(imag),
 * Decomposition-like and fbm + decomposition inside. With every setting where
 * it starts they must draw what they drew before there were any, to the bit,
 * since every position saved until now has them so. There the code that was
 * there before runs as it was: its checksums are pinned here, the same at
 * both precisions, and the pictures were compared with 2.5's, pixel for
 * pixel. Away from it each setting must do what it
 * says, and the rings of the inside go through k*f as every periodic argument
 * of the inside does.
 *
 * formulas.cpp is included rather than linked, since color_output and
 * incolor_output are its own and static.
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

/* What an inside mode gives over a grid of orbits. The grid is a half step
 * off both axes, where a half-plane and a sector meet, so that two ways of
 * cutting the same half-planes agree on every point of it. */
static unsigned long long inside(int m)
{
    cfractalc.incoloringmode = m;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = -60; iy < 60; iy++)
        for (int ix = -60; ix < 60; ix++) {
            number_t zre = ((number_t)ix + (number_t)0.5) / 23;
            number_t zim = ((number_t)iy + (number_t)0.5) / 29;
            number_t pre = (number_t)0.3 + zre / 7;
            number_t pim = (number_t)-0.4 + zim / 5;
            h ^= incolor_output(zre, zim, pre, pim, 100);
            h *= 1099511628211ULL;
        }
    return h;
}

/* And an outside one, over escape points and counts. */
static unsigned long long outside(OutColormodeClass::OutColormode m)
{
    cfractalc.coloringmode = m;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = -60; iy < 60; iy++)
        for (int ix = -60; ix < 60; ix++) {
            number_t zre = ((number_t)ix + (number_t)0.5) / 9;
            number_t zim = ((number_t)iy + (number_t)0.5) / 11;
            unsigned int iter = 1 + (unsigned int)((ix + 60) * 7 + iy + 60) % 97;
            h ^= color_output(zre, zim, iter);
            h *= 1099511628211ULL;
        }
    return h;
}

static void inset(int sectors, number_t rotation, int of, int radial)
{
    cfractalc.indecompsectors = sectors;
    cfractalc.indecomprotation = rotation;
    cfractalc.indecompof = of;
    cfractalc.indecompradial = radial;
}

static void outset(int sectors, number_t rotation, int radial)
{
    cfractalc.outdecompsectors = sectors;
    cfractalc.outdecomprotation = rotation;
    cfractalc.outdecompradial = radial;
}

/* And how the colour ones walk the palette round a turn, both sides at once. */
static void walk(number_t repeats, int way, number_t speed)
{
    cfractalc.indecomprepeats = cfractalc.outdecomprepeats = repeats;
    cfractalc.indecompway = cfractalc.outdecompway = way;
    cfractalc.indecompspeed = cfractalc.outdecompspeed = speed;
}

static const OutColormodeClass::OutColormode binary =
    OutColormodeType::ColOut_binary_decomposition;
static const OutColormodeClass::OutColormode colour =
    OutColormodeType::ColOut_color_decomposition;

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
     * speed change at one, the frequency at one and the function ident */
    cfractalc.incolorspeed = cfractalc.outcolorspeed = 1;
    cfractalc.incolorfuninitial = cfractalc.incolorfuninitialspeed = 1;
    cfractalc.incolorfunspeed = cfractalc.incolorfun2speed = 1;
    cfractalc.outcolorfuninitial = cfractalc.outcolorfuninitialspeed = 1;
    cfractalc.outcolorfunspeed = cfractalc.outcolorfun2speed = 1;
    cfractalc.incolorfrequency = 1;
    cfractalc.infbmintensity = 4;
    cfractalc.infbmfrequency = 8;
    cfractalc.infbmoctaves = 4;
    cfractalc.infbmroughness = (number_t)0.5;
    cfractalc.infbmseed = 1;
    cfractalc.indecomprepeats = cfractalc.outdecomprepeats = 1;
    cfractalc.indecompspeed = cfractalc.outdecompspeed = 1;
    color_change_prepare(&cfractalc);

    /* --- where they start, what there was ------------------------------ */
    {
        check(decomp_plain(&cfractalc, 0) && decomp_plain(&cfractalc, 1),
              "every setting at nought is where they start");
        inset(0, 0, 1, 0);
        check(!decomp_plain(&cfractalc, 1), "the angle of z-c is not");
        inset(0, 0, 0, 0);
        outset(0, 360, 0);
        check(!decomp_plain(&cfractalc, 0), "a rotation is not");
        outset(0, 0, 0);

        static const unsigned long long was[5] = {
            0x1373a4e4a39e79a3ULL, /* binary decomposition */
            0x31917341cab94b5eULL, /* color decomposition */
            0xc5f773ac54b997afULL, /* Decomposition-like */
            0x426819746f8dfba3ULL, /* sign(imag) */
            0x81d380a45720b21cULL, /* fbm + decomposition */
        };
        unsigned long long now[5] = {outside(binary), outside(colour),
                                     inside(2), inside(19),
                                     inside(INCOLORING_FBM_DECOMP)};
        static const char *const name[5] = {
            "binary decomposition", "color decomposition",
            "Decomposition-like", "sign(imag)", "fbm + decomposition"};
        for (int k = 0; k < 5; k++) {
            sprintf(what, "%s with its settings where they start draws as it "
                          "did (%#llx)",
                    name[k], now[k]);
            check(now[k] == was[k], what);
        }
    }

    /* --- each setting does what it says ---------------------------------- */
    {
        /* A full turn is no turn: the half-planes cut the other way, as
         * sectors, land where they always were. */
        unsigned long long b = outside(binary), s = inside(19);
        outset(0, 360, 0);
        inset(0, 360, 0, 0);
        check(outside(binary) == b,
              "binary decomposition turned a full turn is the half-planes");
        check(inside(19) == s, "and so is sign(imag)");
        outset(0, 0, 0);
        inset(0, 0, 0, 0);

        /* Two sectors of a colour would meet: an odd number is the even one
         * above it. */
        outset(3, 0, 0);
        unsigned long long three = outside(binary);
        outset(4, 0, 0);
        check(three == outside(binary),
              "three sectors of binary decomposition are four");
        outset(0, 0, 0);

        /* Every setting moves every mode it belongs to. */
        struct {
            const char *text;
            int in;
            int sectors;
            number_t rotation;
            int of, radial;
        } moves[] = {
            {"sectors", 0, 8, 0, 0, 0},  {"a rotation", 0, 0, 30, 0, 0},
            {"bands", 0, 0, 0, 0, 3},    {"sectors", 1, 6, 0, 0, 0},
            {"a rotation", 1, 0, 30, 0, 0}, {"the angle of z-c", 1, 0, 0, 1, 0},
            {"the angle of z*c", 1, 0, 0, 2, 0}, {"rings", 1, 0, 0, 0, 8},
        };
        for (auto &m : moves) {
            if (m.in) {
                unsigned long long a = inside(2), c = inside(19),
                                   d = inside(INCOLORING_FBM_DECOMP);
                inset(m.sectors, m.rotation, m.of, m.radial);
                sprintf(what, "%s moves the three decompositions of the inside",
                        m.text);
                check(inside(2) != a && inside(19) != c &&
                          inside(INCOLORING_FBM_DECOMP) != d,
                      what);
                inset(0, 0, 0, 0);
            } else {
                unsigned long long a = outside(binary), c = outside(colour);
                outset(m.sectors, m.rotation, m.radial);
                sprintf(what, "%s moves both decompositions of the outside",
                        m.text);
                check(outside(binary) != a && outside(colour) != c, what);
                outset(0, 0, 0);
            }
        }

        /* The walk moves every decomposition, each part of it on its own. */
        {
            struct {
                const char *text;
                number_t repeats;
                int way;
                number_t speed;
            } walks[] = {{"two and a half repeats", (number_t)2.5, 0, 1},
                         {"walking left", 1, 1, 1},
                         {"walking both ways", 1, 2, 1},
                         {"a speed of two", 1, 0, 2}};
            unsigned long long oc = outside(colour), i2 = inside(2),
                               i25 = inside(INCOLORING_FBM_DECOMP),
                               ob = outside(binary), i19 = inside(19);
            for (auto &w : walks) {
                walk(w.repeats, w.way, w.speed);
                sprintf(what, "%s moves the five decompositions", w.text);
                check(outside(colour) != oc && inside(2) != i2 &&
                          inside(INCOLORING_FBM_DECOMP) != i25 &&
                          outside(binary) != ob && inside(19) != i19,
                      what);
            }
            walk(1, 0, 1);

            /* The speed is how far along the palette a repeat goes: at
             * nought every point stays where the walk starts, whichever way
             * it would have gone. */
            walk(3, 0, 0);
            unsigned long long still = outside(colour);
            walk(3, 1, 0);
            check(outside(colour) == still,
                  "at a speed of nought the way does not matter");

            /* For the binary ones the speed is how far apart the colours of
             * alternate sectors are: at nought, none -- binary decomposition
             * is the plain count, sign(imag) one tone. */
            walk(3, 2, 0);
            check(outside(binary) ==
                      outside(OutColormodeType::ColOut_iter),
                  "binary decomposition at a speed of nought is the count");
            unsigned long long one = inside(19);
            walk(1, 0, 0);
            check(inside(19) == one,
                  "sign(imag) at a speed of nought is one tone, however it "
                  "is cut");
            walk(1, 0, 1);
        }

        /* The rings go through k*f: eight to a unit at a frequency of two
         * are sixteen at one, exactly, since the frequency and the number of
         * rings are powers of two. */
        inset(8, 0, 0, 8);
        cfractalc.incolorfrequency = 2;
        unsigned long long k2 = inside(19);
        cfractalc.incolorfrequency = 1;
        inset(8, 0, 0, 16);
        check(inside(19) == k2,
              "eight rings at a frequency of two are sixteen at one");
        /* and the function bends them */
        cfractalc.incolorperiodicfun = 6;
        check(inside(19) != k2, "the function bends the rings");
        cfractalc.incolorperiodicfun = 0;
        inset(0, 0, 0, 0);

        /* The settings of one side are nothing to the other, nor to the
         * modes that are no decomposition. */
        unsigned long long zmag = inside(1), iter =
            outside(OutColormodeType::ColOut_iter);
        unsigned long long b2 = outside(binary), i2 = inside(2);
        inset(6, 30, 1, 8);
        check(outside(binary) == b2,
              "the inside's settings leave the outside alone");
        check(inside(1) == zmag, "and zmag");
        inset(0, 0, 0, 0);
        outset(6, 30, 3);
        check(inside(2) == i2, "the outside's leave the inside alone");
        check(outside(OutColormodeType::ColOut_iter) == iter, "and iter");
        outset(0, 0, 0);
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    return failures != 0;
}
