/* The ways a palette can be made.
 *
 * Three of them scatter colours between black and white anchors, which is the
 * look XaoS has always had; four more pick their colours in relation to each
 * other. What has to hold of all seven is the same:
 *
 *  - a palette is made from a seed, and the same algorithm and seed must give
 *    the same palette, because that pair is all a saved position records of
 *    its colours -- get it wrong and a position comes back in the wrong ones;
 *  - a palette must actually use its entries, and no two algorithms may be
 *    the same algorithm under two numbers;
 *  - and every one of them must have some dark and some light in it. A palette
 *    with no range shows a fractal as one flat wash, which is the failure the
 *    new four came closest to: colours drawn from one narrow arc of the hue
 *    circle came out four shades of the same thing until the light and the
 *    dark were settled by position rather than by the dice.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "filter.h"

const char *qt_gettext(const char *, const char *t) { return t; }

static int failures = 0;

static void check(int ok, const char *what)
{
    printf("%-6s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        failures++;
}

/* One palette, as the entries an algorithm and a seed produce. */
struct made {
    int size;
    unsigned char rgb[4096][3];
};

static struct made *building;

static int alloccolor(struct palette *pal, int init, int r, int g, int b)
{
    if (init)
        pal->size = 0;
    if (pal->size >= pal->maxentries || pal->size >= 4096)
        return -1;
    building->rgb[pal->size][0] = (unsigned char)r;
    building->rgb[pal->size][1] = (unsigned char)g;
    building->rgb[pal->size][2] = (unsigned char)b;
    pal->pixels[pal->size] = (unsigned int)pal->size;
    pal->size += 1;
    return pal->size - 1;
}

static unsigned int pixels[4096];

static void make(struct made *out, int algorithm, int seed)
{
    struct palette pal;
    memset(&pal, 0, sizeof(pal));
    pal.maxentries = 256;
    pal.end = 255;
    pal.pixels = pixels;
    pal.alloccolor = alloccolor;
    memset(out, 0, sizeof(*out));
    building = out;
    mkpalette(&pal, seed, algorithm);
    out->size = pal.size;
}

static int identical(const struct made *a, const struct made *b)
{
    return a->size == b->size &&
           !memcmp(a->rgb, b->rgb, (size_t)a->size * 3);
}

/* How much of the hue circle a palette actually visits, as a fraction of it.
 *
 * Every entry with enough colour in it to have a hue is placed on the circle,
 * and what is measured is the circle less its largest empty arc. One hue gives
 * nearly nothing; two opposite ones give a half; a cycle gives one. The four
 * new ways to make a palette were meant to be schemes of several colours and
 * two of them were quietly one colour each, which this is here to notice. */
static double hue_spread(const struct made *m)
{
    double hue[4096];
    int n = 0;
    for (int i = 0; i < m->size && n < 4096; i++) {
        int r = m->rgb[i][0], g = m->rgb[i][1], b = m->rgb[i][2];
        int hi = r > g ? (r > b ? r : b) : (g > b ? g : b);
        int lo = r < g ? (r < b ? r : b) : (g < b ? g : b);
        if (hi - lo < 40)
            continue; /* too near grey to have a hue worth placing */
        double d = hi - lo, h;
        if (hi == r)
            h = (g - b) / d;
        else if (hi == g)
            h = 2 + (b - r) / d;
        else
            h = 4 + (r - g) / d;
        h /= 6;
        if (h < 0)
            h += 1;
        hue[n++] = h;
    }
    if (n < 2)
        return 0;
    /* sort, then the widest gap between neighbours around the circle */
    for (int i = 1; i < n; i++) {
        double v = hue[i];
        int j = i - 1;
        while (j >= 0 && hue[j] > v) {
            hue[j + 1] = hue[j];
            j--;
        }
        hue[j + 1] = v;
    }
    double widest = hue[0] + 1 - hue[n - 1];
    for (int i = 1; i < n; i++)
        if (hue[i] - hue[i - 1] > widest)
            widest = hue[i] - hue[i - 1];
    return 1 - widest;
}

/* The lightest and the darkest entry, as the sum of their three channels. */
static void range(const struct made *m, int *lo, int *hi)
{
    *lo = 766;
    *hi = -1;
    for (int i = 0; i < m->size; i++) {
        int v = m->rgb[i][0] + m->rgb[i][1] + m->rgb[i][2];
        if (v < *lo)
            *lo = v;
        if (v > *hi)
            *hi = v;
    }
}

int main(void)
{
    static struct made a, b;
    static const int seeds[6] = {1, 99, 777, 4242, 12345, 65535};
    char what[128];

    check(PALGORITHMS == 7, "there are seven ways to make a palette");

    /* --- the same seed gives the same palette ----------------------------- */
    {
        int steady = 1;
        for (int alg = 0; alg < PALGORITHMS; alg++)
            for (int s = 0; s < 6; s++) {
                make(&a, alg, seeds[s]);
                make(&b, alg, seeds[s]);
                if (!identical(&a, &b))
                    steady = 0;
            }
        check(steady, "a seed makes the same palette every time it is used");
    }

    /* --- and a different seed a different one ----------------------------- */
    {
        int varies = 1;
        for (int alg = 0; alg < PALGORITHMS; alg++) {
            make(&a, alg, seeds[0]);
            make(&b, alg, seeds[1]);
            if (identical(&a, &b))
                varies = 0;
        }
        check(varies, "and a different seed a different palette");
    }

    /* --- no algorithm is another one under a second number ---------------- */
    {
        int distinct = 1;
        for (int i = 0; i < PALGORITHMS && distinct; i++)
            for (int j = i + 1; j < PALGORITHMS && distinct; j++) {
                int same = 1;
                for (int s = 0; s < 6; s++) {
                    make(&a, i, seeds[s]);
                    make(&b, j, seeds[s]);
                    if (!identical(&a, &b))
                        same = 0;
                }
                if (same) {
                    sprintf(what, "%d and %d", i + 1, j + 1);
                    distinct = 0;
                }
            }
        check(distinct, "and no two of them are the same algorithm twice");
    }

    /* --- every palette is a palette --------------------------------------- */
    for (int alg = 0; alg < PALGORITHMS; alg++) {
        int filled = 1, ranged = 1, worst = 766;
        for (int s = 0; s < 6; s++) {
            make(&a, alg, seeds[s]);
            if (a.size < 64)
                filled = 0;
            int lo, hi;
            range(&a, &lo, &hi);
            if (hi - lo < worst)
                worst = hi - lo;
            /* a third of the way from black to white, which is far less than
             * any of the seven actually manages and far more than a wash */
            if (hi - lo < 255)
                ranged = 0;
        }
        sprintf(what, "algorithm %d fills its entries", alg + 1);
        check(filled, what);
        sprintf(what, "algorithm %d has dark and light in it (%d of 765)",
                alg + 1, worst);
        check(ranged, what);
    }

    /* --- and the four that are schemes have more than one colour in them ---
     *
     * A spectrum, two inks, three hues spread round the circle and two
     * opposite ones: none of those is one colour, and two of them were. The
     * first three of the seven are free to be whatever the dice say. */
    {
        static const struct {
            int alg;
            double least;
            const char *what;
        } scheme[4] = {{3, 0.55, "the spectrum goes most of the way round"},
                       {4, 0.10, "the two inks are two colours and not one"},
                       {5, 0.25, "the three spread hues are three"},
                       {6, 0.20, "and the opposite pair are opposite"}};
        for (int k = 0; k < 4; k++) {
            double worst = 2;
            for (int s = 0; s < 6; s++) {
                make(&a, scheme[k].alg, seeds[s]);
                double spread = hue_spread(&a);
                if (spread < worst)
                    worst = spread;
            }
            sprintf(what, "%s (%.2f of the circle)", scheme[k].what, worst);
            check(worst >= scheme[k].least, what);
        }
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    else
        printf("\nok     seven ways to make a palette, and all of them do\n");
    return failures != 0;
}
