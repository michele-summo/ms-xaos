/* The ways a palette can be made.
 *
 * Three of them are XaoS's own, colours between black and white anchors; four
 * more keep that skeleton and take their colours from elsewhere. What has to
 * hold of all seven is the same:
 *
 *  - a palette is made from a seed, and the same algorithm and seed must give
 *    the same palette, because that pair is all a saved position records of
 *    its colours -- get it wrong and a position comes back in the wrong ones;
 *  - a palette must actually use its entries, and no two algorithms may be
 *    the same algorithm under two numbers;
 *  - and every one of them must have some dark and some light in it, and more
 *    than one colour, in the part of it a picture actually shows. That last is
 *    the one that was missed: the palettes were made here two hundred and
 *    fifty-six entries long, where the program makes them sixty-five thousand
 *    long, and four ways that spread their colours over the whole length came
 *    out one colour each on the screen and several each here.
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

/* Whether a colour is one of the 256 the classic ring holds: a byte and the
 * two after it on the bottom byte of the generator, bright enough to be a
 * colour stop rather than a blend towards black. */
static int classic_ring(const unsigned char *c)
{
    int g = (109 * c[0] + 57) & 255;
    int b = (109 * g + 57) & 255;
    int mx = c[0] > c[1] ? (c[0] > c[2] ? c[0] : c[2]) : (c[1] > c[2] ? c[1] : c[2]);
    return mx > 40 && c[1] == g && c[2] == b;
}
/* One palette as the program makes it: truecolor, sixty-five thousand entries,
 * made after the default one -- which is what sets the length the next one is
 * laid out over, some three thousand segments of it. Only the first 4096
 * entries are kept, which is several times what a picture uses. */
static void make_app(struct made *out, int algorithm, int seed)
{
    struct palette pal;
    memset(&pal, 0, sizeof(pal));
    pal.type = TRUECOLOR;
    pal.maxentries = 65536;
    pal.end = 65535;
    pal.pixels = pixels;
    pal.alloccolor = alloccolor;
    memset(out, 0, sizeof(*out));
    building = out;
    mkdefaultpalette(&pal);
    mkpalette(&pal, seed, algorithm);
    out->size = pal.size < 4096 ? pal.size : 4096;
}

/* What the first n entries hold: the range of their brightness, as a sum of
 * the three channels; how many twelfths of the hue circle hold a twentieth or
 * more of the entries with a hue worth the name; and whether any of them is
 * warm (red to yellow) and any cool (cyan to blue). */
static void shown(const struct made *m, int n, int *range, int *families,
                  int *warm, int *cool)
{
    int lo = 766, hi = -1, bins[12] = {0}, hued = 0;
    *warm = *cool = 0;
    if (n > m->size)
        n = m->size;
    for (int i = 0; i < n; i++) {
        int r = m->rgb[i][0], g = m->rgb[i][1], b = m->rgb[i][2];
        int v = r + g + b;
        if (v < lo)
            lo = v;
        if (v > hi)
            hi = v;
        int mx = r > g ? (r > b ? r : b) : (g > b ? g : b);
        int mn = r < g ? (r < b ? r : b) : (g < b ? g : b);
        if (mx < 50 || mx - mn < mx / 4)
            continue;
        double d = mx - mn, h;
        if (mx == r)
            h = (g - b) / d;
        else if (mx == g)
            h = 2 + (b - r) / d;
        else
            h = 4 + (r - g) / d;
        h *= 60;
        if (h < 0)
            h += 360;
        bins[(int)(h / 30) % 12]++;
        hued++;
        if (h < 60 || h >= 330)
            *warm = 1;
        if (h >= 180 && h < 250)
            *cool = 1;
    }
    *range = hi - lo;
    *families = 0;
    for (int k = 0; k < 12; k++)
        if (hued && bins[k] * 20 >= hued)
            (*families)++;
}

static int identical(const struct made *a, const struct made *b)
{
    return a->size == b->size &&
           !memcmp(a->rgb, b->rgb, (size_t)a->size * 3);
}

/* How far the palette travels in brightness against how far it reaches: two if
 * it swells once from dark to light and back, and one more for every band it
 * turns at after that.
 *
 * This is what "the gradient is too slow" measures. A palette that swells once
 * across its whole length changes colour over a hundred entries where one that
 * alternates band by band changes over thirty, and a fractal drawn in the first
 * has no edges to its rings. Two of the four new ways swelled once and were
 * told so in as many words.
 *
 * Absolute numbers are no use here: how many bands there are to turn at is
 * chosen before any of the seven is asked, and at four segments not one of them
 * can do better than two. What can be asked is that the new ways are no slower
 * than the three that were always here, seed for seed. */
static double variation(const struct made *m)
{
    long travel = 0;
    int lo = 766, hi = -1;
    for (int i = 0; i < m->size; i++) {
        int v = m->rgb[i][0] + m->rgb[i][1] + m->rgb[i][2];
        if (v < lo)
            lo = v;
        if (v > hi)
            hi = v;
        if (i) {
            int u = m->rgb[i - 1][0] + m->rgb[i - 1][1] + m->rgb[i - 1][2];
            travel += v > u ? v - u : u - v;
        }
    }
    return hi > lo ? (double)travel / (hi - lo) : 0;
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
        int filled = 1;
        for (int s = 0; s < 6; s++) {
            make(&a, alg, seeds[s]);
            if (a.size < 64)
                filled = 0;
        }
        sprintf(what, "algorithm %d fills its entries", alg + 1);
        check(filled, what);
    }

    /* --- short palettes, as short as the older ones make them -------------
     *
     * Made two hundred and fifty-six entries long, a palette is three stops or
     * fewer a third of the time -- that is how mkpalette sizes its segments
     * when it is first asked, and it always has -- and three stops, the last
     * of them the first again, can miss a light or a dark whatever the
     * algorithm. The three older ones do, one palette in fourteen for 1 and 2
     * and one in three for 3. This used to ask all seven for a third of the
     * way from black to white on six seeds, which the older three happened to
     * pass; over a thousand seeds, what can fairly be asked of the newer ones
     * is that none of them does worse than the worst of the older three, give
     * or take five in a hundred.
     *
     * Not of 7, which is colours at random and nothing else, as it was asked
     * to be: it used to hold black and white every third stop, and that rhythm
     * of very dark and very light was exactly what made it not random. Random
     * colours make no promise of a dark and a light in three stops. */
    {
        int good[PALGORITHMS];
        for (int alg = 0; alg < PALGORITHMS; alg++) {
            good[alg] = 0;
            for (int sd = 0; sd < 1000; sd++) {
                make(&a, alg, 1 + sd * 7919);
                int lo, hi;
                range(&a, &lo, &hi);
                if (hi - lo >= 255)
                    good[alg]++;
            }
        }
        int floor = good[0];
        for (int alg = 1; alg < 3; alg++)
            if (good[alg] < floor)
                floor = good[alg];
        for (int alg = 3; alg < PALGORITHMS - 1; alg++) {
            sprintf(what,
                    "algorithm %d has dark and light as often as the older ones "
                    "(%d of 1000, the worst of them %d)",
                    alg + 1, good[alg], floor);
            check(good[alg] >= floor - 50, what);
        }
    }

    /* --- and none of them fades slowly from one end to the other -----------
     *
     * The three that were always here alternate light and dark segment by
     * segment, and that is where their banding comes from. Two of the four new
     * ones swelled once across the whole palette instead and were told, in as
     * many words, that the gradient was too slow. 7 again excepted: random
     * colours turn where the dice turn them. */
    for (int alg = 3; alg < PALGORITHMS - 1; alg++) {
        double worst = 1e9;
        for (int s = 0; s < 6; s++) {
            /* what the three that were always here manage on this seed */
            double old = 1e9;
            for (int k = 0; k < 3; k++) {
                make(&a, k, seeds[s]);
                double v = variation(&a);
                if (v < old)
                    old = v;
            }
            make(&a, alg, seeds[s]);
            double ratio = old > 0 ? variation(&a) / old : 1;
            if (ratio < worst)
                worst = ratio;
        }
        sprintf(what,
                "algorithm %d turns as often as the older ones do (%.2f of "
                "them)",
                alg + 1, worst);
        check(worst >= 0.95, what);
    }

    /* --- and in what a picture shows, they are palettes ---------------------
     *
     * The first six hundred entries of a palette made as the program makes it,
     * which is more than most pictures ever reach. Every one of the seven must
     * go from dark to light there and hold more than one colour: the three
     * that were always here do, on every seed tried, and the four that came
     * after them were one colour each until they were rebuilt. Warm over night
     * must have both of its windows showing. */
    {
        int worst_range[PALGORITHMS], fewest[PALGORITHMS], nightearth = 1;
        for (int alg = 0; alg < PALGORITHMS; alg++) {
            worst_range[alg] = 9999;
            fewest[alg] = 99;
        }
        for (int alg = 0; alg < PALGORITHMS; alg++)
            for (int sd = 0; sd < 60; sd++) {
                make_app(&a, alg, 1 + sd * 7919);
                int r, f, warm, cool;
                shown(&a, 600, &r, &f, &warm, &cool);
                if (r < worst_range[alg])
                    worst_range[alg] = r;
                if (f < fewest[alg])
                    fewest[alg] = f;
                if (alg == 4 && !(warm && cool))
                    nightearth = 0;
            }
        for (int alg = 0; alg < PALGORITHMS; alg++) {
            sprintf(what,
                    "algorithm %d shows dark and light where a picture looks "
                    "(%d of 765) and more than one colour (%d families)",
                    alg + 1, worst_range[alg], fewest[alg]);
            check(worst_range[alg] >= 250 && fewest[alg] >= 2, what);
        }
        check(nightearth, "and warm over night shows both warm and night");
    }

    /* --- smog is not the classic ring -------------------------------------
     *
     * 1 and 2 hold 256 colours, a byte and the two after it on the generator's
     * bottom byte. 4 walks a ring too, of hues, and must not come back to
     * theirs. What is read here is the entries rather than the stops, and an
     * entry between two stops is a blend that lands on one of those 256 now and
     * then by chance, where 2 shows two thousand in thirty palettes -- so what is
     * asked is that 4 shows a hundredth of what 2 does at most. */
    {
        int shared = 0, seen = 0;
        for (int sd = 0; sd < 30; sd++) {
            make_app(&a, 3, 1 + sd * 7919);
            make_app(&b, 1, 1 + sd * 7919);
            for (int i = 0; i < a.size; i++)
                shared += classic_ring(a.rgb[i]);
            for (int i = 0; i < b.size; i++)
                seen += classic_ring(b.rgb[i]);
        }
        sprintf(what, "smog is not the classic ring (%d entries on it, "
                      "against %d in 2)", shared, seen);
        check(seen > 0 && shared * 100 <= seen, what);
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    else
        printf("\nok     seven ways to make a palette, and all of them do\n");
    return failures != 0;
}
