/* The colour speed functions.
 *
 * The nine there were must go on giving what they gave, to the bit, since the
 * number of one is what a saved position holds. The three added must do what
 * they are for: log log and exp as the names say, triangle as a straightened
 * sine. The speed change after ITER, which acts on whichever of them is
 * chosen, must do nothing at all when off or at an exponent of one, and
 * otherwise what its two modes say -- jump untouched below ITER and its pace
 * turning there at once, smooth untouched below ITER less the radius, never
 * turning back, without a corner, and the jump moved on by a constant past
 * ITER plus the radius. A position from 1.8.0 in continuous or in its smooth
 * must come out near what it was. And a value past what an int holds has to
 * land in the palette cycle like any other, where it used to land on one
 * flat entry.
 *
 * formulas.cpp is included rather than linked, since the functions are its
 * own and static.
 */
#include <cmath>
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
static char what[160];

static void check(int ok, const char *text)
{
    printf("%s   %s\n", ok ? "ok  " : "FAIL", text);
    if (!ok)
        failures++;
}

/* A function alone: no speed, no shift, the outside's numbers. */
static double F(int func, double v)
{
    return color_f_func((number_t)v, 1.0f, func, 0, 0);
}

enum { OFF = 0, JUMP = 1, SMOOTH = 2 };

/* 1.8.0's jump -- one group, an initial exponent of one -- in the long double
 * it is worked in: slopes are taken from this, a float's rounding being
 * coarser than a small step. */
static long double P(long double v)
{
    return color_power(v, (long double)cfractalc.outcolorfuniter,
                       (long double)cfractalc.outcolorfunexp);
}

/* 1.8.0's smooth and continuous, as it drew them, for setting the modes
 * they are read into beside them: the corner rounded in the logarithm, with
 * the k that kept the count at half of ITER within two per cent; and the pace
 * going on from ITER unbroken. */
static double sharpness_180(double p)
{
    const double limit = log(1.02), d = fabs(p - 1);
    auto moved = [d](double k) { return d / k * log1p(exp2(-k)); };
    if (!(d > 0) || moved(4) <= limit)
        return 4;
    double lo = 4, hi = 8;
    while (moved(hi) > limit && hi < 4096)
        hi *= 2;
    for (int i = 0; i < 60; i++) {
        double mid = (lo + hi) / 2;
        (moved(mid) > limit ? lo : hi) = mid;
    }
    return hi;
}
static double smooth_180(double v, double I, double p)
{
    const double k = sharpness_180(p);
    return v * pow(1 + pow(v / I, k), (p - 1) / k);
}
static double continuous_180(double v, double I, double p)
{
    if (v <= I)
        return v;
    const double x = log(v / I);
    return fabs(p) < 1e-12 ? I * (1 + x) : I + I / p * expm1(p * x);
}

/* The change as the outside's plan works it, whatever its groups. */
static long double Q(long double v) { return color_change(v, cfractalc.outplan); }

/* The outside's initial exponent and both groups, in the order of the
 * dialog, and the plan made from them, as set_fractalc would. */
static void groups(double initial, double iter, double p, int mode, double r,
                   double iter2, double p2, int mode2, double r2)
{
    cfractalc.outcolorfuninitial = initial;
    cfractalc.outcolorfuniter = iter;
    cfractalc.outcolorfunexp = p;
    cfractalc.outcolorfunmode = mode;
    cfractalc.outcolorfunradius = r;
    cfractalc.outcolorfun2iter = iter2;
    cfractalc.outcolorfun2exp = p2;
    cfractalc.outcolorfun2mode = mode2;
    cfractalc.outcolorfun2radius = r2;
    color_change_prepare(&cfractalc);
}

/* One group: an initial exponent of one, the second group off. */
static void after(double iter, double p, int mode,
                  double r = COLORFUN_RADIUS)
{
    groups(1, iter, p, mode, r, COLORFUN2_ITER, COLORFUN2_EXP, OFF,
           COLORFUN2_RADIUS);
}

/* The exponent the count has at v: its log-log slope, from Q. */
static double local_exponent(long double v)
{
    long double h = v * 1e-7L;
    return (double)((logl(Q(v + h)) - logl(Q(v - h))) /
                    (logl(v + h) - logl(v - h)));
}

/* The slopes either side of v, from Q. */
static double slope_left(long double v)
{
    long double h = v * 1e-7L;
    return (double)((Q(v) - Q(v - h)) / h);
}
static double slope_right(long double v)
{
    long double h = v * 1e-7L;
    return (double)((Q(v + h) - Q(v)) / h);
}

static unsigned long long fnv(unsigned long long h, unsigned int x, int bytes)
{
    for (int b = 0; b < bytes; b++) {
        h ^= (x >> (8 * b)) & 0xFF;
        h *= 1099511628211ULL;
    }
    return h;
}

int main(void)
{
    const double pi = 3.14159265358979323846;

    /* --- the nine there were --------------------------------------------
     *
     * Checksums of what each returned before the four were added, over
     * values from minus to plus sixty-four thousand, at a speed of 1.5 and a
     * shift of three; the same at both precisions. */
    {
        static const unsigned long long was[9] = {
            0x86f94c9dad52da04ULL, /* ident */
            0x89d81777f3ee2d0bULL, /* log */
            0x1cfb0e8703868bafULL, /* atan */
            0xc89a928f6df9518cULL, /* sin */
            0x737ecbbb1ccbf2ffULL, /* cos */
            0x145fd807c436e286ULL, /* square */
            0x0bcaaa57389558bbULL, /* cube */
            0x831fbed1d3c2ee05ULL, /* square root */
            0xfa4bfa4d03db8e0aULL, /* cube root */
        };
        for (int f = 0; f < 9; f++) {
            unsigned long long h = 1469598103934665603ULL;
            for (int i = -4000; i <= 4000; i++) {
                number_t v = (number_t)i * i * i / 1000000 + (number_t)i / 7;
                float r = color_f_func(v, 1.5f, f, 3, 0);
                unsigned int bits;
                memcpy(&bits, &r, 4);
                h = fnv(h, bits, 4);
            }
            sprintf(what, "%s gives what it gave", colorfun[f]);
            check(h == was[f], what);
        }
    }

    /* And the palette entries a picture is drawn with, where the value stays
     * inside an int -- which is where the old conversion was right, and where
     * the wrap that replaced it has to be the same. */
    static unsigned int pix[257];
    for (int k = 0; k < 257; k++)
        pix[k] = k * 2654435761u;
    cpalette.pixels = pix;
    cpalette.size = 200;
    cpalette.maxentries = 256;
    cpalette.type = C256;
    cfractalc.maxiter = 100000;
    cfractalc.outcolorspeed = 1.5f;
    cfractalc.outcolorshift = 7;
    cfractalc.coloringmode = OutColormodeType::ColOut_iter_plus_real;
    {
        static const unsigned long long was[9] = {
            0x09162b3dabbd603fULL, 0x415aa9f6ce56c664ULL, 0xd58647bfc82ccb29ULL,
            0x24a39108ef3ac5d4ULL, 0x93ce72292ae57784ULL, 0x2997877a97e99914ULL,
            0xf82f56aeb277f090ULL, 0x53ea6974bedf79d8ULL, 0x5e9d70857b8e0f66ULL,
        };
        for (int f = 0; f < 9; f++) {
            cfractalc.outcolorfun = f;
            unsigned long long h = 1469598103934665603ULL;
            for (unsigned it = 0; it < 160; it++)
                for (int k = -3; k <= 3; k++) {
                    unsigned int r = color_output((number_t)k / 3, 0.25L, it);
                    h ^= r;
                    h *= 1099511628211ULL;
                }
            sprintf(what, "%s colours the counts as it did", colorfun[f]);
            check(h == was[f], what);
        }
    }

    /* --- past an int ---------------------------------------------------
     *
     * cube passes 2^31 256ths of a band at 204 iterations. From there every
     * pixel took the same entry; now each lands where its value says, the
     * value wrapped round the cycle as it is below 2^31. */
    {
        for (int k = 0; k < 257; k++)
            pix[k] = k;
        cpalette.size = 256;
        cfractalc.outcolorspeed = 1.0f;
        cfractalc.outcolorshift = 0;
        cfractalc.coloringmode = OutColormodeType::ColOut_iter;
        cfractalc.outcolorfun = 6;
        const unsigned int cycle = 255 << 8;
        int right = 1, seen[256] = {0}, distinct = 0;
        for (unsigned it = 204; it < 1000; it++) {
            float r = color_f_func((number_t)it, 1.0f, 6, 0, 0);
            unsigned int at = (unsigned int)fmodl(
                truncl((long double)r * 256), (long double)cycle);
            unsigned int got = color_output(0, 0, it);
            if (got != 1 + (at >> 8))
                right = 0;
            if (!seen[got & 255]++)
                distinct++;
        }
        sprintf(what, "cube past 204 iterations lands where its value says "
                      "(%d entries of 255 seen, not one)",
                distinct);
        check(right && distinct > 100, what);
        cfractalc.outcolorfun = 0;
    }

    /* --- log log, exp, triangle ----------------------------------------- */
    {
        int odd = 1, rising = 1;
        double last = -1e300;
        for (int i = -300; i <= 300; i++) {
            double v = i * 13.7;
            if (fabs(F(9, v) + F(9, -v)) > 1e-6)
                odd = 0;
            if (F(9, v) < last)
                rising = 0;
            last = F(9, v);
        }
        check(odd && rising, "log log is odd and rising");
        check(fabs(F(9, exp(1.0) - 1) - log(2.0)) < 1e-6,
              "and log(1 + log(1 + v)): ln 2 at e - 1");
        sprintf(what, "and 2.7 bands at a million iterations (%.3f)",
                F(9, 1e6));
        check(fabs(F(9, 1e6) - 2.70) < 0.01, what);
    }
    {
        int odd = 1;
        for (int i = -80; i <= 80; i++) {
            double v = i / 7.0;
            if (fabs(F(10, v) + F(10, -v)) > 1e-5 * (1 + fabs(F(10, v))))
                odd = 0;
        }
        check(odd, "exp is odd");
        check(fabs(F(10, 1) - (exp(1.0) - 1)) < 1e-6, "and e - 1 at one");
        double big = F(10, 5000);
        check(std::isfinite(big) && big > 3e38,
              "and a value past what a float holds is held at the largest "
              "one, not infinity");
    }
    {
        int odd = 1, period = 1, bounded = 1;
        for (int i = -400; i <= 400; i++) {
            double v = i / 9.0;
            double t = F(11, v);
            if (fabs(t + F(11, -v)) > 1e-6)
                odd = 0;
            if (fabs(t - F(11, v + 2 * pi)) > 1e-5)
                period = 0;
            if (fabs(t) > 1 + 1e-6)
                bounded = 0;
        }
        check(odd && period && bounded,
              "triangle is odd, repeats every two pi and stays within one");
        check(fabs(F(11, pi / 2) - 1) < 1e-6 && fabs(F(11, pi)) < 1e-6 &&
                  fabs(F(11, 3 * pi / 2) + 1) < 1e-6,
              "and turns where the sine turns");
        double slope = (F(11, 0.3) - F(11, 0.1)) / 0.2;
        sprintf(what, "and climbs at a steady 2/pi between the turns (%.6f)",
                slope);
        check(fabs(slope - 2 / pi) < 1e-5, what);
    }

    /* --- the speed change after ITER ------------------------------------
     *
     * Not a function of its own but a change to the value any function is
     * handed. Off it does nothing; and at an exponent of one it does nothing
     * in any mode -- to the bit, for every function, since that is what keeps
     * every picture drawn before it came as it was. */
    check(colorfun_modes[0] && !strcmp(colorfun_modes[0], "off") &&
              !strcmp(colorfun_modes[1], "jump") &&
              !strcmp(colorfun_modes[2], "smooth") && colorfun_modes[3] == NULL,
          "the speed change offers off, jump and smooth");
    check(!strcmp(colorfun_after_modes[2], "continuous") &&
              !strcmp(colorfun_after_modes[3], "smooth") &&
              colorfun_after_modes[4] == NULL,
          "and still reads the four 1.8.0 named");
    {
        int same = 1;
        for (int f = 0; f < COLORFUN - 1; f++)
            for (int mode = OFF; mode <= SMOOTH; mode++) {
                for (int i = -300; i <= 300; i++) {
                    number_t v = (number_t)i * i * i / 10000 + (number_t)i / 7;
                    after(50, 1, mode);
                    float a = color_f_func(v, 1.5f, f, 3, 0);
                    after(50, 2, OFF);
                    float b = color_f_func(v, 1.5f, f, 3, 0);
                    if (memcmp(&a, &b, sizeof a))
                        same = 0;
                }
            }
        check(same, "off, or at an exponent of one, every function gives to "
                    "the bit what it gives without it");
    }

    const double T = 50, R = 10;
    const double exps[] = {0.01, 0.5, 1.5, 2, 100};
    for (double p : exps) {
        /* jump: untouched up to ITER, and the pace turning to p there */
        after(T, p, JUMP);
        int kept = 1;
        for (int i = 1; i <= 100; i++) {
            double v = T * i / 100;
            if (F(0, v) != (float)v)
                kept = 0;
        }
        long double h = T * 1e-7L;
        double s1 = (double)((P(T + h) - P(T)) / h);
        sprintf(what,
                "exponent %g, jump: untouched up to ITER, the pace %.4g past "
                "it",
                p, s1);
        check(kept && fabs(s1 / p - 1) < 0.02, what);

        /* smooth: untouched up to ITER - R, no corner anywhere, the pace
         * always between the two, and past ITER + R the jump moved on by a
         * constant */
        after(T, p, SMOOTH, R);
        int before = 1, rising = 1;
        for (int i = 1; i <= 100; i++) {
            double v = (T - R) * i / 100;
            if (fabsl(Q(v) / v - 1) > 1e-15L)
                before = 0;
        }
        long double last = Q(T - R - 1);
        for (int i = 0; i <= 400; i++) {
            long double v = T - R - 1 + (2 * R + 2) * i / 400.0L, q = Q(v);
            if (!(q >= last))
                rising = 0;
            last = q;
        }
        int corners = 1;
        for (double at : {T - R, T, T + R}) {
            double l = slope_left(at), r = slope_right(at);
            if (fabs(l / r - 1) > 0.01)
                corners = 0;
        }
        long double moved = 0;
        int far = 1;
        for (int i = 0; i <= 100; i++) {
            long double v = T + R + 3 * T * i / 100.0L;
            after(T, p, JUMP);
            long double j = Q(v);
            after(T, p, SMOOTH, R);
            long double d = Q(v) - j;
            if (i == 0)
                moved = d;
            else if (fabsl(d - moved) > 1e-9L * fabsl(j) + 1e-9L)
                far = 0;
        }
        sprintf(what,
                "exponent %g, smooth: untouched up to ITER - R, rising, no "
                "corner, and the jump past ITER + R (moved by %.3g)",
                p, (double)moved);
        check(before && rising && corners && far, what);
        /* both odd, and finite however large */
        int odd = 1, finite = 1;
        for (int mode = JUMP; mode <= SMOOTH; mode++) {
            after(T, p, mode, R);
            for (int i = 1; i < 200; i++) {
                double v = i * 3.1;
                if (F(0, v) != -F(0, -v))
                    odd = 0;
            }
            if (!std::isfinite(F(0, 1e6)))
                finite = 0;
        }
        sprintf(what, "exponent %g: odd and finite in both", p);
        check(odd && finite, what);
    }
    {
        /* where both paces are steady, smooth moves nothing on: an exponent
         * of nought after the count untouched comes back to the jump exactly
         * past ITER + R, and at ITER is a quarter of the radius short of it,
         * the corner cut as a parabola cuts it */
        after(T, 0, SMOOTH, R);
        double at = (double)Q(T), past = (double)Q(T + 2 * R);
        sprintf(what,
                "steady paces: %.6g at ITER, a quarter of the radius under "
                "it, and %.6g past the rounding, ITER again",
                at, past);
        check(fabs(at - (T - R / 4)) < 1e-9 && fabs(past - T) < 1e-9, what);
        /* and a radius of nought is the jump */
        int same = 1;
        for (double p : {0.5, 2.0})
            for (int i = 1; i <= 300; i++) {
                long double v = i * 0.7L;
                after(T, p, JUMP);
                long double a = Q(v);
                after(T, p, SMOOTH, 0);
                if (fabsl(Q(v) / a - 1) > 1e-15L)
                    same = 0;
            }
        check(same, "smooth at a radius of nought is the jump");
    }
    {
        /* On a function other than ident: the square root up to ITER, and
         * past it the square root of the changed value -- the function
         * chosen, driven faster from there. */
        int before = 1, past = 1;
        after(T, 2, JUMP);
        for (int i = 1; i <= 300; i++) {
            double v = i * 0.5;
            double got = F(7, v);
            double want = v <= T ? sqrt(v) : sqrt((double)P(v));
            if (v <= T) {
                after(T, 2, OFF);
                if (got != F(7, v))
                    before = 0;
                after(T, 2, JUMP);
            } else if (fabs(got - want) > 1e-6 * want) {
                past = 0;
            }
        }
        check(before && past, "square root with the change: itself up to ITER, "
                              "the square root of the changed value past it");
        /* and a periodic one runs faster past ITER: sin of the changed value */
        int periodic = 1;
        after(T, 1.5, SMOOTH, R);
        for (int i = 1; i <= 300; i++) {
            double v = i * 0.5;
            if (fabs(F(3, v) - sin((double)Q(v))) > 1e-6)
                periodic = 0;
        }
        check(periodic, "and sin with it is the sine of the changed value");
    }
    {
        /* ITER at nought: a plain power from the start, in either mode */
        after(0, 2, JUMP);
        double a = F(0, 3);
        after(0, 0.5, SMOOTH, R);
        double b = F(0, 4);
        check(fabs(a - 9) < 1e-5 && fabs(b - 2) < 1e-6,
              "ITER at nought is a plain power: 3 squared is 9, 4 to the half "
              "is 2");
        /* a rounding wider than ITER stops short of nought, and is still
         * rising and finite there */
        after(10, 2, SMOOTH, 50);
        int fine = 1;
        long double last = 0;
        for (int i = 1; i <= 400; i++) {
            long double v = i * 0.1L, q = Q(v);
            if (!std::isfinite((double)q) || !(q > last))
                fine = 0;
            last = q;
        }
        check(fine, "a radius past ITER stops short of nought, rising all "
                    "the way");
    }

    /* --- 1.8.0's continuous and smooth ---------------------------------
     *
     * Read into these two near enough: drawn side by side over four ITERs,
     * the count keeps within a few hundredths of ITER of what 1.8.0 drew. */
    {
        double worst_s = 0, worst_c = 0;
        for (double p : {0.0, 0.01, 0.1, 0.5}) {
            number_t e = p, r = COLORFUN_RADIUS;
            int m = 3;
            colorfun_from_legacy(100, &e, &m, &r);
            groups(1, 100, (double)e, m, (double)r, COLORFUN2_ITER,
                   COLORFUN2_EXP, OFF, COLORFUN2_RADIUS);
            for (int i = 1; i <= 400; i++) {
                double v = i, gap = fabs((double)Q(v) - smooth_180(v, 100, p));
                if (gap > worst_s)
                    worst_s = gap;
            }
            e = p, r = COLORFUN_RADIUS, m = 2;
            colorfun_from_legacy(100, &e, &m, &r);
            groups(1, 100, (double)e, m, (double)r, COLORFUN2_ITER,
                   COLORFUN2_EXP, OFF, COLORFUN2_RADIUS);
            for (int i = 1; i <= 400; i++) {
                double v = i,
                       gap = fabs((double)Q(v) - continuous_180(v, 100, p));
                if (gap > worst_c)
                    worst_c = gap;
            }
        }
        sprintf(what,
                "1.8.0's smooth and continuous, exponents up to a half, read "
                "into smooth: within %.3g and %.3g of ITER of it",
                worst_s / 100, worst_c / 100);
        check(worst_s < 5 && worst_c < 11, what);
        /* off and jump as they were */
        number_t e = 1.5, r = 7;
        int m = 1;
        colorfun_from_legacy(100, &e, &m, &r);
        int m0 = 0;
        number_t e0 = 2, r0 = 7;
        colorfun_from_legacy(100, &e0, &m0, &r0);
        check(m == JUMP && e == 1.5 && r == 7 && m0 == OFF && e0 == 2,
              "1.8.0's off and jump read as they were");
    }

    /* --- an initial exponent, and a second group -------------------------
     *
     * One group in jump with an initial exponent of one is still 1.8.0's
     * change, worked by color_power. The stretches that work everything else
     * have to agree with it there: a second group at ITER nought and exponent
     * one is the count untouched from the start, which sends the first
     * through them. */
    {
        int agree = 1;
        double worst = 0;
        for (double p : {0.01, 0.5, 1.5, 2.0, 100.0}) {
            after(T, p, JUMP);
            long double want[400];
            for (int i = 0; i < 400; i++)
                want[i] = P((long double)(i + 1) * T / 100);
            groups(1, T, p, JUMP, R, 0, 1, JUMP, R);
            if (cfractalc.outplan.kind != 2)
                agree = 0;
            for (int i = 0; i < 400; i++) {
                long double v = (long double)(i + 1) * T / 100;
                long double got = Q(v);
                if (!std::isfinite((double)want[i]))
                    continue;
                double off = (double)fabsl(got / want[i] - 1);
                if (off > worst)
                    worst = off;
                if (!(off < 1e-12) || Q(-v) != -got)
                    agree = 0;
            }
        }
        sprintf(what,
                "the stretches give 1.8.0's jump where it is the same "
                "(worst %.1e of it)",
                worst);
        check(agree, what);
    }
    {
        /* a group that is off changes nothing, whatever it holds -- to the
         * bit; and with both off nothing is done, whatever the initial
         * exponent */
        int same = 1;
        for (int mode = JUMP; mode <= SMOOTH; mode++)
            for (int i = 1; i <= 300; i++) {
                double v = i * 0.7;
                after(T, 1.7, mode, R);
                float a = (float)F(0, v);
                groups(1, T, 1.7, mode, R, 30, 5, OFF, 3);
                float b = (float)F(0, v);
                groups(1, 80, 9, OFF, 4, T, 1.7, mode, R);
                float c = (float)F(0, v);
                if (memcmp(&a, &b, sizeof a) || memcmp(&a, &c, sizeof a))
                    same = 0;
            }
        check(same, "a group that is off changes nothing, to the bit");
        int nothing = 1;
        for (int i = -300; i <= 300; i++) {
            double v = i * 0.7;
            after(T, 2, OFF);
            float a = (float)F(5, v);
            groups(1, T, 3, OFF, R, 20, 4, OFF, R);
            float b = (float)F(5, v);
            if (memcmp(&a, &b, sizeof a))
                nothing = 0;
        }
        check(nothing && cfractalc.outplan.kind == 0,
              "both groups off and the initial exponent at one: nothing at "
              "all, to the bit");
    }
    {
        /* A group that is off is left out, ITER and all -- the three cases
         * as they were asked for, with the initial exponent at 0.5, the first
         * group at 30 and the second at 60, jump. Up to the ITER of the one
         * that is on, the plain power; past it, what that group adds. */
        groups(0.5, 30, 5, OFF, R, 60, 2, JUMP, R);
        double a1 = local_exponent(20), a2 = local_exponent(45),
               a3 = (double)Q(100);
        const double a3want = sqrt(60.0) + 60 * (pow(100 / 60.0, 2) - 1);
        sprintf(what,
                "first off, second on: v^0.5 up to the second's ITER (%.6g at "
                "20, %.6g at 45), then its 2 from there (100 to %.6g)",
                a1, a2, a3);
        check(fabs(a1 - 0.5) < 1e-6 && fabs(a2 - 0.5) < 1e-6 &&
                  fabs(a3 / a3want - 1) < 1e-12,
              what);
        groups(0.5, 30, 3, JUMP, R, 60, 7, OFF, R);
        double b1 = local_exponent(20), b2 = (double)Q(45),
               b3 = (double)Q(100);
        const double b2want = sqrt(30.0) + 30 * (pow(1.5, 3) - 1),
                     b3want = sqrt(30.0) + 30 * (pow(100 / 30.0, 3) - 1);
        sprintf(what,
                "first on, second off: v^0.5 up to the first's ITER (%.6g), "
                "then its 3 on and on (45 to %.6g, 100 to %.6g)",
                b1, b2, b3);
        check(fabs(b1 - 0.5) < 1e-6 && fabs(b2 / b2want - 1) < 1e-12 &&
                  fabs(b3 / b3want - 1) < 1e-12,
              what);
        groups(0.5, 30, 3, OFF, R, 60, 7, OFF, R);
        double c1 = local_exponent(20), c2 = local_exponent(45),
               c3 = local_exponent(100);
        double c100 = (double)Q(100);
        sprintf(what,
                "both off: 0.5 for every count (%.6g, %.6g, %.6g), the plain "
                "power -- 100 to %.6g",
                c1, c2, c3, c100);
        check(fabs(c1 - 0.5) < 1e-6 && fabs(c2 - 0.5) < 1e-6 &&
                  fabs(c3 - 0.5) < 1e-6 && fabs(c100 - 10) < 1e-9 &&
                  Q(-100) == -Q(100),
              what);
    }

    const double initials[] = {0.5, 2, 3};
    for (double b : initials) {
        /* up to ITER the count is the plain power v^initial in jump, and up
         * to ITER - R in smooth, wherever ITER is -- as with both groups off */
        int below = 1;
        for (int mode = JUMP; mode <= SMOOTH; mode++) {
            groups(b, T, 1.5, mode, R, COLORFUN2_ITER, COLORFUN2_EXP, OFF, R);
            const double upto = mode == JUMP ? T : T - R;
            for (int i = 1; i <= 100; i++) {
                double v = upto * i / 100;
                if (fabs((double)Q(v) / pow(v, b) - 1) > 1e-12)
                    below = 0;
            }
        }
        sprintf(what, "initial %g: up to ITER the count is v^%g", b, b);
        check(below, what);
        /* the pace coming to ITER is the power's; past it jump turns to the
         * group's exponent at once, and smooth has no corner anywhere */
        const double come = b * pow(T, b - 1);
        groups(b, T, 1.5, JUMP, R, COLORFUN2_ITER, COLORFUN2_EXP, OFF, R);
        double jl = slope_left(T), jr = slope_right(T);
        groups(b, T, 1.5, SMOOTH, R, COLORFUN2_ITER, COLORFUN2_EXP, OFF, R);
        int corners = 1;
        for (double at : {T - R, T, T + R})
            if (fabs(slope_left(at) / slope_right(at) - 1) > 0.01)
                corners = 0;
        sprintf(what,
                "initial %g: jump turns from %.4g to %.4g at ITER, smooth has "
                "no corner",
                b, jl, jr);
        check(fabs(jl / come - 1) < 0.01 && fabs(jr / 1.5 - 1) < 0.01 &&
                  corners,
              what);
    }

    {
        /* Each stretch has its own exponent and no other: v^p0 up to the
         * first ITER, then from each ITER what its group adds, at the pace
         * of its exponent there -- not p0 p1, not p0 p1 p2. */
        groups(2, 50, 3, JUMP, R, 100, 0.5, JUMP, R);
        double e0 = local_exponent(25), s1 = slope_right(50),
               s2 = slope_right(100);
        double q75 = (double)Q(75), q400 = (double)Q(400);
        sprintf(what,
                "each stretch its own exponent: v^2 up to 50 (%.6g), then 3 "
                "from 50 (pace %.5g, 75 to %.6g), then 0.5 from 100 (pace "
                "%.5g, 400 to %.6g)",
                e0, s1, q75, s2, q400);
        check(fabs(e0 - 2) < 1e-6 && fabs(s1 / 3 - 1) < 1e-4 &&
                  fabs(q75 - (2500 + 50 * (pow(1.5, 3) - 1))) < 1e-9 &&
                  fabs(s2 / 0.5 - 1) < 1e-4 && fabs(q400 - 2950) < 1e-9,
              what);

        /* what the stretches were asked to be: squared up to ITER 50; an
         * exponent of one from there, so that the count goes on as v + a, a
         * being how far past 50 the square had taken it; squared again from
         * ITER 100 */
        {
            const double a = 2500 - 50;
            groups(2, 50, 1, JUMP, R, 100, 2, JUMP, R);
            double u30 = (double)Q(30), u80 = (double)Q(80),
                   p80 = slope_right(80), u150 = (double)Q(150),
                   p100 = slope_right(100);
            sprintf(what,
                    "v^2 up to 50 (30 to %.6g), v + %g from there (80 to "
                    "%.6g, pace %.5g), squared again from 100 (pace %.5g, "
                    "150 to %.6g)",
                    u30, a, u80, p80, p100, u150);
            check(fabs(u30 - 900) < 1e-9 && fabs(u80 - (80 + a)) < 1e-9 &&
                      fabs(p80 - 1) < 1e-4 && fabs(p100 / 2 - 1) < 1e-4 &&
                      fabs(u150 - (100 + a + 100 * (1.5 * 1.5 - 1))) < 1e-9,
                  what);
            /* smooth goes on as v + a too once past its radius, at a pace of
             * one -- a moved on by what the rounding took off the square */
            groups(2, 50, 1, SMOOTH, R, 100, 2, OFF, R);
            double p70 = slope_right(70), d70 = (double)Q(70) - 70,
                   d90 = (double)Q(90) - 90;
            sprintf(what,
                    "smooth after a square goes on as v + %.6g past its "
                    "radius (pace %.5g)",
                    d70, p70);
            check(fabs(p70 - 1) < 1e-4 && fabs(d70 - d90) < 1e-9 &&
                      fabs(d70 - a) < R * R,
                  what);
        }

        /* the example the change was asked for first: the second group at
         * ITER 50, squaring from there; the first at 100, one and a half
         * from there */
        groups(1, 100, 1.5, JUMP, R, 50, 2, JUMP, R);
        double a = (double)Q(30), b = (double)Q(75), c = (double)Q(150);
        sprintf(what,
                "two groups, a stretch each: 30 stays 30 (%.4g), 75 is "
                "50 x 1.5^2 (%.4g), 150 is 200 + 100 (1.5^1.5 - 1) (%.5g)",
                a, b, c);
        check(fabs(a - 30) < 1e-9 && fabs(b - 112.5) < 1e-9 &&
                  fabs(c - (200 + 100 * (pow(1.5, 1.5) - 1))) < 1e-9,
              what);
        /* the lower ITER goes first whichever group holds it, and the
         * initial exponent runs up to it */
        int swapped = 1;
        for (int mode = JUMP; mode <= SMOOTH; mode++)
            for (int i = 1; i <= 300; i++) {
                groups(0.5, 100, 1.5, mode, 20, 50, 2, mode, 5);
                long double x = Q(i * 0.9);
                groups(0.5, 50, 2, mode, 5, 100, 1.5, mode, 20);
                if (fabsl(Q(i * 0.9) / x - 1) > 1e-15L)
                    swapped = 0;
            }
        groups(0.5, 100, 1.5, JUMP, R, 50, 2, JUMP, R);
        double i25 = (double)Q(25);
        sprintf(what,
                "the lower ITER goes first, whichever group holds it, with "
                "the initial exponent up to it (25 to %.6g)",
                i25);
        check(swapped && fabs(i25 - 5) < 1e-9, what);
        /* at the same ITER the second replaces the first */
        groups(1, 100, 3, JUMP, R, 100, 2, JUMP, R);
        double r = (double)Q(200);
        sprintf(what, "at the same ITER the second replaces the first (%.4g)",
                r);
        check(fabs(r - 400) < 1e-9, what);
        /* a group at ITER nought: a plain power from the start, then the
         * other's corner */
        groups(1, 100, 1.5, JUMP, R, 0, 2, JUMP, R);
        double p1 = (double)Q(10), p2 = (double)Q(400);
        sprintf(what,
                "a group at ITER nought squares from the start (10 to %.4g) "
                "and the other turns at its ITER (400 to %.6g)",
                p1, p2);
        check(fabs(p1 - 100) < 1e-9 &&
                  fabs(p2 - (10000 + 100 * (pow(4, 1.5) - 1))) < 1e-6,
              what);
        /* smooth at both: no corner at either ITER or either end of a
         * rounding, and past both the jump moved on by a constant */
        groups(1, 100, 1.5, SMOOTH, 15, 50, 2, SMOOTH, 8);
        int corners = 1;
        for (double at : {42.0, 50.0, 58.0, 85.0, 100.0, 115.0})
            if (fabs(slope_left(at) / slope_right(at) - 1) > 0.01)
                corners = 0;
        long double s1500 = Q(1500), s2000 = Q(2000);
        groups(1, 100, 1.5, JUMP, 15, 50, 2, JUMP, 8);
        long double moved = s1500 - Q(1500);
        double drift = (double)fabsl(s2000 - Q(2000) - moved);
        sprintf(what,
                "smooth at both: no corner, and the jump past both moved on "
                "by %.4g (%.1e off it)",
                (double)moved, drift);
        check(corners && drift < 1e-9, what);
        /* two roundings that would overlap are cut down alike, and meet
         * without a corner */
        groups(1, 50, 2, SMOOTH, 40, 80, 0.5, SMOOTH, 40);
        double meet_l = slope_left(65), meet_r = slope_right(65);
        double untouched = (double)Q(34);
        sprintf(what,
                "radii of 40 and 40 between ITERs 30 apart cut to 15 each: "
                "34 untouched (%.6g), no corner where they meet (%.4g, %.4g)",
                untouched, meet_l, meet_r);
        check(fabs(untouched - 34) < 1e-12 && fabs(meet_l / meet_r - 1) < 0.01,
              what);
        /* every mixture odd and finite */
        int odd = 1, finite = 1;
        for (int m1 = JUMP; m1 <= SMOOTH; m1++)
            for (int m2 = JUMP; m2 <= SMOOTH; m2++) {
                groups(0.5, 100, 2, m1, 30, 40, 0.7, m2, 30);
                for (int i = 1; i < 200; i++) {
                    double v = i * 3.1;
                    if (F(0, v) != -F(0, -v))
                        odd = 0;
                }
                if (!std::isfinite(F(0, 1e6)))
                    finite = 0;
            }
        check(odd && finite, "two groups in every mixture of modes: odd and "
                             "finite");
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    return failures != 0;
}
