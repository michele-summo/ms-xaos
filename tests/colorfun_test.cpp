/* The colour speed functions.
 *
 * The nine there were must go on giving what they gave, to the bit, since the
 * number of one is what a saved position holds. The three added must do what
 * they are for: log log and exp as the names say, triangle as a straightened
 * sine. The speed change after ITER, which acts on whichever of them is
 * chosen, must do nothing at all when off or at an exponent of one, and
 * otherwise what its three modes say -- untouched below ITER in jump and
 * continuous, the slope jumping in the first and continuous in the second,
 * smooth rounded and still within two per cent at half of ITER whatever the
 * exponent. And a value past what an int holds has to land in the palette
 * cycle like any other, where it used to land on one flat entry.
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

/* The change itself, in the long double it is worked in: slopes are taken
 * from this, a float's rounding being coarser than a small step. */
static long double P(long double v)
{
    return color_power(v, (long double)cfractalc.outcolorfuniter,
                       (long double)cfractalc.outcolorfunexp,
                       cfractalc.outcolorfunmode,
                       (long double)cfractalc.outcolorfunsharp);
}

static void after(double iter, double p, int mode)
{
    cfractalc.outcolorfuniter = iter;
    cfractalc.outcolorfunexp = p;
    cfractalc.outcolorfunmode = mode;
    cfractalc.outcolorfunsharp = color_power_sharpness(p);
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
    check(colorfun_after_modes[0] && !strcmp(colorfun_after_modes[0], "off") &&
              !strcmp(colorfun_after_modes[1], "jump") &&
              !strcmp(colorfun_after_modes[2], "continuous") &&
              !strcmp(colorfun_after_modes[3], "smooth") &&
              colorfun_after_modes[4] == NULL,
          "the speed change offers off, jump, continuous and smooth");
    {
        int same = 1;
        for (int f = 0; f < COLORFUN - 1; f++)
            for (int mode = 0; mode <= 3; mode++) {
                after(50, 1, mode);
                for (int i = -300; i <= 300; i++) {
                    number_t v = (number_t)i * i * i / 10000 + (number_t)i / 7;
                    after(50, 1, mode);
                    float a = color_f_func(v, 1.5f, f, 3, 0);
                    after(50, 2, 0);
                    float b = color_f_func(v, 1.5f, f, 3, 0);
                    if (memcmp(&a, &b, sizeof a))
                        same = 0;
                }
            }
        check(same, "off, or at an exponent of one, every function gives to "
                    "the bit what it gives without it");
    }

    const double T = 50;
    const double exps[] = {0.01, 0.5, 1.5, 2, 100};
    for (double p : exps) {
        /* jump and continuous leave everything up to ITER untouched */
        for (int mode = 1; mode <= 2; mode++) {
            after(T, p, mode);
            int kept = 1;
            for (int i = 1; i <= 100; i++) {
                double v = T * i / 100;
                if (F(0, v) != (float)v)
                    kept = 0;
            }
            sprintf(what, "exponent %g, %s: untouched up to ITER", p,
                    colorfun_after_modes[mode]);
            check(kept, what);
        }
        /* the slope past ITER: p at once in jump, one in continuous */
        long double h = T * 1e-7L;
        after(T, p, 1);
        double s1 = (double)((P(T + h) - P(T)) / h);
        after(T, p, 2);
        double s2 = (double)((P(T + h) - P(T)) / h);
        sprintf(what,
                "exponent %g: past ITER the slope jumps to it in jump (%.4g) "
                "and stays one in continuous (%.4g)",
                p, s1, s2);
        check(fabs(s1 / p - 1) < 0.02 && fabs(s2 - 1) < 0.02, what);
        /* smooth: within two per cent at half of ITER, no corner at ITER */
        after(T, p, 3);
        double half = F(0, T / 2);
        double left = (double)((P(T) - P(T - h)) / h);
        double right = (double)((P(T + h) - P(T)) / h);
        sprintf(what,
                "exponent %g, smooth: %.2f%% off at half of ITER, and no "
                "corner at ITER (%.4g against %.4g)",
                p, 100 * (half - T / 2) / (T / 2), left, right);
        check(fabs(half - T / 2) / (T / 2) <= 0.0201 &&
                  fabs(left / right - 1) < 0.01,
              what);
        /* every mode odd, and finite however large */
        int odd = 1, finite = 1;
        for (int mode = 1; mode <= 3; mode++) {
            after(T, p, mode);
            for (int i = 1; i < 200; i++) {
                double v = i * 3.1;
                if (F(0, v) != -F(0, -v))
                    odd = 0;
            }
            if (!std::isfinite(F(0, 1e6)))
                finite = 0;
        }
        sprintf(what, "exponent %g: odd and finite in every mode", p);
        check(odd && finite, what);
    }
    {
        /* smooth is jump rounded: well past ITER the two agree */
        int agree = 1;
        for (double p : {0.5, 1.5, 2.0}) {
            after(T, p, 1);
            double a = F(0, 4 * T);
            after(T, p, 3);
            double b = F(0, 4 * T);
            if (fabs(b / a - 1) > 0.01)
                agree = 0;
        }
        check(agree, "smooth comes to jump well past ITER");
    }
    {
        /* On a function other than ident: the square root up to ITER, and
         * past it the square root of the changed value -- the function
         * chosen, driven faster from there. */
        int before = 1, past = 1;
        after(T, 2, 2);
        for (int i = 1; i <= 300; i++) {
            double v = i * 0.5;
            double got = F(7, v);
            double want = v <= T ? sqrt(v) : sqrt((double)P(v));
            if (v <= T) {
                after(T, 2, 0);
                if (got != F(7, v))
                    before = 0;
                after(T, 2, 2);
            } else if (fabs(got - want) > 1e-6 * want) {
                past = 0;
            }
        }
        check(before && past, "square root with the change: itself up to ITER, "
                              "the square root of the changed value past it");
        /* and a periodic one runs faster past ITER: sin of the changed value */
        int periodic = 1;
        after(T, 1.5, 3);
        for (int i = 1; i <= 300; i++) {
            double v = i * 0.5;
            if (fabs(F(3, v) - sin((double)P(v))) > 1e-6)
                periodic = 0;
        }
        check(periodic, "and sin with it is the sine of the changed value");
    }
    {
        /* ITER at nought: a plain power from the start */
        after(0, 2, 2);
        double a = F(0, 3);
        after(0, 0.5, 1);
        double b = F(0, 4);
        check(fabs(a - 9) < 1e-5 && fabs(b - 2) < 1e-6,
              "ITER at nought is a plain power: 3 squared is 9, 4 to the half "
              "is 2");
        /* continuous at an exponent of nought: its limit, a logarithm */
        after(T, 0, 2);
        double c = F(0, T * exp(1.0));
        check(fabs(c - 2 * T) < 1e-4,
              "continuous at an exponent of nought goes on as a logarithm");
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    return failures != 0;
}
