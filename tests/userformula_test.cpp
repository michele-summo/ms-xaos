/* The user formula, driven the way the renderer drives it.
 *
 * Every other test of the parser stops at the parser. This one sets a formula
 * and an initialization on a context, hands them to the thread-local parsers
 * the way sffe_setlocal does, and iterates a grid through the loop the picture
 * is actually drawn by -- which is where the variables a formula may name are
 * either registered or not, and where an initialization that names one of the
 * missing ones is quietly thrown away instead of failing.
 *
 * Pictures are compared as checksums of the iteration counts: two formulas
 * that must mean the same thing must give the same one, and two that must
 * differ must not.
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "config.h"
#include "filter.h"
#include "fractal.h"
#include "xthread.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "phist.h"

/* --- what formulas.cpp needs from the rest of XaoS ------------------------ */
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

/* -------------------------------------------------------------------------- */

static int failures = 0;

static void check(int ok, const char *what)
{
    printf("%-6s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        failures++;
}

static const struct formula *userformula(void)
{
    for (const struct formula *f = formulas; f->magic == FORMULAMAGIC; f++)
        if (f->flags & SFFE_FRACTAL)
            return f;
    return NULL;
}

/* One picture, as a checksum of its iteration counts.
 *
 * The user formula carries no STARTZERO, so in mandelbrot mode calculate()
 * hands the loop the pixel as both the starting z and the point -- which is
 * why x, the starting z, is the plain coordinate in either mode while c is the
 * point in one and the constant in the other. Zero says the formula would not
 * parse, which no check below expects.
 */
static unsigned long long picture(const char *formula, const char *initial)
{
    if (sffe_parse(&cfractalc.userformula, formula) != 0) {
        printf("FAIL   the formula will not parse: %s\n", formula);
        failures++;
        return 0;
    }
    if (sffe_parse(&cfractalc.userinitial, initial) != 0 && initial[0]) {
        printf("FAIL   the initialization will not parse: %s\n", initial);
        failures++;
        return 0;
    }
    sffe_setlocal(&cfractalc);

    iterationfunc fn = userformula()->calculate;
    unsigned long long h = 1469598103934665603ULL;
    for (int iy = 0; iy < 40; iy++)
        for (int ix = 0; ix < 40; ix++) {
            number_t px = (number_t)-1.8 + (number_t)ix * (number_t)0.075;
            number_t py = (number_t)-1.1 + (number_t)iy * (number_t)0.055;
            h ^= fn(px, py, px, py);
            h *= 1099511628211ULL;
        }
    return h;
}

int main(void)
{
    static unsigned int pix[257];
    for (int k = 0; k < 257; k++)
        pix[k] = (unsigned int)(k * 2654435761u);
    cpalette.pixels = pix;
    cpalette.size = 256;
    cpalette.maxentries = 256;
    cpalette.start = 0;
    cpalette.end = 255;

    cfractalc.maxiter = 300;
    cfractalc.bailout = 4;
    cfractalc.incolorspeed = 1.0f;
    cfractalc.outcolorspeed = 1.0f;
    cfractalc.newtonconvergence = 1E-6;
    cfractalc.periodicity_limit = 1e-8;
    cfractalc.range = 2;
    /* the p variables start at the point, which is the setting that makes them
     * worth reading on the first passes */
    cfractalc.pndefault = 1;

    /* The two parsers a context carries, set up as make_fractalc sets them up:
     * they only validate what the user typed, so their variables point at one
     * shared place that nothing reads. */
    static cmplx dummy;
    cfractalc.userformula = sffe_alloc();
    cfractalc.userformula->resolve = sffe_resolve_p;
    sffe_regvar(&cfractalc.userformula, &dummy, "z");
    sffe_regvar(&cfractalc.userformula, &dummy, "c");
    sffe_regvar(&cfractalc.userformula, &dummy, "n");
    sffe_regvar(&cfractalc.userformula, &dummy, "x");
    cfractalc.userinitial = sffe_alloc();
    cfractalc.userinitial->resolve = sffe_resolve_p;
    sffe_regvar(&cfractalc.userinitial, &dummy, "z");
    sffe_regvar(&cfractalc.userinitial, &dummy, "c");
    sffe_regvar(&cfractalc.userinitial, &dummy, "n");
    sffe_regvar(&cfractalc.userinitial, &dummy, "x");

    if (userformula() == NULL) {
        printf("FAIL   there is no user formula in the table\n");
        return 1;
    }

    /* --- the loop draws something, and something that depends on what it was
     * given ------------------------------------------------------------- */
    unsigned long long square = picture("z^2+c", "");
    unsigned long long cube = picture("z^3+c", "");
    check(square != 0 && square != cube,
          "two formulas draw two pictures");

    /* --- what an initialization may name ---------------------------------
     *
     * The dialog accepted an initialization naming z or x and the thread that
     * had to compute it did not, so the initialization was thrown away and the
     * picture came out as though there had been none. Both are the value z
     * starts at, so both say what saying nothing says -- and a picture that
     * comes out the same is only worth checking beside one that comes out
     * different, or an initialization still being thrown away would pass.
     */
    check(picture("z^2+c", "z") == square,
          "an initialization of z is where z starts, and changes nothing");
    check(picture("z^2+c", "x") == square,
          "and so is one of x, which is the same value there");
    check(picture("z^2+c", "z*2") != square,
          "while one that alters it draws a different picture");
    check(picture("z^2+c", "z*2") == picture("z^2+c", "x*2"),
          "z and x alter it alike");

    /* --- and what the formula may name ------------------------------------ */
    check(picture("z^2+c+x*0", "") == square,
          "the formula may name x as well");
    check(picture("z^2+c+n*0", "") == square, "and n");

    /* --- the history of z --------------------------------------------------
     *
     * p1 is the value z had on the pass before. Multiplying it by nothing must
     * leave the picture alone whatever the name is written as, and a formula
     * that actually reads one must draw something that depends on which.
     */
    check(picture("p*0+z^2+c", "") == square,
          "naming p leaves the picture alone when nothing is done with it");
    check(picture("p1*0+z^2+c", "") == square, "and p1 is the same name");
    check(picture("p9999*0+z^2+c", "") == square, "and so is p9999");

    unsigned long long back1 = picture("z^2+c+p1*0.01", "");
    unsigned long long back2 = picture("z^2+c+p2*0.01", "");
    unsigned long long back6 = picture("z^2+c+p6*0.01", "");
    check(back1 != square && back1 != back2 && back2 != back6,
          "how far back a formula reaches changes what it draws");

    /* Further back than the picture has passes, every p stands at what the
     * history started from and stays there, so these two must agree -- and
     * must not agree with one that does come round. */
    unsigned long long deep = picture("z^2+c+p9999*0.01", "");
    check(deep == picture("z^2+c+p400*0.01", "") && deep != back6,
          "beyond the passes the picture has, they all stand where they began");

    /* An initialization may reach into it too, which is the parser the p
     * variables were registered on the wrong side of. */
    check(picture("z^2+c", "p1") == square,
          "an initialization may name p1, which before any pass is where z starts");


    /* --- a figure has to reach the picture --------------------------------
     *
     * These functions hand back a number between zero and one, and a number
     * between zero and one cannot leave a bailout of four. Written alone, the
     * formula therefore sets z to a value that never escapes: every pixel runs
     * to the iteration limit, every pixel gets the same colour, and the picture
     * is the bailout shape in one flat tone with the figure nowhere in it. That
     * is what a bare call draws, and it is worth having said so in a test as
     * well as in the guide, because it looks like the figure is broken.
     *
     * Add it to z and the value drives the escape instead: a point in the
     * gasket climbs to the bailout in a few passes and one in a deep hole takes
     * many, so the levels of the figure come out as the bands of the picture.
     * That is the same thing the noise functions ask for and are documented
     * asking for.
     */
    {
        /* how many different iteration counts a formula draws with */
        struct {
            const char *formula;
            int least;
            const char *what;
        } reach[6] = {
            {"sierpinskyt()", 0, "a bare gasket draws one flat tone"},
            {"sierpinskyc()", 0, "and so does a bare carpet"},
            {"snowflake()", 0, "and a bare snowflake"},
            {"z+sierpinskyt()", 12, "added to z, the gasket draws its levels"},
            {"z+sierpinskyc()", 12, "and the carpet draws its own"},
            {"z+snowflake()", 4, "and the snowflake its fringe"},
        };
        for (int k = 0; k < 6; k++) {
            if (sffe_parse(&cfractalc.userformula, reach[k].formula) != 0) {
                printf("FAIL   will not parse: %s\n", reach[k].formula);
                failures++;
                continue;
            }
            sffe_setlocal(&cfractalc);
            iterationfunc fn = userformula()->calculate;
            /* a coarse grid over the figure, and how many answers it gives */
            unsigned seen[64];
            int n = 0;
            for (int iy = 0; iy < 24 && n < 64; iy++)
                for (int ix = 0; ix < 24 && n < 64; ix++) {
                    number_t px = (number_t)-15 / 10 +
                                  (number_t)ix * (number_t)3 / 23;
                    number_t py = (number_t)-15 / 10 +
                                  (number_t)iy * (number_t)3 / 23;
                    cmplxset(sffe_position, px, py);
                    unsigned v = fn(px, py, px, py);
                    int known = 0;
                    for (int q = 0; q < n; q++)
                        if (seen[q] == v)
                            known = 1;
                    if (!known)
                        seen[n++] = v;
                }
            char line[160];
            sprintf(line, "%s (%d tones)", reach[k].what, n);
            /* a bare call gives one or two, which is a flat picture; added to
             * z it has to give many more than that */
            check(reach[k].least ? n >= reach[k].least : n <= 2, line);
        }
    }

    if (failures)
        printf("\n%d problem(s)\n", failures);
    else
        printf("\nok     the user formula loop draws what it is told to\n");
    return failures != 0;
}
