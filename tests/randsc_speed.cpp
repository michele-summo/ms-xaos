/* What the noise functions cost, and the one thing about it that has to keep
 * holding.
 *
 * A test that asserts nanoseconds asserts which machine it is running on, so
 * this one measures a ratio: the extra time a degradation costs over the same
 * call without one, counted in multiplications. Both halves move together
 * when the machine is fast or slow, or busy, or built at a different
 * precision, and what is left is how much work the degradation is doing --
 * which is the thing that has been wrong twice and is worth pinning.
 *
 * What it should cost is two multiplications, one per component, because the
 * running size is multiplied by the degradation once per pass. Raising the
 * degradation to the power of the pass instead costs a great deal more and
 * grows with the pass number: squaring takes a step per bit of it, the plain
 * loop a step per pass, and a library pow twenty to a hundred multiplications,
 * more at quad. Any of those coming back shows up here.
 *
 * Measured on this machine: carrying reads 1.7 at 64 bits of mantissa and 0.5
 * at 113, squaring 7.2 and 6.1. The numbers are steady to about a tenth over
 * repeated runs.
 */

#include <chrono>
#include <cstdio>

#include "config.h"
#include "number_math.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "misc-f.h"

const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

static int failures = 0;

static void check(int ok, const char *what)
{
    printf("%-6s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        failures++;
}

static sffe *compile(const char *expression)
{
    sffe *parser = sffe_alloc();
    if (!parser || sffe_parse(&parser, expression)) {
        printf("FAIL   cannot parse %s\n", expression);
        failures++;
        return NULL;
    }
    return parser;
}

/* Nanoseconds per evaluation, walking a pixel from pass zero as the engine
 * does. The best of several runs, a minimum being the reading least disturbed
 * by whatever else the machine is doing. */
static double per_eval(sffe *f, unsigned int passes)
{
    const int pixels = 40000 / (int)passes;
    double best = 1e30;
    for (int round = 0; round < 3; round++) {
        volatile number_t sink = 0;
        auto start = std::chrono::steady_clock::now();
        for (int px = 0; px < pixels; px++) {
            GSL_SET_COMPLEX(&sffe_position, (number_t)(px % 97) / 31,
                            (number_t)(px % 89) / 29);
            for (unsigned int n = 0; n < passes; n++) {
                sffe_iteration = n;
                sink = GSL_REAL(sffe_eval(f));
            }
        }
        auto stop = std::chrono::steady_clock::now();
        double ns = std::chrono::duration<double, std::nano>(stop - start).count() /
                    ((double)pixels * passes);
        if (ns < best)
            best = ns;
        (void)sink;
    }
    return best;
}

int main(void)
{
    sffe *plain = compile("{0.3,0.3}");
    sffe *timesone = compile("{0.3,0.3}*{1.1,1.1}");
    sffe *unit = compile("randsch({13,0};{0.3,0.3};{1,1})");
    sffe *fade = compile("randsch({13,0};{0.3,0.3};{0.97,0.97})");
    if (failures)
        return 1;

    /* One measurement thrown away first: the very first reading on a cold
     * cache reads high enough to blunt the comparison below. */
    per_eval(plain, 100);
    per_eval(timesone, 100);
    per_eval(unit, 100);
    per_eval(fade, 100);

    /* The unit the answer is given in: one multiplication. Dividing by a bare
     * evaluation instead would compare arithmetic against parser overhead,
     * and those two do not keep step across precisions -- a multiplication is
     * hardware at 64 bits of mantissa and software at 113, where it costs
     * some twenty times as much while the parser costs the same. */
    double nothing0 = per_eval(plain, 400);
    double onemul = per_eval(timesone, 400) - nothing0;
    printf("       a multiplication costs %.1f ns, an evaluation that does"
           " nothing %.1f\n\n",
           onemul, nothing0);

    printf("       %8s %10s %10s %10s %10s\n", "passes", "nothing", "d=1",
           "d=0.97", "marginal");
    double marginal[3];
    unsigned int passes[3] = {100, 400, 1600};
    for (int i = 0; i < 3; i++) {
        double nothing = per_eval(plain, passes[i]);
        double one = per_eval(unit, passes[i]);
        double faded = per_eval(fade, passes[i]);
        marginal[i] = (faded - one) / onemul;
        printf("       %8u %10.1f %10.1f %10.1f %10.2f\n", passes[i], nothing,
               one, faded, marginal[i]);
    }

    char what[128];

    /* Four sits between what this costs and what the alternatives cost: twice
     * the reading above, and well under the six or seven squaring gives. */
    sprintf(what, "a degradation costs a few multiplications (%.1f)",
            marginal[2]);
    check(marginal[2] < 4.0, what);

    /* And it does not grow: sixteen times the passes, the same cost. Squaring
     * grows by a fifth over that range, a pow by more. */
    sprintf(what, "and does not grow with the passes (%.1f -> %.1f)",
            marginal[0], marginal[2]);
    check(marginal[2] < marginal[0] * 1.5 + 0.3, what);

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures ? 1 : 0;
}
