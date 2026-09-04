/* What reaching back into the history of z costs, and the one thing about it
 * that has to keep holding: that it does not depend on how far back you reach.
 *
 * A test that asserts nanoseconds asserts which machine it is running on, so
 * this one measures a ratio. p1 and p9999 are the same amount of work -- one
 * value written into the ring, one read out of it -- and differ only in the
 * number the reading is offset by, so a pass that names p9999 must cost what a
 * pass that names p1 costs. Shifting an array along instead costs a copy per
 * place per pass: it is what this replaced, it is why the history used to stop
 * at six, and at nine thousand places it would show up here as a factor of
 * hundreds rather than a few percent.
 *
 * Measured on this machine: p9999 against p1 reads 0.98 at 64 bits of mantissa
 * and 0.99 at 113, and a formula naming none of them 0.90 and 0.96. Steady to
 * a hundredth over repeated runs.
 */

#include <chrono>
#include <cstdio>
#include <cstring>

#include "config.h"
#include "number_math.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "phist.h"

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

static cmplx var_z, var_c;

static sffe *compile(const char *expression)
{
    sffe *parser = sffe_alloc();
    parser->resolve = sffe_resolve_p;
    sffe_regvar(&parser, &var_z, "z");
    sffe_regvar(&parser, &var_c, "c");
    if (sffe_parse(&parser, expression)) {
        printf("FAIL   cannot parse %s\n", expression);
        failures++;
        return NULL;
    }
    return parser;
}

/* Nanoseconds per pass, walking a pixel from the first pass as the engine
 * does: reset the history, then evaluate and push, over and over. The best of
 * several runs, a minimum being the reading least disturbed by whatever else
 * the machine is doing. */
static double per_pass(sffe *f, unsigned int passes)
{
    const int pixels = 200000 / (int)passes;
    double best = 1e30;
    for (int round = 0; round < 3; round++) {
        volatile number_t sink = 0;
        auto start = std::chrono::steady_clock::now();
        for (int px = 0; px < pixels; px++) {
            GSL_SET_COMPLEX(&var_c, (number_t)(px % 97) / 311,
                            (number_t)(px % 89) / 293);
            GSL_SET_COMPLEX(&var_z, 0, 0);
            sffe_phist_reset(GSL_REAL(var_c), GSL_IMAG(var_c));
            for (unsigned int n = 0; n < passes; n++) {
                sffe_iteration = n;
                cmplx next = sffe_eval(f);
                sffe_phist_push(GSL_REAL(var_z), GSL_IMAG(var_z));
                var_z = next;
                sink = GSL_REAL(var_z);
            }
        }
        auto stop = std::chrono::steady_clock::now();
        double ns =
            std::chrono::duration<double, std::nano>(stop - start).count() /
            ((double)pixels * passes);
        if (ns < best)
            best = ns;
        (void)sink;
    }
    return best;
}

int main(void)
{
    /* The same arithmetic in all four, so that what is left between them is
     * the history and nothing else. "n" would be another variable to read;
     * these differ only in which place is read. */
    sffe *none = compile("z*0.5+c+c*0");
    sffe *shallow = compile("z*0.5+c+p1*0");
    sffe *middle = compile("z*0.5+c+p64*0");
    sffe *deep = compile("z*0.5+c+p9999*0");
    if (failures)
        return 1;

    struct {
        const char *what;
        sffe *f;
        double ns;
    } row[4] = {{"names no p", none, 0},
                {"names p1", shallow, 0},
                {"names p64", middle, 0},
                {"names p9999", deep, 0}};

    /* One measurement thrown away first: the very first reading on a cold
     * cache reads high enough to blunt the comparison below. */
    for (int i = 0; i < 4; i++) {
        sffe_ptaps_build(row[i].f, NULL);
        per_pass(row[i].f, 100);
    }

    printf("       %-14s %8s %10s %10s\n", "formula", "depth", "ns/pass",
           "places");
    for (int i = 0; i < 4; i++) {
        sffe_ptaps_build(row[i].f, NULL);
        row[i].ns = per_pass(row[i].f, 800);
        printf("       %-14s %8u %10.1f %10u\n", row[i].what, sffe_ph.depth,
               row[i].ns, sffe_ph.tapcount);
    }
    printf("\n");

    char what[160];

    /* The claim: the depth costs nothing. Ten thousand places against one, and
     * the shift this replaced would have been reading the whole ring. */
    double ratio = row[3].ns / row[1].ns;
    sprintf(what,
            "reaching back nine thousand passes costs what reaching back one "
            "costs (%.2f)",
            ratio);
    check(ratio < 1.25, what);

    /* And the step from one to sixty-four, which is where an array shifted
     * along by one would already have been sixty-four times the work. */
    double mid = row[2].ns / row[1].ns;
    sprintf(what, "and so does reaching back sixty-four (%.2f)", mid);
    check(mid < 1.25, what);

    /* A formula that names none of them must not pay for the machinery: with
     * no places to fill there is nothing to write and nothing to read, which
     * is what the shift could not say. */
    double idle = row[0].ns / row[1].ns;
    sprintf(what, "and a formula that names none pays less than one that does (%.2f)",
            idle);
    check(idle < 1.02, what);

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    else
        printf("\nok     how far back costs nothing\n");
    return failures != 0;
}
