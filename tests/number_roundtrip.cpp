/* A saved position is only as good as the round trip through text: XaoS writes
 * the view coordinates with save_float and reads them back with xstrtonum, and
 * if either side carries fewer digits than number_t holds, reloading a file
 * lands somewhere slightly different from where it was saved. Deep zoom is
 * exactly where that shows.
 *
 * This checks the two halves agree at whatever precision the build was
 * compiled for, using the same NUMBER_DIGITS the real formatting code does.
 */
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "config.h"
#include "misc-f.h"

#ifdef USE_FLOAT128
#include <quadmath.h>
#endif

/* Same shape as save_float in src/ui-hlp/save.cpp. */
static void format_number(char *out, size_t n, number_t v)
{
    char fs[10];
#ifdef USE_FLOAT128
    snprintf(fs, sizeof(fs), "%%.%iQG", NUMBER_DIGITS);
    quadmath_snprintf(out, n, fs, (__float128)v);
#else
#ifdef USE_LONG_DOUBLE
    snprintf(fs, sizeof(fs), "%%.%iLG", NUMBER_DIGITS);
    snprintf(out, n, fs, (long double)v);
#else
    snprintf(fs, sizeof(fs), "%%.%iG", NUMBER_DIGITS);
    snprintf(out, n, fs, (double)v);
#endif
#endif
}

static int failures = 0;

static void roundtrip(const char *what, number_t v)
{
    char text[256];
    char *end;
    format_number(text, sizeof(text), v);
    number_t back = xstrtonum(text, &end);
    if (back != v) {
        printf("FAIL  %s did not survive the round trip, wrote \"%s\"\n", what,
               text);
        failures++;
    }
}

int main(void)
{
    printf("number_t carries %d significant digits\n", NUMBER_DIGITS);

    roundtrip("zero", (number_t)0);
    roundtrip("one", (number_t)1);
    roundtrip("minus one", (number_t)-1);

    /* A view centre from a deep zoom: enough digits that a double would round
     * it away. Built by division so the compiler cannot fold it to a literal
     * of some other type. */
    number_t deep = (number_t)1 / (number_t)3;
    roundtrip("1/3", deep);
    roundtrip("1/3 shifted small", deep * (number_t)1e-20);
    roundtrip("1/3 shifted large", deep * (number_t)1e20);

    /* The last representable step: if the format drops a digit, adding one ulp
     * and printing gives back the value we started from. */
    number_t one = 1;
    number_t eps = 1;
    while ((number_t)(one + eps / 2) != one)
        eps /= 2;
    roundtrip("1 + one ulp", one + eps);
    if (xstrtonum("1", NULL) + eps == (number_t)1) {
        printf("FAIL  epsilon search produced no representable step\n");
        failures++;
    }

    /* And the round trip has to distinguish neighbours, not merely return
     * something close. */
    {
        char a[256], b[256];
        format_number(a, sizeof(a), one);
        format_number(b, sizeof(b), one + eps);
        if (!strcmp(a, b)) {
            printf("FAIL  1 and 1+ulp print identically as \"%s\": the format "
                   "carries fewer digits than number_t\n",
                   a);
            failures++;
        } else {
            printf("1 and 1+ulp print differently, as they must\n");
        }
    }

    printf("%s\n", failures ? "FAILED" : "ok");
    return failures ? 1 : 0;
}
