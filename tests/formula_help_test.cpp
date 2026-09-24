/* The help window lists the functions a user formula may call. That list is a
 * second copy of the parser's own table, and a second copy drifts: a function
 * added to sfcmplxfunc and not described here would simply never appear in the
 * reference, silently, and one renamed would be listed under a name that no
 * longer parses.
 *
 * So compare the two. Every name the parser accepts must be described, and
 * every name described must be one the parser accepts.
 *
 * The Tilings tab is a third copy, of randsctile this time, and is compared
 * the same way: a caption for every tiling the parser draws, and a picture
 * of each that is still what the parser draws.
 */

#include <cstdio>
#include <cstring>

#include "config.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "formulahelp.h"
#include "randsctile-thumbnails.h"

const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

static int failures = 0;

/* How many arguments a "takes" line describes: the parts between commas,
 * with a trailing "..." meaning any number more. */
static int described_arity(const char *args, int *open_ended)
{
    *open_ended = 0;
    if (args == NULL || *args == 0)
        return 0;
    int parts = 1;
    for (const char *c = args; *c; c++)
        if (*c == ',')
            parts++;
    if (strstr(args, "...") != NULL) {
        *open_ended = 1;
        parts--; /* the "..." is not an argument of its own */
    }
    return parts;
}

static number_t tiling_at(sffe *p, double x, double y)
{
    GSL_SET_COMPLEX(&sffe_position, (number_t)x, (number_t)y);
    sffe_iteration = 0;
    return GSL_REAL(sffe_eval(p));
}

#define TILINGS_DIR XAOS_SOURCE_DIR "/src/ui/images/tilings/"
#define REDRAW "run cmake --build <build directory> --target randsctile-thumbnails"

/* The Tilings tab: its captions against the parser, its pictures against
 * what the parser draws now. */
static void check_tilings(void)
{
    /* The tilings the parser has are the numbers before the first that
     * draws nothing, which is how the program that draws the pictures counts
     * them too. */
    enum { MOST = 256 };
    uint64_t print[MOST + 1];
    int drawn = 0;
    for (int k = 1; k <= MOST; k++) {
        char call[64];
        snprintf(call, sizeof call, RANDSCTILE_THUMB_CALL, k);
        sffe *p = sffe_alloc();
        if (!p || sffe_parse(&p, call)) {
            printf("FAIL   cannot parse %s\n", call);
            failures++;
            return;
        }
        if (tiling_at(p, 0.3, 0.7) == 0) {
            sffe_free(&p);
            break;
        }
        print[k] = randsctile_thumb_fingerprint(
            [p](double x, double y) { return (double)tiling_at(p, x, y); });
        sffe_free(&p);
        drawn = k;
    }

    const struct formula_help_row *rows;
    int listed = formula_help_tilings(&rows);
    if (listed != drawn) {
        printf("FAIL   the Values tab lists %d tilings under randsctile, and "
               "the parser draws %d\n",
               listed, drawn);
        failures++;
    }
    for (int k = 1; k <= listed; k++)
        if (rows[k - 1].summary == NULL || rows[k - 1].summary[0] == '\0') {
            printf("FAIL   tiling %d is listed with no description\n", k);
            failures++;
        }

    /* The fingerprints the pictures were drawn with, one line a tiling. */
    FILE *f = fopen(TILINGS_DIR "fingerprints.txt", "r");
    if (!f) {
        printf("FAIL   no " TILINGS_DIR "fingerprints.txt: " REDRAW "\n");
        failures++;
        return;
    }
    int stored = 0, stale = 0;
    char line[256];
    while (fgets(line, sizeof line, f)) {
        int k;
        unsigned long long was;
        if (line[0] == '#' || sscanf(line, "%d %llx", &k, &was) != 2)
            continue;
        if (k != stored + 1) {
            printf("FAIL   fingerprints.txt has tiling %d after %d\n", k,
                   stored);
            failures++;
            break;
        }
        stored = k;
        if (k <= drawn && was != print[k]) {
            printf("FAIL   the thumbnail of tiling %d is not what randsctile "
                   "draws now\n",
                   k);
            stale++;
        }
    }
    fclose(f);
    if (stored != drawn) {
        printf("FAIL   there are thumbnails of %d tilings, and the parser draws "
               "%d\n",
               stored, drawn);
        failures++;
    }
    if (stale || stored != drawn) {
        printf("       " REDRAW "\n");
        failures += stale;
    }

    /* And each picture is there, and is in the resource file that puts it in
     * the binary. */
    f = fopen(TILINGS_DIR "tilings.qrc", "rb");
    static char qrc[1 << 16];
    size_t size = f ? fread(qrc, 1, sizeof qrc - 1, f) : 0;
    if (f)
        fclose(f);
    qrc[size] = 0;
    for (int k = 1; k <= drawn; k++) {
        char name[64], path[512];
        snprintf(name, sizeof name, "<file>tiling-%02d.png</file>", k);
        snprintf(path, sizeof path, TILINGS_DIR "tiling-%02d.png", k);
        FILE *png = fopen(path, "rb");
        if (!png || !strstr(qrc, name)) {
            printf("FAIL   the thumbnail of tiling %d is %s: " REDRAW "\n", k,
                   png ? "not in tilings.qrc" : "missing");
            failures++;
        }
        if (png)
            fclose(png);
    }
    if (!failures)
        printf("ok     %d tilings, each captioned and drawn as randsctile "
               "draws it\n",
               drawn);
}

int main(void)
{
    /* The first sffnctsfirst entries are the operators, which sffe_function
     * never reaches by name -- but they are worth describing, and the help
     * table lists them, so the whole table is compared. */
    for (int i = 0; i < sffnctscount; i++) {
        const char *name = sfcmplxfunc[i].name;
        int found = 0;
        for (const struct formula_help_row *r = formula_help_functions;
             r->name || r->section; r++)
            if (r->name && !strcmp(r->name, name)) {
                found = 1;
                break;
            }
        if (!found) {
            printf("FAIL   \"%s\" is in the parser and not in the help table\n",
                   name);
            failures++;
        }
    }

    /* And each must be described as taking what it takes. The arity is the
     * parser's own, so a function that gains or loses an argument cannot go on
     * being described with the one it used to have. */
    for (int i = 0; i < sffnctscount; i++) {
        const char *name = sfcmplxfunc[i].name;
        for (const struct formula_help_row *r = formula_help_functions;
             r->name || r->section; r++) {
            if (!r->name || strcmp(r->name, name))
                continue;
            int open_ended = 0;
            int described = described_arity(r->args, &open_ended);
            int actual = sfcmplxfunc[i].parcnt;
            if (actual == SFFE_VARIADIC) {
                if (!open_ended && described == 0) {
                    printf("FAIL   \"%s\" takes any number of arguments and "
                           "says nothing about them\n",
                           name);
                    failures++;
                }
            } else if (described != actual && described != 0) {
                /* zero stands for the operators, which are shown whole as
                 * "a + b" and have nothing to list */
                printf("FAIL   \"%s\" takes %d, described as taking %d\n", name,
                       actual, described);
                failures++;
            }
            break;
        }
    }

    for (const struct formula_help_row *r = formula_help_functions;
         r->name || r->section; r++) {
        if (!r->name)
            continue;
        int found = 0;
        for (int i = 0; i < sffnctscount; i++)
            if (!strcmp(sfcmplxfunc[i].name, r->name)) {
                found = 1;
                break;
            }
        if (!found) {
            printf("FAIL   \"%s\" is in the help table and not in the parser\n",
                   r->name);
            failures++;
        }
        if (r->summary == NULL || r->summary[0] == '\0') {
            printf("FAIL   \"%s\" is listed with no description\n", r->name);
            failures++;
        }
    }

    /* A section with no rows under it would be a heading over nothing. */
    for (const struct formula_help_row *r = formula_help_functions;
         r->name || r->section; r++)
        if (r->section && !(r + 1)->name) {
            printf("FAIL   section \"%s\" has nothing under it\n", r->section);
            failures++;
        }

    check_tilings();

    if (failures)
        printf("\n%d problem(s)\n", failures);
    else
        printf("ok     %d functions, each described exactly once\n",
               sffnctscount);
    return failures != 0;
}
