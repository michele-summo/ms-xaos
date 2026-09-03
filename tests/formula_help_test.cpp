/* The help window lists the functions a user formula may call. That list is a
 * second copy of the parser's own table, and a second copy drifts: a function
 * added to sfcmplxfunc and not described here would simply never appear in the
 * reference, silently, and one renamed would be listed under a name that no
 * longer parses.
 *
 * So compare the two. Every name the parser accepts must be described, and
 * every name described must be one the parser accepts.
 */

#include <cstdio>
#include <cstring>

#include "config.h"
#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include "formulahelp.h"

const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

static int failures = 0;

/* How many arguments a "takes" line describes: the parts between semicolons,
 * with a trailing "..." meaning any number more. */
static int described_arity(const char *args, int *open_ended)
{
    *open_ended = 0;
    if (args == NULL || *args == 0)
        return 0;
    int parts = 1;
    for (const char *c = args; *c; c++)
        if (*c == ';')
            parts++;
    if (strstr(args, "...") != NULL) {
        *open_ended = 1;
        parts--; /* the "..." is not an argument of its own */
    }
    return parts;
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

    if (failures)
        printf("\n%d problem(s)\n", failures);
    else
        printf("ok     %d functions, each described exactly once\n",
               sffnctscount);
    return failures != 0;
}
