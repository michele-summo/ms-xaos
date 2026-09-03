#ifndef FORMULAHELP_H
#define FORMULAHELP_H

#include <cstddef> /* NULL */

/* The reference shown by Help -> User formula reference.
 *
 * The tables are plain data and mention no Qt, so that a test can link them
 * next to the parser and check that every function the parser accepts is
 * described here and nothing is described that the parser does not have.
 */

struct formula_help_row {
    const char *name;    /* NULL on a section heading */
    /* What the call takes, in order: "a; b" for two, "seed; size" where the
     * position means something, empty where there is nothing to take. A test
     * counts these against the arity the parser holds, so a function that
     * gains or loses an argument cannot go on being described with the old
     * one. NULL on a section heading. */
    const char *args;
    const char *summary; /* NULL on a section heading */
    const char *section; /* set only on a section heading */
};

/* Terminated by a row with all three NULL. */
extern const struct formula_help_row formula_help_functions[];
extern const struct formula_help_row formula_help_variables[];
extern const struct formula_help_row formula_help_notation[];
/* Not compared against the parser: these are argument values, not names. */
extern const struct formula_help_row formula_help_values[];
/* Not compared against the parser: these are argument values, not names. */
extern const struct formula_help_row formula_help_values[];

struct uih_context;
void ui_formulahelp(struct uih_context *uih);

#endif // FORMULAHELP_H
