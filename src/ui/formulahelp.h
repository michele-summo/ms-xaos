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

/* --- the line under the formula bar --------------------------------------
 *
 * What the cursor is standing in: the call it is inside, and which of that
 * call's arguments. Written without Qt, like the tables above, so that the
 * walk over the text can be tested next to the parser rather than by opening
 * a dialog and looking. See formulahint.cpp.
 */
struct formula_hint {
    /* The call the cursor is in, or NULL when there is nothing to show. */
    const struct formula_help_row *fn;
    /* Which argument the cursor stands in, counting from nought, or -1 when
     * the call is named but no argument is current. */
    int arg;
};

struct formula_hint formula_hint_at(const char *text, int cursor);

/* Where one argument stands inside a "takes" string, so that it can be shown
 * apart from the rest. Returns 0 when that argument is not there. */
int formula_hint_arg_span(const char *args, int arg, int *start, int *len);

/* The two of them worded into the line the bar shows, with the argument the
 * cursor is in emphasised, and empty when there is nothing to show.
 *
 * QString is named and not included: this header is linked by a test that has
 * no Qt, and a declaration nothing there calls costs that test nothing. The
 * definition is in formulahelp.cpp, which has Qt and is not in the test.
 */
class QString;
QString formula_hint_text(const QString &text, int cursor);

struct uih_context;
void ui_formulahelp(struct uih_context *uih);

#endif // FORMULAHELP_H
