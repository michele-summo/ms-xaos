#ifndef FORMULAHELP_H
#define FORMULAHELP_H

#include <climits>
#include <cstddef> /* NULL */
#include <functional>
#include <string>
#include <vector>

/* The reference shown by Help -> User formula reference.
 *
 * The tables are plain data and mention no Qt, so that a test can link them
 * next to the parser and check that every function the parser accepts is
 * described here and nothing is described that the parser does not have.
 */

/* The descriptions and section headings are marked for lupdate in the
 * "FormulaHelp" context, and the dialog translates them when it shows them;
 * names and argument lists are syntax and stay as written. Qt defines the
 * marker as the text itself, and so does this, for the test that links the
 * tables with no Qt. */
#ifndef QT_TRANSLATE_NOOP
#define QT_TRANSLATE_NOOP(scope, x) x
#endif

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

/* The tilings randsctile draws, as the Values tab lists them: the numbered
 * rows under its heading that begins "randsctile:", tiling k at rows[k - 1].
 * Returns how many, and points *rows at the first. The Tilings tab captions
 * its pictures with these, and a test checks they are every tiling the
 * parser has -- no more, no fewer. */
int formula_help_tilings(const struct formula_help_row **rows);

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

/* --- the formula bar: going on to a new line ---------------------------------
 *
 * A formula too long for the bar goes on to the next line between names, never
 * inside one: a name or a number is a run of letters, digits, '_' and '.', and
 * a line may end anywhere else. Only a name wider than the whole bar is cut,
 * and then as late as it fits. Every line holds at least one character, so an
 * empty formula is one empty line and the bar never stops short.
 *
 * The text is UTF-16, as Qt holds it, and is never cut between the two halves
 * of a surrogate pair. width(start, end) is how wide text[start, end) is drawn,
 * and avail how wide a line may be. Returns where each line ends, the last one
 * at length -- or, given most, where the first most lines end and no more, the
 * last of them short of length if the text goes on: measuring the rest of a
 * long formula for lines nobody will draw is what most saves. Written without
 * Qt, like the walk above, so that the rule can be tested on its own;
 * formulaedit.cpp is what lays the bar out by it.
 */
std::vector<int> formula_wrap(const char16_t *text, int length, double avail,
                              const std::function<double(int, int)> &width,
                              int most = INT_MAX);

/* A formula as the history window lists it: whole if it fits in avail, and
 * otherwise as much of its start as fits followed by suffix -- which the
 * window makes "..." and how many characters the formula has, so that two
 * long formulas that begin alike can still be told apart. width(text) is how
 * wide text is drawn. The start ends between names, as a line of the bar
 * does, without the spaces before the cut; only a name wider than the whole
 * room is cut, never between the halves of a surrogate pair. It is empty when
 * not even the suffix fits.
 */
std::u16string formula_history_label(
    const std::u16string &formula, const std::u16string &suffix, double avail,
    const std::function<double(const std::u16string &)> &width);

/* A formula as View -> Status shows it, avail wide: one line if it fits, and
 * otherwise two, going on between names as the bar does. If two do not hold it
 * either, the second is as much of the rest as fits followed by suffix -- which
 * the status makes "..." and how many characters the formula has, as the
 * history window does. Returns the one or two lines. width(text) is how wide
 * text is drawn.
 */
std::vector<std::u16string> formula_status_lines(
    const std::u16string &formula, const std::u16string &suffix, double avail,
    const std::function<double(const std::u16string &)> &width);

/* --- the Utilities tab: a random poly ----------------------------------------
 *
 * poly(z, k1, ..., km) with m drawn between a shortest and a longest length,
 * and every coefficient a complex number whose magnitude is drawn between a
 * smallest and a largest one and whose argument is drawn at random. Written
 * without Qt, like the walk above, so that a test can hand every formula it
 * makes to the parser. See formulapoly.cpp.
 */

/* poly takes at most 255 arguments -- the parser counts them in a byte -- and
 * z is one of them. */
#define FORMULA_POLY_LONGEST 254
/* The largest magnitude a coefficient may be given: written out in figures,
 * which is what the parser reads, anything larger is a long row of digits
 * saying nothing a fractal can use. */
#define FORMULA_POLY_LARGEST 1000000

enum formula_poly_reading {
    FORMULA_POLY_OK = 0,
    FORMULA_POLY_EMPTY,       /* nothing written, or only spaces */
    FORMULA_POLY_NOT_NUMBER,  /* something that does not read as a number */
    FORMULA_POLY_NOT_WHOLE,   /* a length, but with a fraction */
    FORMULA_POLY_TOO_SHORT,   /* a length below 1 */
    FORMULA_POLY_TOO_LONG,    /* a length above FORMULA_POLY_LONGEST */
    FORMULA_POLY_NEGATIVE,    /* a magnitude below nought */
    FORMULA_POLY_TOO_LARGE    /* a magnitude above FORMULA_POLY_LARGEST */
};

/* One of the two lengths as written in its field: FORMULA_POLY_OK and the
 * length in *length, or why there is none. */
int formula_poly_length(const char *text, int *length);

/* One of the two magnitudes as written in its field: FORMULA_POLY_OK and the
 * magnitude in *magnitude -- the fallback when the field is empty -- or why
 * there is none. */
int formula_poly_magnitude(const char *text, double fallback,
                           double *magnitude);

/* A call of a length from shortest to longest, its coefficients of a
 * magnitude from smallest to largest, all four already read; the same seed
 * gives the same call. The length drawn is put in *length. */
std::string formula_poly_random(int shortest, int longest, double smallest,
                                double largest, unsigned long long seed,
                                int *length);

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
