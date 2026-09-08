/* What the cursor is standing in, for the line under the formula bar.
 *
 * The bar takes one line of text and the reference lives behind a menu, so a
 * user writing randsc(7,{1,1},,,,0.4) has to remember what the fourth place
 * means or go and look. This reads the text around the cursor and says which
 * call it is inside and which argument of it, and the bar shows that call's
 * arguments with the one the cursor is in emphasised.
 *
 * It is deliberately not a parser. A formula is half written most of the time
 * it is looked at -- brackets unclosed, arguments missing, a name that is not
 * a function yet -- and a hint that only appears once the formula is correct
 * would be silent exactly when it is wanted. So this walks backwards over the
 * text looking for an open bracket nothing has closed, and takes the name in
 * front of it. Anything it cannot make sense of is nothing to show, never an
 * error.
 *
 * Brackets and braces share one depth count, which a real parser could not do:
 * "f({a)" would be read as balanced. That is the right trade here -- getting
 * it wrong shows the wrong signature for a moment on text that means nothing
 * anyway, where telling the two apart would cost the simplicity that makes
 * this readable.
 */

#include <cctype>
#include <cstddef>
#include <cstring>

#include "formulahelp.h"

/* A name is what the parser will read as one: a letter or an underscore, then
 * letters, digits and underscores. */
static int hint_isname(char c)
{
    return isalnum((unsigned char)c) || c == '_';
}

static int hint_isnamestart(char c)
{
    return isalpha((unsigned char)c) || c == '_';
}

/* The name ending just before position i, spaces between it and i skipped.
 * Returns its length and writes where it starts; zero when there is none. */
static int hint_name_before(const char *text, int i, int *start)
{
    while (i > 0 && isspace((unsigned char)text[i - 1]))
        i--;
    int end = i;
    while (i > 0 && hint_isname(text[i - 1]))
        i--;
    if (i == end || !hint_isnamestart(text[i]))
        return 0;
    *start = i;
    return end - i;
}

/* The row describing that name, or NULL. The parser lowers the case of a
 * formula before it reads it, so LOGN and logn are one function and the table
 * holds the lower one. */
static const struct formula_help_row *hint_lookup(const char *text, int start,
                                                  int len)
{
    for (const struct formula_help_row *r = formula_help_functions;
         r->name || r->section; r++) {
        if (!r->name || (int)strlen(r->name) != len)
            continue;
        int k = 0;
        while (k < len && tolower((unsigned char)text[start + k]) ==
                              tolower((unsigned char)r->name[k]))
            k++;
        if (k == len)
            return r;
    }
    return NULL;
}

struct formula_hint formula_hint_at(const char *text, int cursor)
{
    struct formula_hint hint;
    hint.fn = NULL;
    hint.arg = -1;
    if (!text)
        return hint;

    int len = (int)strlen(text);
    if (cursor < 0)
        cursor = 0;
    if (cursor > len)
        cursor = len;

    /* Backwards for a bracket nothing has closed. A closing bracket met on the
     * way stands for a call that is finished and has nothing to do with where
     * the cursor is, so it is skipped whole. An opening brace met at depth
     * nought is a complex number the cursor is inside -- transparent here,
     * since what is wanted is the call around it. */
    int open = -1;
    int depth = 0;
    for (int i = cursor - 1; i >= 0; i--) {
        char c = text[i];
        if (c == ')' || c == '}') {
            depth++;
        } else if (c == '(' || c == '{') {
            if (depth > 0)
                depth--;
            else if (c == '(') {
                open = i;
                break;
            }
        }
    }

    if (open < 0) {
        /* No call open. There is still one thing worth showing: a cursor that
         * has just finished a name whose bracket is already there -- "randsc_("
         * -- which is a user about to write the first argument. The call is
         * named, no argument is current. */
        int j = cursor;
        while (j < len && isspace((unsigned char)text[j]))
            j++;
        if (j < len && text[j] == '(') {
            int start;
            int n = hint_name_before(text, cursor, &start);
            if (n)
                hint.fn = hint_lookup(text, start, n);
        }
        return hint;
    }

    int start;
    int n = hint_name_before(text, open, &start);
    if (!n)
        return hint; /* a bracket that groups rather than calls */
    hint.fn = hint_lookup(text, start, n);
    if (!hint.fn)
        return hint;

    /* Which argument: the separators between that bracket and the cursor that
     * belong to this call rather than to something written inside it. A comma
     * inside a nested call, or between the two parts of a complex number,
     * belongs to that and is not counted -- {1,1} as one argument is the shape
     * half the arguments in these formulas have. */
    int arg = 0;
    depth = 0;
    for (int i = open + 1; i < cursor; i++) {
        char c = text[i];
        if (c == '(' || c == '{')
            depth++;
        else if (c == ')' || c == '}') {
            if (depth > 0)
                depth--;
        } else if (c == ',' && depth == 0)
            arg++;
    }
    hint.arg = arg;
    return hint;
}

int formula_hint_arg_span(const char *args, int arg, int *start, int *len)
{
    if (!args || arg < 0)
        return 0;

    /* The arguments are written separated by commas and none of them holds one
     * of its own, so this is a split and nothing more. */
    int i = 0, n = (int)strlen(args), k = 0;
    int from = 0, last_from = -1, last_len = 0;
    for (;; i++) {
        if (i == n || args[i] == ',') {
            int s = from, e = i;
            while (s < e && isspace((unsigned char)args[s]))
                s++;
            while (e > s && isspace((unsigned char)args[e - 1]))
                e--;
            if (e > s) {
                if (k == arg) {
                    *start = s;
                    *len = e - s;
                    return 1;
                }
                last_from = s;
                last_len = e - s;
                k++;
            }
            if (i == n)
                break;
            from = i + 1;
        }
    }

    /* Past the last one. A call that ends in "..." takes as many more as one
     * likes -- poly is written "z, k1, k2, ..." -- so the tail stays current
     * however many are written, which is what it means. Anything else has run
     * out of arguments and nothing is emphasised. */
    if (last_from >= 0 && last_len == 3 && !strncmp(args + last_from, "...", 3)) {
        *start = last_from;
        *len = last_len;
        return 1;
    }
    return 0;
}
