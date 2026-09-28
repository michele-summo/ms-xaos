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

/* --- going on to a new line ---------------------------------------------- */

/* Part of a name or a number, which a line may not end inside. */
static bool in_name(char16_t ch)
{
    return (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
           (ch >= '0' && ch <= '9') || ch == '_' || ch == '.';
}

/* The second half of a surrogate pair, which cannot begin a line. */
static bool low_surrogate(char16_t ch) { return ch >= 0xDC00 && ch <= 0xDFFF; }

std::vector<int> formula_wrap(const char16_t *text, int length, double avail,
                              const std::function<double(int, int)> &width,
                              int most)
{
    std::vector<int> ends;
    if (length <= 0) {
        ends.push_back(0);
        return ends;
    }
    int start = 0;
    while (start < length && (int)ends.size() < most) {
        /* The furthest end that fits at all, and the furthest that fits and
         * falls between names. A line takes one character however wide, and
         * then goes on while the next still fits. */
        int fits = start, between = start;
        for (int end = start + 1; end <= length; end++) {
            if (end < length && low_surrogate(text[end]))
                continue;
            if (fits > start && width(start, end) > avail)
                break;
            fits = end;
            if (end == length || !in_name(text[end - 1]) || !in_name(text[end]))
                between = end;
        }
        const int end = between > start ? between : fits;
        ends.push_back(end);
        start = end;
    }
    return ends;
}

/* --- the history window -------------------------------------------------- */

std::u16string formula_history_label(
    const std::u16string &formula, const std::u16string &suffix, double avail,
    const std::function<double(const std::u16string &)> &width)
{
    if (width(formula) <= avail)
        return formula;
    /* The longest start that fits with the suffix after it. A longer start is
     * never narrower, so the answer is found by halving. */
    size_t lo = 0, hi = formula.size();
    while (lo < hi) {
        const size_t mid = (lo + hi + 1) / 2;
        if (width(formula.substr(0, mid) + suffix) <= avail)
            lo = mid;
        else
            hi = mid - 1;
    }
    /* Cut between names, as the bar goes on to a new line: a number cut after
     * its point would run into the dots ("0...."), and a cut name reads as
     * another. Only a name wider than the whole room is cut, and never
     * between the halves of a surrogate pair. */
    size_t cut = lo;
    while (cut > 0 && cut < formula.size() &&
           ((in_name(formula[cut - 1]) && in_name(formula[cut])) ||
            low_surrogate(formula[cut])))
        cut--;
    if (cut == 0) {
        cut = lo;
        if (cut > 0 && cut < formula.size() && low_surrogate(formula[cut]))
            cut--;
    }
    while (cut > 0 && formula[cut - 1] == u' ')
        cut--;
    return formula.substr(0, cut) + suffix;
}

/* --- the status ----------------------------------------------------------- */

std::vector<std::u16string> formula_status_lines(
    const std::u16string &formula, const std::u16string &suffix, double avail,
    const std::function<double(const std::u16string &)> &width)
{
    /* Two lines at most are laid out, so a formula of thousands of characters
     * costs two lines' measuring rather than all of its own. */
    const std::vector<int> ends = formula_wrap(
        formula.data(), (int)formula.size(), avail,
        [&](int start, int end) {
            return width(formula.substr(start, end - start));
        },
        2);
    std::vector<std::u16string> lines;
    lines.push_back(formula.substr(0, ends[0]));
    if (ends.size() > 1) {
        const std::u16string rest = formula.substr(ends[0]);
        lines.push_back(ends[1] == (int)formula.size()
                            ? rest
                            : formula_history_label(rest, suffix, avail, width));
    }
    return lines;
}
