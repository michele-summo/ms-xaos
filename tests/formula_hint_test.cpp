/* The line under the formula bar: what the cursor is standing in.
 *
 * A walk over half-written text, which is the only text it will ever see, so
 * what has to hold is that it never lies and never insists: the right call and
 * the right argument where there is one, and nothing at all where there is
 * not. The cursor is written as an underscore in the cases below and taken out
 * before the text is handed over, so that each case reads the way it looks in
 * the bar.
 */

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "formulahelp.h"

static int failures = 0;

static void check(int condition, const char *what)
{
    if (condition) {
        printf("ok     %s\n", what);
    } else {
        printf("FAIL   %s\n", what);
        failures++;
    }
}

/* the text with the underscore taken out, and where it stood */
static int uncursor(const char *marked, char *out, size_t size)
{
    size_t j = 0;
    int cursor = -1;
    for (size_t i = 0; marked[i] && j + 1 < size; i++) {
        if (marked[i] == '_' && cursor < 0) {
            cursor = (int)j;
            continue;
        }
        out[j++] = marked[i];
    }
    out[j] = 0;
    return cursor < 0 ? (int)j : cursor;
}

/* One case: the text with the cursor marked, the call it should name (NULL for
 * none) and the argument it should be standing in (-1 for none). */
static void one(const char *marked, const char *fn, int arg)
{
    char text[256], what[512];
    int cursor = uncursor(marked, text, sizeof text);
    struct formula_hint h = formula_hint_at(text, cursor);

    const char *got = h.fn ? h.fn->name : "(nothing)";
    sprintf(what, "\"%s\" -> %s, argument %d", marked, got, h.arg);
    int right = fn ? (h.fn && !strcmp(h.fn->name, fn) && h.arg == arg)
                   : (h.fn == NULL);
    check(right, what);
}

int main(void)
{
    /* --- the four the bar was asked for ----------------------------------*/
    one("randsc(a,b) + _c", NULL, 0);      /* outside every call: nothing */
    one("randsc_(a,b) + c", "randsc", -1); /* named, no argument yet */
    one("randsc(_a,b) + c", "randsc", 0);
    one("randsc(a_,b) + c", "randsc", 0);
    one("randsc(a,_b) + c", "randsc", 1);
    one("randsc(a,b_) + c", "randsc", 1);

    /* --- half-written text, which is what it will mostly see -------------*/
    one("randsc(_", "randsc", 0);
    one("randsc(_a,b", "randsc", 0);
    one("randsc(a,b,_", "randsc", 2);
    one("randsc(_)", "randsc", 0);
    one("randsc()_", NULL, 0); /* the call is closed and behind the cursor */
    one("_randsc(a)", NULL, 0);

    /* --- a call written inside another -----------------------------------*/
    one("trap(randsc(7,_1),3)", "randsc", 1);
    one("trap(randsc(7,1),_3)", "trap", 1);
    one("trap(randsc(7,1)_,3)", "trap", 0);

    /* --- a complex number is one argument, not two -----------------------*/
    one("randsc(7,{1,_1})", "randsc", 1);
    one("randsc(7,{1,1},_2)", "randsc", 2);
    one("randsc(7,{1,_1", "randsc", 1);

    /* --- spaces, and the case the parser does not care about -------------*/
    one("randsc ( a , _b )", "randsc", 1);
    one("RANDSC(a,_b)", "randsc", 1);
    one("  randsc(_a)", "randsc", 0);

    /* --- brackets that group rather than call ----------------------------*/
    one("(a+_b)*z", NULL, 0);
    one("z*(a+_b)", NULL, 0);

    /* --- a name that is not a function -----------------------------------*/
    one("nosuchthing(_a)", NULL, 0);
    one("z(_a)", NULL, 0);

    /* --- an empty bar, and a cursor at either end ------------------------*/
    one("_", NULL, 0);
    one("z^2+c_", NULL, 0);

    /* --- the argument names come out of the help table -------------------*/
    {
        struct formula_hint h = formula_hint_at("randsc(7,{1,1},{0.5,0.5},6", 26);
        int start = 0, len = 0;
        int got = h.fn && formula_hint_arg_span(h.fn->args, h.arg, &start, &len);
        char what[256];
        sprintf(what, "the fourth argument of randsc is \"%.*s\"", got ? len : 0,
                got ? h.fn->args + start : "");
        check(got && len == 16 && !strncmp(h.fn->args + start,
                                           "[kaleidoscope=1]", 16),
              what);
    }

    /* Past the last argument of a call that takes a fixed number: nothing to
     * point at, and pointing at the last one would be a lie. */
    {
        int start = 0, len = 0;
        check(!formula_hint_arg_span("a, b", 2, &start, &len),
              "a third argument of a call that takes two is nothing");
        check(formula_hint_arg_span("a, b", 1, &start, &len) && len == 1 &&
                  start == 3,
              "and the second is where it is written");
    }

    /* A call that ends in "..." goes on taking them, so the tail stays what is
     * current however many are written: poly is "z, k1, k2, ...". */
    {
        int start = 0, len = 0;
        int got = formula_hint_arg_span("z, k1, k2, ...", 7, &start, &len);
        check(got && len == 3 && !strncmp("z, k1, k2, ..." + start, "...", 3),
              "the eighth coefficient of poly is the tail");
    }

    /* --- going on to a new line ------------------------------------------
     *
     * Every character one unit wide, so that a line of avail units holds
     * avail characters and each case can be read off the text. */
    {
        struct {
            const char *text;
            double avail;
            const char *lines; /* the lines, parted by | */
            const char *what;
        } wraps[] = {
            {"z^2+c", 20, "z^2+c", "a formula that fits stays on one line"},
            {"", 20, "", "an empty formula is one empty line"},
            {"sin(z)*cos(z)+tan(z)", 8, "sin(z)*|cos(z)+|tan(z)",
             "a line ends between names, as late as one fits"},
            {"parchment(z,6)", 6, "parchm|ent(z,|6)",
             "a name wider than the whole bar is cut, as late as it fits"},
            {"randsc(13,{0.6,0.6},{1,1})", 11, "randsc(13,{|0.6,0.6},{1|,1})",
             "a number is a name too: 0.6 is never cut at the point"},
            {"poly(z, 1, 2, c)", 10, "poly(z, 1,| 2, c)",
             "a space is between names like anything else that is not one"},
            {"z_p3+p12", 4, "z_p3|+p12",
             "a suffixed variable is one name: z_p3 is not cut at the "
             "underscore"},
            {"abc", 0.5, "a|b|c",
             "a line narrower than any character still takes one"},
        };
        for (const auto &w : wraps) {
            std::u16string text;
            for (const char *c = w.text; *c; c++)
                text += (char16_t)*c;
            std::vector<int> ends = formula_wrap(
                text.data(), (int)text.size(), w.avail,
                [](int start, int end) { return (double)(end - start); });
            std::string got;
            int start = 0;
            for (size_t k = 0; k < ends.size(); k++) {
                if (k)
                    got += '|';
                got.append(w.text + start, ends[k] - start);
                start = ends[k];
            }
            char what[256];
            snprintf(what, sizeof what, "%s (\"%s\")", w.what, got.c_str());
            check(got == w.lines, what);
        }

        /* The two halves of a character outside the basic plane stay
         * together, whatever the width. */
        std::u16string pair = u"a\U0001D49Eb";
        std::vector<int> ends = formula_wrap(
            pair.data(), (int)pair.size(), 1,
            [](int start, int end) { return (double)(end - start); });
        bool whole = true;
        for (int e : ends)
            if (e == 2)
                whole = false;
        check(whole, "a surrogate pair is never cut in two");
    }

    /* --- the history window ----------------------------------------------
     *
     * Again one unit a character. The suffix is what the window gives: the
     * three dots and the length. */
    {
        auto units = [](const std::u16string &t) { return (double)t.size(); };
        auto label = [&units](const char *text, double avail) {
            std::u16string formula, suffix;
            for (const char *c = text; *c; c++)
                formula += (char16_t)*c;
            std::string length = std::to_string(formula.size());
            for (char c : "... (" + length + " characters)")
                suffix += (char16_t)c;
            std::u16string got =
                formula_history_label(formula, suffix, avail, units);
            return std::string(got.begin(), got.end());
        };
        check(label("z^2+c", 40) == "z^2+c", "a formula that fits is listed whole");
        check(label("z^2+c", 5) == "z^2+c", "exactly as wide as the room still fits");
        std::string got = label("randsc(13,{0.6,0.6},{1,1})+sin(z)*cos(z)+c", 30);
        char what[256];
        snprintf(what, sizeof what,
                 "one too long is its start, three dots and its length, and "
                 "no wider than the room (\"%s\")",
                 got.c_str());
        check(got == "randsc(13,{... (42 characters)", what);
        check(label("randsc(13,{0.6,0.6},{1,1})+sin(z)*cos(z)+c", 10) ==
                  "... (42 characters)",
              "with room for less than the suffix, the suffix alone");

        /* Between names: room for "poly(z, {0." ends before the number,
         * rather than in "0...." */
        check(label("poly(z, {0.5,0.25}, {0.125,0.0625})+c", 30) ==
                  "poly(z, {... (37 characters)",
              "the start ends before a number it would cut");
        check(label("poly(z, {0.5,0.25}, {0.125,0.0625})+c", 27) ==
                  "poly(z,... (37 characters)",
              "and without the space before the cut");
        check(label("randsctileabcdefghijklmnopqrstuv", 25) ==
                  "randsc... (32 characters)",
              "a name wider than the whole room is cut");

        std::u16string pair = u"ab\U0001D49Ecde";
        std::u16string cut = formula_history_label(pair, u"...", 6, units);
        check(cut == u"ab...",
              "the start is never cut between the halves of a surrogate pair");
    }

    /* --- the status ------------------------------------------------------
     *
     * One unit a character again, and "..." for the suffix: the length the
     * status puts after it is the window's business, not the rule's. */
    {
        auto u16 = [](const char *text) {
            std::u16string s;
            for (const char *c = text; *c; c++)
                s += (char16_t)*c;
            return s;
        };
        auto units = [](const std::u16string &t) { return (double)t.size(); };
        auto status = [&](const char *formula, double avail) {
            std::vector<std::string> got;
            for (const std::u16string &line :
                 formula_status_lines(u16(formula), u"...", avail, units))
                got.push_back(std::string(line.begin(), line.end()));
            return got;
        };

        const std::u16string text = u16("sin(z)+cos(z)+tan(z)+c");
        auto count = [](int start, int end) { return (double)(end - start); };
        std::vector<int> all =
            formula_wrap(text.data(), (int)text.size(), 8, count);
        std::vector<int> two =
            formula_wrap(text.data(), (int)text.size(), 8, count, 2);
        check(all.size() == 3 && two.size() == 2 && two[0] == all[0] &&
                  two[1] == all[1],
              "most stops after that many lines, each ending where it would");

        check(status("z^2+c", 10) == std::vector<std::string>{"z^2+c"},
              "a formula that fits is one line");
        check(status("sin(z)+cos(z)+c", 8) ==
                  std::vector<std::string>{"sin(z)+", "cos(z)+c"},
              "one too long for a line goes on to a second, between names");
        check(status("sin(z)+cos(z)+tan(z)+c", 8) ==
                  std::vector<std::string>{"sin(z)+", "cos(z..."},
              "one too long for two ends the second with the suffix");
    }

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures != 0;
}
