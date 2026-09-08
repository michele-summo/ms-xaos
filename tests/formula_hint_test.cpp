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

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures != 0;
}
