/* The undo history in undo.cpp.
 *
 * Undo used to keep a ring of positions and go back two places in it, on the
 * understanding that the last one had just been recorded by the undo command
 * itself. Anything that recorded a position at another moment -- the mouse,
 * at the release, while the zoom was still slowing down -- or that recorded
 * one without changing anything -- a recalculation -- put undo one place out:
 * it went back to a picture half way through a zoom, or to the one it was
 * already showing.
 *
 * These checks pin what the history now promises: a step only where the state
 * differs, undo and redo moving one step whatever was recorded in between, and
 * redo surviving the state undo went to being recorded again.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "undo.h"

static int failures = 0;

static void check(int ok, const char *what)
{
    printf("%-6s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        failures++;
}

static char *state(const char *text)
{
    char *s = (char *)malloc(strlen(text) + 1);
    strcpy(s, text);
    return s;
}

static int is(const char *got, const char *want)
{
    return got != NULL && !strcmp(got, want);
}

int main()
{
    static struct undo_history h; /* all zeros: empty */
    char name[32];

    check(undo_back(&h) == NULL && undo_forward(&h) == NULL &&
              !undo_canredo(&h),
          "an empty history goes nowhere");

    check(undo_record(&h, state("A")) == 1, "the first state is a step");
    check(undo_record(&h, state("A")) == 0 && h.count == 1,
          "the same state again is not");
    check(undo_back(&h) == NULL && h.current == 0,
          "with one state there is nothing to undo");

    undo_record(&h, state("B"));
    undo_record(&h, state("B"));
    undo_record(&h, state("C"));
    check(h.count == 3, "three states, the repeat left out");
    check(is(undo_back(&h), "B"), "undo goes to the state before");
    check(is(undo_back(&h), "A"), "and to the one before that");
    check(undo_back(&h) == NULL && h.current == 0,
          "and stops at the first");
    check(is(undo_forward(&h), "B") && is(undo_forward(&h), "C"),
          "redo comes back one state at a time");
    check(undo_forward(&h) == NULL && !undo_canredo(&h),
          "and stops at the last");

    /* What undo went to is recorded again by the next command that changes
     * nothing, and by the next click that moves nothing: redo must not be
     * lost to it. */
    undo_back(&h);
    check(undo_record(&h, state("B")) == 0 && undo_canredo(&h) &&
              is(undo_forward(&h), "C"),
          "recording the state undo went to keeps redo");

    undo_back(&h);
    undo_back(&h);
    check(undo_record(&h, state("D")) == 1 && !undo_canredo(&h) &&
              h.count == 2,
          "a new state after an undo cuts off redo");
    check(is(undo_back(&h), "A") && is(undo_forward(&h), "D"),
          "and follows the state it was recorded from");

    undo_clear(&h);
    check(h.count == 0 && undo_back(&h) == NULL, "a cleared history is empty");
    check(undo_record(&h, state("A")) == 1 && h.current == 0,
          "and starts again");

    for (int i = 1; i < UNDOLEVEL + 10; i++) {
        sprintf(name, "S%d", i);
        undo_record(&h, state(name));
    }
    sprintf(name, "S%d", UNDOLEVEL + 9);
    check(h.count == UNDOLEVEL && is(h.states[h.current], name),
          "a full history keeps the newest states");
    int back = 0;
    while (undo_back(&h))
        back++;
    check(back == UNDOLEVEL - 1 && is(h.states[0], "S10"),
          "and loses the oldest");

    undo_clear(&h);
    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures != 0;
}
