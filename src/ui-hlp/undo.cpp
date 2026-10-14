#include <cstdlib>
#include <cstring>

#include "undo.h"

/* The states after the current one, which a new step makes unreachable. */
static void cut_redo(struct undo_history *h)
{
    while (h->count > h->current + 1)
        free(h->states[--h->count]);
}

int undo_record(struct undo_history *h, char *state)
{
    if (h->count && !strcmp(state, h->states[h->current])) {
        free(state);
        return 0;
    }
    if (h->count) {
        cut_redo(h);
        if (h->count == UNDOLEVEL) {
            free(h->states[0]);
            memmove(h->states, h->states + 1,
                    (UNDOLEVEL - 1) * sizeof(h->states[0]));
            h->count--;
        }
    }
    h->states[h->count] = state;
    h->current = h->count++;
    return 1;
}

const char *undo_back(struct undo_history *h)
{
    if (!h->count || h->current == 0)
        return NULL;
    return h->states[--h->current];
}

const char *undo_forward(struct undo_history *h)
{
    if (!undo_canredo(h))
        return NULL;
    return h->states[++h->current];
}

int undo_canredo(const struct undo_history *h)
{
    return h->count && h->current < h->count - 1;
}

void undo_clear(struct undo_history *h)
{
    while (h->count)
        free(h->states[--h->count]);
    h->current = 0;
}
