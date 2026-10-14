#ifndef UNDO_H
#define UNDO_H

/* The undo history: the states the picture has been left in, each as the
 * text of a saved position, oldest first, and which of them it is in now.
 *
 * A state is kept only if it differs from the one the picture is in, so
 * recording the same picture twice -- a command that changes nothing, a click
 * that moves nothing -- makes no step to undo. Undo and redo move along the
 * history without changing it; a new state recorded after an undo cuts off
 * what could have been redone, as in an editor. When the history is full the
 * oldest state goes.
 *
 * A history that is all zeros is empty, so one inside a context from calloc
 * needs no setting up. See ui_helper.cpp for when states are recorded. */

#define UNDOLEVEL 256

struct undo_history {
    char *states[UNDOLEVEL];
    int count;   /* states kept */
    int current; /* the one the picture is in, when count is not nought */
};

/* Takes the state, which comes from malloc. 1 if it is a new step; 0 if it is
 * the state the picture is in already, and then it is freed. */
int undo_record(struct undo_history *h, char *state);

/* The state undo or redo goes to, which becomes the current one; NULL, and
 * nothing moves, at either end. The text stays the history's. */
const char *undo_back(struct undo_history *h);
const char *undo_forward(struct undo_history *h);

/* Whether something has been undone that redo could bring back. */
int undo_canredo(const struct undo_history *h);

void undo_clear(struct undo_history *h);

#endif
