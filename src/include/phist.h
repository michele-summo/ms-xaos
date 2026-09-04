/* The p variables: how far back a user formula may look.
 *
 * Split out of formulas.cpp so that the few lines the iteration loops run per
 * pass can be inlined into them while a test can still drive them, which is
 * the only way to say plainly that what replaced the shift reads back the way
 * the shift did.
 */
#ifndef PHIST_H
#define PHIST_H

#include "config.h"

#ifdef USE_SFFE

#include "sffe.h"
#include "cmplx.h"

/* --- the p variables: how far back a formula may look ---------------------
 *
 * p1 is the value z had on the pass before this one, p2 the one before that,
 * and so on; p is another name for p1. They used to stop at p6 and were kept
 * by shifting a six-place array along by one every pass, which is the obvious
 * way to keep a short history and the wrong way to keep a long one: the shift
 * costs a copy per place per pass whether or not the formula names any of
 * them, so p9999 would have cost 320 KB of copying per pass -- twenty-four
 * times the whole cost of iterating z^2+c, measured.
 *
 * So the history is a ring, written once per pass, and only the places the
 * formula actually names are read out of it. sffe hands out a place for a name
 * when it first meets it, so only those exist; the ring is only as deep as the
 * furthest one reaches; and a pass costs one store plus one copy per name
 * written, whatever the depth. Flat from p1 to p9999, measured, and a formula
 * that names no p at all pays nothing at all, where the shift used to take an
 * eighth of its time.
 *
 * All of it in one structure because all of it is thread-local: reached one
 * member at a time it was seven thread-local lookups and three pointers
 * chased per pass, which cost more than the shift it replaced. Reached
 * through one address it costs less.
 */
struct sffe_ptap {
    sfNumber *slot;    /* what the parser bound the name to */
    unsigned int back; /* how many passes back it stands for, p1 being one */
};

struct sffe_phistory {
    cmplx *ring;                /* depth places, newest written at head - 1 */
    struct sffe_ptap *taps;     /* the places the formulas name */
    unsigned int depth;         /* as far back as the deepest of them looks */
    unsigned int head;          /* where the next pass will be written */
    unsigned int tapcount;
    /* What the history stands at before anything has been computed, and how
     * much of the ring has been written since. A place reaching further back
     * than that has nothing behind it yet and stands at the start value --
     * which is what makes starting a pixel cost nothing, where filling the
     * ring with the start value cost a store per place and was, at nine
     * thousand of them, six times the whole cost of the pixel. */
    cmplx start;
    unsigned int filled;
};

extern thread_local struct sffe_phistory sffe_ph;

/* The index in a name, or zero if the name is not one of ours. */
unsigned int sffe_pindex(const char *name);

/* Hands the parser a place for p, p1, p2 ... on sight, so that a formula may
 * name any of them without all of them being registered first. */
sfNumber *sffe_resolve_p(sffe *parser, const char *name);

/* Which places the two parsed formulas read, and how deep a history feeds
 * them. Either may be null. Once per formula, not once per pass. */
void sffe_ptaps_build(sffe *formula, sffe *initial);

/* Fills every place the formulas name from the ring. p1 is the value z had
 * going into the pass before this one, and the ring holds the last depth of
 * those with the newest at head - 1. */
static inline void sffe_ptaps_read(struct sffe_phistory *h)
{
    for (unsigned int k = 0; k < h->tapcount; k += 1) {
        unsigned int back = h->taps[k].back;
        if (back > h->filled) {
            *h->taps[k].slot = h->start;
            continue;
        }
        unsigned int at =
            h->head >= back ? h->head - back : h->head + h->depth - back;
        *h->taps[k].slot = h->ring[at];
    }
}

/* A new pixel: nothing has been computed yet, so every place stands for the
 * value the whole history starts at -- the point itself, or zero, as the menu
 * says. */
static inline void sffe_phist_reset(number_t re, number_t im)
{
    struct sffe_phistory *h = &sffe_ph;
    cmplxset(h->start, re, im);
    h->head = 0;
    h->filled = 0;
    sffe_ptaps_read(h);
}

/* One pass done: what z was going into it becomes p1, and everything the
 * formula names moves back one. */
static inline void sffe_phist_push(number_t re, number_t im)
{
    struct sffe_phistory *h = &sffe_ph;
    if (!h->depth) {
        return;
    }
    cmplxset(h->ring[h->head], re, im);
    h->head = h->head + 1 == h->depth ? 0 : h->head + 1;
    if (h->filled < h->depth) {
        h->filled += 1;
    }
    sffe_ptaps_read(h);
}

#endif /* USE_SFFE */
#endif /* PHIST_H */
