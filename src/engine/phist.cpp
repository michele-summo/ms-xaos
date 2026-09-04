/* The p variables: the history of z a user formula may reach back into.
 * See phist.h for what they are and why they are kept this way.
 */
#include <cstdlib>
#include <cstdio>

#include "config.h"

#ifdef USE_SFFE

#include "phist.h"

thread_local struct sffe_phistory sffe_ph = {};

/* The index in a name, or zero if the name is not one of ours. "p" alone is
 * p1; anything past NUM_P_MAX is refused rather than silently taken for
 * something else, and a leading zero is not how one writes a number. */
unsigned int sffe_pindex(const char *name)
{
    if (name[0] != 'p') {
        return 0;
    }
    if (name[1] == '\0') {
        return 1;
    }
    if (name[1] == '0') {
        return 0;
    }
    unsigned int k = 0;
    for (const char *d = name + 1; *d; d++) {
        if (*d < '0' || *d > '9') {
            return 0;
        }
        k = k * 10 + (unsigned int)(*d - '0');
        if (k > NUM_P_MAX) {
            return 0;
        }
    }
    return k;
}

/* Asked by the parser for a name nothing registered answers to. The place is
 * the parser's own, made once and kept until the parser is freed, so it
 * outlives this parse and every later one.
 *
 * Registered under the name the index spells rather than the one written, so
 * that "p" and "p1" are one place and not two saying the same thing. */
sfNumber *sffe_resolve_p(sffe *parser, const char *name)
{
    unsigned int k = sffe_pindex(name);
    if (!k) {
        return NULL;
    }
    char canonical[8];
    snprintf(canonical, sizeof(canonical), "p%u", k);
    sffe *self = parser;
    sfvariable *var = sffe_regvar(&self, NULL, canonical);
    return var ? var->value : NULL;
}

/* Works out which places a parsed formula reads and how deep the ring must be
 * to feed them. Called for both formulas, the second adding to the first. */
static void sffe_ptaps_add(sffe *parser)
{
    struct sffe_phistory *h = &sffe_ph;
    if (parser == NULL) {
        return;
    }
    for (unsigned int i = 0; i < parser->varCount; i += 1) {
        unsigned int back = sffe_pindex(parser->variables[i].name);
        if (!back || !sffe_reads(parser, parser->variables[i].value)) {
            continue;
        }
        /* two parsers may name the same depth, and each has its own place */
        unsigned int k = 0;
        while (k < h->tapcount &&
               h->taps[k].slot != parser->variables[i].value) {
            k += 1;
        }
        if (k == h->tapcount) {
            struct sffe_ptap *grown = (struct sffe_ptap *)realloc(
                h->taps, (h->tapcount + 1) * sizeof(struct sffe_ptap));
            if (!grown) {
                h->tapcount = 0;
                h->depth = 0;
                return;
            }
            h->taps = grown;
            h->taps[h->tapcount].slot = parser->variables[i].value;
            h->taps[h->tapcount].back = back;
            h->tapcount += 1;
        }
        if (back > h->depth) {
            h->depth = back;
        }
    }
}

/* Both formulas have been parsed: work out the taps afresh, since a formula
 * that no longer names p5 must stop paying for it. */
void sffe_ptaps_build(sffe *formula, sffe *initial)
{
    struct sffe_phistory *h = &sffe_ph;
    h->tapcount = 0;
    h->depth = 0;
    sffe_ptaps_add(formula);
    sffe_ptaps_add(initial);
    cmplx *grown =
        (cmplx *)realloc(h->ring, (h->depth ? h->depth : 1) * sizeof(cmplx));
    if (grown) {
        h->ring = grown;
    } else {
        /* no history means no places to fill from it */
        h->depth = 0;
        h->tapcount = 0;
    }
    h->head = 0;
    h->filled = 0;
    cmplxset(h->start, 0, 0);
}

#endif /* USE_SFFE */
