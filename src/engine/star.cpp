#include "config.h"
#include "config.h"
#include <climits>
#include <cstdlib>
#include "filter.h"
#include "xthread.h"

/* The stars are drawn from one random sequence, a number for every pixel and
 * a second for every pixel the fractal reached -- which is what one thread
 * drawing the rows in order takes from it.
 *
 * The threads used to share that sequence: each reseeded it from rand() and
 * they all drew from it at once, so the stars landed differently every time a
 * picture was drawn with more than one. Now each row starts from the place in
 * the sequence that one thread would have reached at it -- counted first, row
 * by row, and jumped to -- and draws from a state of its own: the same stars
 * whatever the number of threads, and with one the stars they always were. */
struct starfielddata {
    struct palette *palette;
    struct palette *savedpalette;
    /* For the picture being drawn: how many numbers each row takes, and then
     * the state each row starts from. */
    unsigned int *starts;
    int rows;
};

#define MYLONG_MAX 0xffffff /*this is enough for me. */
#define STAR_A 1103515245u
#define STAR_C 12345u
static inline unsigned int myrandom(unsigned int &state)
{
    state = ((state * STAR_A) + STAR_C) & MYLONG_MAX;
    return state;
}

/* The state k numbers on from state: the step is x -> A x + C, so k of them
 * are one such map, which squaring finds in log2 k steps. Worked modulo 2^32
 * and cut to 24 bits at the end, which is what cutting after every step
 * gives, the low bits of a product depending only on the low bits of its
 * factors. */
static unsigned int star_jump(unsigned int state, unsigned int k)
{
    unsigned int a = STAR_A, c = STAR_C, ka = 1, kc = 0;
    while (k) {
        if (k & 1) {
            ka = ka * a;
            kc = kc * a + c;
        }
        c = c * a + c;
        a = a * a;
        k >>= 1;
    }
    return (state * ka + kc) & MYLONG_MAX;
}

/* How many numbers each row of the picture takes: one a pixel, and one more
 * for each the fractal reached. */
static void starcount(void *data, struct taskinfo * /*task*/, int r1, int r2)
{
    struct filter *f = (struct filter *)data;
    struct starfielddata *s = (struct starfielddata *)f->data;
    for (int y = r1; y < r2; y++) {
        const pixel8_t *src = f->childimage->currlines[y];
        unsigned int n = (unsigned int)f->childimage->width;
        for (int x = 0; x < f->childimage->width; x++)
            n += src[x] != 0;
        s->starts[y] = n;
    }
}

#define IMAGETYPE SMALLITER
#include "c256.h"
#define do_starfield do_starfield8
#include "stard.h"
#include "hicolor.h"
#define do_starfield do_starfield16
#include "stard.h"
#include "true24.h"
#define do_starfield do_starfield24
#include "stard.h"
#include "truecolor.h"
#define do_starfield do_starfield32
#include "stard.h"
static int requirement(struct filter *f, struct requirements *r)
{
    f->req = *r;
    r->nimages = 1;
    r->flags &= ~IMAGEDATA;
    r->supportedmask = C256 | TRUECOLOR24 | TRUECOLOR | TRUECOLOR16 | GRAYSCALE;
    return (f->next->action->requirement(f->next, r));
}

static int initialize(struct filter *f, struct initdata *i)
{
    struct starfielddata *s = (struct starfielddata *)f->data;
    inhermisc(f, i);
    if (s->savedpalette == NULL)
        s->savedpalette = clonepalette(i->image->palette);
    mkstarfieldpalette(i->image->palette);
    if (!inherimage(f, i, TOUCHIMAGE, 0, 0, s->palette, 0, 0)) {
        return 0;
    }
    setfractalpalette(f, s->savedpalette);
    return (f->previous->action->initialize(f->previous, i));
}

static struct filter *getinstance(const struct filteraction *a)
{
    struct filter *f = createfilter(a);
    struct starfielddata *i = (struct starfielddata *)calloc(1, sizeof(*i));
    i->savedpalette = NULL;
    i->starts = NULL;
    i->rows = 0;
    i->palette = createpalette(0, 65536, IMAGETYPE, 0, 65536, NULL, NULL, NULL,
                               NULL, NULL);
    f->data = i;
    f->name = QT_TRANSLATE_NOOP("Filter", "Starfield");
    return (f);
}

static void destroyinstance(struct filter *f)
{
    struct starfielddata *i = (struct starfielddata *)f->data;
    if (i->savedpalette != NULL)
        destroypalette(i->savedpalette);
    destroypalette(i->palette);
    free(i->starts);
    destroyinheredimage(f);
    free(f);
    free(i);
}

static int doit(struct filter *f, int flags, int time)
{
    struct starfielddata *s = (struct starfielddata *)f->data;
    int val;
    unsigned int state;
    val = f->previous->action->doit(f->previous, flags, time);
    if (s->rows != f->image->height) {
        free(s->starts);
        s->starts = (unsigned int *)malloc(sizeof(*s->starts) *
                                           (size_t)f->image->height);
        if (s->starts == NULL) {
            s->rows = 0;
            return val;
        }
        s->rows = f->image->height;
    }
    /* Seeded once a picture, from the thread that always seeded it when
     * there was one: the picture after this one is seeded as it was too. */
    state = (unsigned int)rand();
    xth_function(starcount, f, f->image->height);
    xth_sync();
    for (int y = 0; y < f->image->height; y++) {
        unsigned int n = s->starts[y];
        s->starts[y] = state;
        state = star_jump(state, n);
    }
    drivercall(*f->image, xth_function(do_starfield8, f, f->image->height),
               xth_function(do_starfield16, f, f->image->height),
               xth_function(do_starfield24, f, f->image->height),
               xth_function(do_starfield32, f, f->image->height));
    xth_sync();
    return val | CHANGED;
}

static void myremovefilter(struct filter *f)
{
    struct starfielddata *s = (struct starfielddata *)f->data;
    if (s->savedpalette != NULL) {
        restorepalette(f->image->palette, s->savedpalette);
        destroypalette(s->savedpalette);
        s->savedpalette = NULL;
    }
}

const struct filteraction starfield_filter = {
    "Starfield",        "starfield",     0,
    getinstance,        destroyinstance, doit,
    requirement,        initialize,      convertupgeneric,
    convertdowngeneric, myremovefilter};
