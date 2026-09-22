#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <climits>
#include "config.h"
#include "filter.h"
#include "misc-f.h"
#include "xerror.h"

#define nprecells                                                              \
    (context->type &                                                           \
             (LARGEITER | SMALLITER | TRUECOLOR | TRUECOLOR16 | TRUECOLOR24)   \
         ? 0                                                                   \
         : (context)->ncells)
#define PREALLOCATED(palette)                                                  \
    ((palette)->type &                                                         \
             (LARGEITER | SMALLITER | TRUECOLOR | TRUECOLOR16 | TRUECOLOR24)   \
         ? 0                                                                   \
         : (palette)->npreallocated)
/*emulate allocation routines using setpalette */
unsigned col_diff[3][512];
static struct palette *context;
static int maxentries;
static int needupdate;
static void genertruecolorinfo(struct truec *t, unsigned int r, unsigned int g,
                               unsigned int b)
{
    int n;

    for (n = 0; !((r >> n) & 1); n++)
        ;
    t->rshift = n;
    for (; ((r >> n) & 1); n++)
        ;
    t->rprec = 8 - (n - t->rshift);
    t->rmask = r;

    for (n = 0; !((g >> n) & 1); n++)
        ;
    t->gshift = n;
    for (; ((g >> n) & 1); n++)
        ;
    t->gprec = 8 - (n - t->gshift);
    t->gmask = g;

    for (n = 0; !((b >> n) & 1); n++)
        ;
    t->bshift = n;
    for (; ((b >> n) & 1); n++)
        ;
    t->bprec = 8 - (n - t->bshift);
    t->bmask = b;

    t->allmask = r | g | b;
    if ((r & b) || (r & g) || (b & g)) {
        x_fatalerror("Internal error:Invalid color masks 1 %x %x %x!\n", r, g,
                     b);
    }
    if ((r < g && g < b) || (b < g && g < r)) {
        t->mask1 = r | b;
        t->mask2 = g;
    } else if ((g < r && r < b) || (b < r && r < g)) {
        t->mask1 = g | b;
        t->mask2 = r;
    } else if ((g < b && b < r) || (r < b && b < g)) {
        t->mask1 = g | r;
        t->mask2 = b;
    }
    t->byteexact = 0;
    t->missingbyte = -1;
    if (!(t->rshift % 8) && !(t->gshift % 8) && !(t->bshift % 8) && !t->rprec &&
        !t->gprec && !t->bprec) {
        t->byteexact = 1;
        {
            unsigned char ch[4];
            *(unsigned int *)ch = t->allmask;
            if (!ch[0])
                t->missingbyte = 0;
            if (!ch[1])
                t->missingbyte = 1;
            if (!ch[2])
                t->missingbyte = 2;
            if (!ch[3])
                t->missingbyte = 3;
        }
    }
#ifdef DEBUG
    printf("Image:\n");
    printf("rshift:%i gshift:%i bshift:%i\n", t->rshift, t->gshift, t->bshift);
    printf("rprec:%i gprec:%i bprec:%i\n", t->rprec, t->gprec, t->bprec);
    printf("rmask:%x gmask:%x bmask:%x\n", t->rmask, t->gmask, t->bmask);
    printf("mask1:%x mask2:%x allmask:%x\n", t->mask1, t->mask2, t->allmask);
    printf("byteexact:%x missingbyte:%i\n", t->byteexact, t->missingbyte);
#endif
}

void bestfit_init(void)
{
    int i;

    for (i = 1; i < 256; i++) {
        int k = i * i;
        col_diff[0][i] = col_diff[0][512 - i] = k * (59 * 59) / 256;
        col_diff[1][i] = col_diff[1][512 - i] = k * (30 * 30) / 256;
        col_diff[2][i] = col_diff[2][512 - i] = k * (11 * 11) / 256;
    }
}

static int allocgenerictruecolor(struct palette *palette, int init, int r,
                                 int g, int b)
{
    unsigned int n;
    switch (palette->type) {
        case LARGEITER:
        case SMALLITER:
            return 1;
        case TRUECOLOR:
        case TRUECOLOR16:
        case TRUECOLOR24:
        default:
            n = ((r >> palette->info.truec.rprec)
                 << palette->info.truec.rshift) |
                ((g >> palette->info.truec.gprec)
                 << palette->info.truec.gshift) |
                ((b >> palette->info.truec.bprec)
                 << palette->info.truec.bshift) |
                palette->info.truec.alpha;
            break;
    }
    if (init)
        palette->size = 0;
    else if (palette->size >= palette->maxentries)
        return -1;
    palette->pixels[palette->size] = n;
    palette->size++;
    return palette->size;
}

static int allocgeneric(struct palette *palette, int init, int r, int g, int b)
{
    int start = palette->npreallocated + palette->start;
    if (init)
        palette->size = 0;
    else if (palette->size >= palette->end - start)
        return -1;
    palette->pixels[palette->size] = palette->size + start;
    palette->rgb[palette->size + start][0] = r;
    palette->rgb[palette->size + start][1] = g;
    palette->rgb[palette->size + start][2] = b;
    palette->size++;
    return (palette->size - 1);
}

int fixedalloccolor(struct palette *palette, int init, int r, int g, int b)
{
    int i, coldif, lowest, bestfit;
    if (init)
        palette->size = 0;
    else if (palette->size >= palette->maxentries)
        return -1;
    lowest = INT_MAX;
    bestfit = 1;
    if (palette->type == FIXEDCOLOR || (palette->type & BITMAPS)) {
        for (i = palette->start; i < palette->end; i++) {
            coldif = col_diff[0][(g - palette->rgb[i][1]) & 0x1ff];
            if (coldif < lowest) {
                coldif += col_diff[1][(r - palette->rgb[i][0]) & 0x1ff];
                if (coldif < lowest) {
                    coldif += col_diff[2][(b - palette->rgb[i][2]) & 0x1ff];
                    if (coldif < lowest) {
                        bestfit = i;
                        if (!coldif)
                            break;
                        lowest = coldif;
                    }
                }
            }
        }
    } else {
        bestfit = (r * 30 + g * 59 + b * 11) * (palette->end - palette->start) /
                      256 / 100 +
                  palette->start;
    }
    palette->pixels[palette->size] = bestfit;
    palette->size++;
    return (palette->size - 1);
}

static void setcolorgeneric(struct palette * /*palette*/, int /*start*/,
                            int /*end*/, rgb_t * /*rgb*/)
{
}

static void allocfinishedgeneric(struct palette *palette)
{
    palette->setpalette(palette, palette->start,
                        palette->size + palette->start + palette->npreallocated,
                        palette->rgb + palette->start);
}

static void cycle_entries(struct palette *c, int direction)
{
    int i;
    int i1, i2, i3;
    rgb_t *co;
    if (direction > 0)
        direction %= c->size - 1;
    else
        direction = -((-direction) % (c->size - 1));
    if (!direction)
        return;
    co = (rgb_t *)malloc(c->end * sizeof(rgb_t));
    memcpy(co, c->rgb, sizeof(*co) * c->end);
    i3 = (c->size - 1 + direction) % (c->size - 1) + 1;
    for (i = 1; i < c->size; i++) {
        i1 = c->pixels[i];
        i2 = c->pixels[i3];
        c->rgb[i1][0] = co[i2][0];
        c->rgb[i1][1] = co[i2][1];
        c->rgb[i1][2] = co[i2][2];
        i3++;
        if (i3 >= c->size)
            i3 = 1;
    }
    free(co);
}

static void cyclecolorsgeneric(struct palette *pal, int direction)
{
    cycle_entries(pal, direction);
    pal->setpalette(pal, pal->pixels[0], pal->size + pal->pixels[0],
                    pal->rgb + pal->pixels[0]);
}

#define TRUECOLORPALETTE 65536
struct palette *createpalette(
    int start, int end, int type, int flags, int maxentries,
    int (*alloccolor)(struct palette *pal, int init, int r, int g, int b),
    void (*setcolor)(struct palette *pal, int start, int end, rgb_t *rgb),
    void (*allocfinished)(struct palette *pal),
    void (*cyclecolors)(struct palette *pal, int direction),
    union paletteinfo *info)
{
    static int versioncount;
    struct palette *palette =
        (struct palette *)calloc(1, sizeof(struct palette));

    if (col_diff[0][1] == 0)
        bestfit_init();
    palette->ncells = 0;
    palette->index = NULL;
    palette->size = 0;
    palette->rgb = NULL;
    if (palette == NULL)
        return NULL;

    switch (type) {
        case LBITMAP:
        case MBITMAP:
        case LIBITMAP:
        case MIBITMAP:
            end = 2;
            start = 0;
            maxentries = 256;
            palette->rgb = (rgb_t *)calloc(end, sizeof(*palette->rgb));
            if (palette->rgb == NULL) {
                free(palette);
                return NULL;
            }
            if (type & (LIBITMAP | MIBITMAP)) {
                palette->rgb[0][0] = 255;
                palette->rgb[0][1] = 255;
                palette->rgb[0][2] = 255;
                palette->rgb[1][0] = 0;
                palette->rgb[1][1] = 0;
                palette->rgb[1][2] = 0;
            } else {
                palette->rgb[0][0] = 0;
                palette->rgb[0][1] = 0;
                palette->rgb[0][2] = 0;
                palette->rgb[1][0] = 255;
                palette->rgb[1][1] = 255;
                palette->rgb[1][2] = 255;
            }
            palette->maxentries = maxentries;
            palette->alloccolor = fixedalloccolor;
            palette->setpalette = NULL;
            palette->allocfinished = NULL;
            palette->cyclecolors = NULL;
            break;
        case FIXEDCOLOR:
            if (!end)
                end = 256;
            if (!maxentries)
                maxentries = 256;
            palette->rgb = (rgb_t *)calloc(end, sizeof(*palette->rgb));
            if (palette->rgb == NULL) {
                free(palette);
                return NULL;
            }
            palette->maxentries = maxentries;
            palette->alloccolor = fixedalloccolor;
            palette->setpalette = NULL;
            palette->allocfinished = NULL;
            palette->cyclecolors = NULL;
            break;
        case GRAYSCALE:
            if (!end)
                end = 256;
            if (!maxentries)
                maxentries = end - start;
            palette->maxentries = 65536;
            palette->alloccolor = fixedalloccolor;
            palette->setpalette = NULL;
            palette->allocfinished = NULL;
            palette->cyclecolors = NULL;
            palette->size = end - start;
            break;
        case C256:

            if (!end)
                end = 256;
            if (!maxentries)
                maxentries = end - start;
            if (cyclecolors == NULL && setcolor != NULL)
                cyclecolors = cyclecolorsgeneric;

            if (alloccolor == NULL) {
                alloccolor = allocgeneric, allocfinished = allocfinishedgeneric;
                if (setcolor == NULL && type == C256) /*non hardware palette */
                    setcolor = setcolorgeneric,
                    cyclecolors = cyclecolorsgeneric;
            }
            palette->rgb = (rgb_t *)calloc(end, sizeof(*palette->rgb));

            if (palette->rgb == NULL) {
                free(palette);
                return NULL;
            }

            palette->maxentries = maxentries;
            palette->alloccolor = alloccolor;
            palette->setpalette = setcolor;
            palette->allocfinished = allocfinished;
            palette->cyclecolors = cyclecolors;

            break;
        default:
            end = TRUECOLORPALETTE;
            start = 0;
            if (type == SMALLITER)
                end = 256;
            start = 0;
            palette->maxentries = end;
            palette->alloccolor = allocgenerictruecolor;
            palette->cyclecolors = NULL;
            palette->setpalette = NULL;
            palette->allocfinished = NULL;
    }
    {
        int ee = palette->maxentries;
        /*if(end>palette->maxentries) ee=end; */
        if (ee < 256)
            palette->pixels =
                (unsigned int *)calloc(256, sizeof(*palette->pixels));
        else
            palette->pixels =
                (unsigned int *)calloc(ee, sizeof(*palette->pixels));
    }
    if (palette->pixels == NULL) {
        free(palette);
        return NULL;
    }
    if (type & (LARGEITER | SMALLITER)) {
        int i;
        palette->size = end;
        palette->flags |= DONOTCHANGE;
        for (i = 0; i < end; i++)
            palette->pixels[i] = i;
    }
    palette->start = start;
    palette->end = end;
    palette->type = type;
    palette->flags |= flags;
    palette->version = (versioncount += 65536);
    if (type == FIXEDCOLOR) {
        int i;
        if (setcolor != NULL)
            setcolor(palette, start, end, palette->rgb);
        for (i = 0; i < end - start; i++)
            palette->pixels[i] = i + start;
    }
    if (type == GRAYSCALE) {
        int i;
        for (i = palette->start; i < end - start; i++)
            palette->pixels[i] = i + start;
    }
    if (info != NULL &&
        (type == TRUECOLOR || type == TRUECOLOR24 || type == TRUECOLOR16)) {
        genertruecolorinfo(&palette->info.truec, info->truec.rmask,
                           info->truec.gmask, info->truec.bmask);
    } else {
        if (type == TRUECOLOR)
            genertruecolorinfo(&palette->info.truec, 255 << 16, 255 << 8, 255);
        if (type == TRUECOLOR24)
            genertruecolorinfo(&palette->info.truec, 255 << 16, 255 << 8, 255);
        if (type == TRUECOLOR16)
            genertruecolorinfo(&palette->info.truec, (255 >> 3) << 11,
                               (255 >> 2) << 5, (255 >> 3));
    }
    if (type == TRUECOLOR)
        palette->info.truec.alpha = ~palette->info.truec.allmask;
    else
        palette->info.truec.alpha = 0;

    return (palette);
}

void destroypalette(struct palette *palette)
{
    free(palette->pixels);
    if (palette->rgb != NULL)
        free(palette->rgb);
    if (palette->index != NULL)
        free(palette->index);
    free(palette);
}

#define MYMIN(x, y) ((x) < (y) ? (x) : (y))
struct palette *clonepalette(struct palette *palette)
{
    struct palette *pal =
        createpalette(palette->start, palette->end, palette->type,
                      palette->flags, palette->maxentries,
                      /*palette->alloccolor, palette->setpalette,
                         palette->allocfinished, palette->cyclecolors */
                      NULL, NULL, NULL, NULL, &palette->info);
    memcpy(pal->pixels, palette->pixels,
           sizeof(*pal->pixels) * MYMIN(palette->end, pal->end));
    if (pal->rgb != NULL) {
        memcpy(pal->rgb, palette->rgb,
               sizeof(*pal->rgb) * MYMIN(palette->end, pal->end));
    }
    pal->size = palette->size;
    return (pal);
}

void preallocpalette(struct palette *pal)
{
    int i;
    int p;
    if (pal->index != NULL)
        free(pal->index), pal->index = NULL;
    pal->npreallocated = 0;
    if (!pal->ncells)
        return;
    pal->index = (unsigned int *)malloc(sizeof(int) * (pal->ncells + 1));
    for (i = 0; i < pal->ncells; i++) {
        if (!i)
            p = pal->pixels[0];
        else
            p = pal->pixels[pal->size];
        pal->alloccolor(pal, i == 0, pal->prergb[i][0], pal->prergb[i][1],
                        pal->prergb[i][2]);
        if (pal->size) {
            pal->index[i] = pal->pixels[pal->size - 1];
            pal->pixels[pal->size - 1] = p;
        }
    }
    pal->npreallocated = pal->size;
    pal->size = 0;
    needupdate = 0;
}

void restorepalette(struct palette *dest, struct palette *src)
{
    int i;
    preallocpalette(dest);
    for (i = 0; i < src->size; i++) {
        int r = 0, g = 0, b = 0;
        switch (src->type) {
            case SMALLITER:
                r = g = b = i;
                break;
            case LARGEITER:
                r = g = b = i / 256;
                break;
            case GRAYSCALE:
                r = g = b = src->pixels[i];
                break;
            case C256:
            case FIXEDCOLOR:
            case MBITMAP:
            case LBITMAP:
            case MIBITMAP:
            case LIBITMAP:
                r = src->rgb[src->pixels[i]][0];
                g = src->rgb[src->pixels[i]][1];
                b = src->rgb[src->pixels[i]][2];
                break;
            case TRUECOLOR:
            case TRUECOLOR16:
            case TRUECOLOR24:
                r = (((src->pixels[i] & src->info.truec.rmask) >>
                      src->info.truec.rshift))
                    << src->info.truec.rprec;
                g = (((src->pixels[i] & src->info.truec.gmask) >>
                      src->info.truec.gshift))
                    << src->info.truec.gprec;
                b = (((src->pixels[i] & src->info.truec.bmask) >>
                      src->info.truec.bshift))
                    << src->info.truec.bprec;
                break;
        }
        if (dest->size >= dest->maxentries - PREALLOCATED(dest))
            break;
        if (dest->alloccolor(dest, (i + dest->npreallocated) == 0, r, g, b) ==
            -1)
            break;
    }
    if (!(dest->flags & FINISHLATER)) {
        if (dest->allocfinished != NULL)
            dest->allocfinished(dest);
    } else
        dest->flags |= UNFINISHED;
    dest->version++;
}

/*Generation code for various palettes */

#define DEFNSEGMENTS (255 / 8)
#define MAXNSEGMENTS (4096)
#define NSEGMENTS ((255 / segmentsize))
static int segmentsize;

static unsigned char colors[MAXNSEGMENTS][3];
static const unsigned char colors1[DEFNSEGMENTS][3] = {
    /*{8, 14, 32}, */
    {0, 0, 0},    {120, 119, 238}, {24, 7, 25},  {197, 66, 28},
    {29, 18, 11}, {135, 46, 71},   {24, 27, 13}, {241, 230, 128},
    {17, 31, 24}, {240, 162, 139}, {11, 4, 30},  {106, 87, 189},
    {29, 21, 14}, {12, 140, 118},  {10, 6, 29},  {50, 144, 77},
    {22, 0, 24},  {148, 188, 243}, {4, 32, 7},   {231, 146, 14},
    {10, 13, 20}, {184, 147, 68},  {13, 28, 3},  {169, 248, 152},
    {4, 0, 34},   {62, 83, 48},    {7, 21, 22},  {152, 97, 184},
    {8, 3, 12},   {247, 92, 235},  {31, 32, 16}};

static int allocate(int r, int g, int b, int init)
{
    unsigned int n;
    if (init)
        preallocpalette(context);
    n = context->pixels[(init ? 0 : context->size) /*+ context->start */];
    if (!init && context->size == maxentries)
        return 0;
    if ((context->alloccolor(context, init && !context->npreallocated, (int)r,
                             (int)g, (int)b)) == -1) {
        return 0;
    }
    if (context->pixels[context->size - 1 /*+ context->start */] != n) {
        needupdate = 1;
    }
    return (1);
}

static int mksmooth(int nsegments, int setsegments)
{
    int i, y;
    float r, g, b, rs, gs, bs;
    int segmentsize1 = segmentsize;

    for (i = 0; i < setsegments; i++) {
        if (i == setsegments - 1 && !(context->flags & UNKNOWNENTRIES)) {
            segmentsize1 = maxentries - context->size - 2;
        }
        r = colors[i % nsegments][0];
        g = colors[i % nsegments][1];
        b = colors[i % nsegments][2];
        rs = ((int)colors[(i + 1) % setsegments % nsegments][0] - r) /
             (unsigned int)segmentsize1;
        gs = ((int)colors[(i + 1) % setsegments % nsegments][1] - g) /
             (unsigned int)segmentsize1;
        bs = ((int)colors[(i + 1) % setsegments % nsegments][2] - b) /
             (unsigned int)segmentsize1;
        for (y = 0; y < segmentsize1; y++) {
            if (!allocate((int)r, (int)g, (int)b, i == 0 && y == 0)) {
                if (!i)
                    context->size = 2;
                context->size = i * segmentsize;
                return 0;
            }
            r += rs;
            g += gs;
            b += bs;
        }
    }
    if (context->flags & UNKNOWNENTRIES)
        context->size = i * segmentsize;
    return 1;
}

static inline void hsv_to_rgb(int h, int s, int v, unsigned char *red,
                              unsigned char *green, unsigned char *blue)
{
    int hue;
    int f, p, q, t;
    h += 256;
    h %= 256;

    if (s == 0) {
        *red = v;
        *green = v;
        *blue = v;
    } else {
        h %= 256;
        if (h < 0)
            h += 256;
        hue = h * 6;

        f = hue & 255;
        p = v * (256 - s) / 256;
        q = v * (256 - (s * f) / 256) >> 8;
        t = v * (256 * 256 - (s * (256 - f))) >> 16;

        switch ((int)(hue / 256)) {
            case 0:
                *red = v;
                *green = t;
                *blue = p;
                break;
            case 1:
                *red = q;
                *green = v;
                *blue = p;
                break;
            case 2:
                *red = p;
                *green = v;
                *blue = t;
                break;
            case 3:
                *red = p;
                *green = q;
                *blue = v;
                break;
            case 4:
                *red = t;
                *green = p;
                *blue = v;
                break;
            case 5:
                *red = v;
                *green = p;
                *blue = q;
                break;
        }
    }
}

static void randomize_segments3(int /*whitemode*/, int nsegments)
{
    int i = 0;
    int h, s, v;

    for (i = 0; i < nsegments; i++) {
        if (!(i % 3)) {
            if (i % 6)
                colors[i][0] = 255, colors[i][1] = 255, colors[i][2] = 255;
            else
                colors[i][0] = 0, colors[i][1] = 0, colors[i][2] = 0;
        } else {
            s = (int)XaoS_random() % 256;
            h = (int)XaoS_random() % (128 - 32);
            v = (int)XaoS_random() % 128;
            if (((i) % 6 > 3) ^ ((i) % 3 == 1))
                /*if(((i)%3==1)) */
                h += 42 + 16;
            else
                h += 42 + 128 + 16, v += 128 + 64;
            hsv_to_rgb(h, s, v, colors[i], colors[i] + 1, colors[i] + 2);
        }
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

static void randomize_segments2(int whitemode, int nsegments)
{
    int i = 0;

    for (i = 0; i < nsegments; i++) {
        if (i % 3 == 2)
            colors[i][0] = whitemode * 255,
            colors[i][1] = whitemode * 255,
            colors[i][2] = whitemode * 255;
        else if (i % 3 == 0)
            colors[i][0] = (!whitemode) * 255,
            colors[i][1] = (!whitemode) * 255,
            colors[i][2] = (!whitemode) * 255;
        else
            colors[i][0] = (int)XaoS_random() % 256,
            colors[i][1] = (int)XaoS_random() % 256,
            colors[i][2] = (int)XaoS_random() % 256;
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

static void randomize_segments(int whitemode, int nsegments)
{
    int i = 0;
    if (whitemode) {
        colors[0][0] = 255, colors[0][1] = 255, colors[0][2] = 255;
        for (i = 0; i < nsegments; i += 2) {
            if (i != 0) {
                colors[i][0] = (int)XaoS_random() % 256,
                colors[i][1] = (int)XaoS_random() % 256,
                colors[i][2] = (int)XaoS_random() % 256;
            }
            if (i + 1 < nsegments)
                colors[i + 1][0] = (int)XaoS_random() % 35,
                           colors[i + 1][1] = (int)XaoS_random() % 35,
                           colors[i + 1][2] = (int)XaoS_random() % 35;
        }
    } else {
        for (i = 0; i < nsegments; i += 2) {
            colors[i][0] = (int)XaoS_random() % 35,
            colors[i][1] = (int)XaoS_random() % 35,
            colors[i][2] = (int)XaoS_random() % 35;
            if (i + 1 < nsegments)
                colors[i + 1][0] = (int)XaoS_random() % 256,
                           colors[i + 1][1] = (int)XaoS_random() % 256,
                           colors[i + 1][2] = (int)XaoS_random() % 256;
        }
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}


/* --- four more, taken from what the first three are made of ---------------
 *
 * These four were made after measuring the first three rather than by reasoning
 * about colour: twenty thousand palettes of each, made as the program makes
 * them, every colour named by the nearest of three dozen reference colours in
 * CIELAB, and the colours and the pairs of colours counted.
 *
 * What that found is what they are built on. The first three get their look
 * from anchors -- black and white, a third to a half of all their stops -- with
 * colours between them, a period of two or three stops so that the pattern
 * shows in the first few bands, which are all most pictures use. And 1 and 2
 * turn out to hold only 256 colours: they take the low byte of the generator,
 * which is a generator of period 256 on its own, so red decides green and blue
 * and every colour is always followed by the same one. That is where their
 * recurring pairs come from, pink beside yellow twenty-four times as often as
 * chance, burgundy beside olive six. 3 draws from two windows instead, warm
 * bright colours over dark greens and teals.
 *
 * The four before these were a spectrum, a duotone, a triad and a
 * complementary pair, reasoned out from the colour circle, and they read as
 * one colour: they spread their hues over the whole length of the palette,
 * which is some three thousand stops in the program and four to nine in the
 * test that was meant to catch it, so a picture saw the first hue alone; and
 * they kept every stop a colour, deep or pale, with no anchor, so nothing in
 * them was dark or light enough to give the colours an edge.
 *
 * So each of these decides everything stop by stop, and each was measured
 * against the first three before it stayed: the Jensen-Shannon divergence
 * between the colours two algorithms show, by name, and between the pairs of
 * neighbours they show, over five thousand palettes. Two of them took a
 * skeleton of the older ones and changed only the colours between its anchors,
 * and measured all but the same as the one they took it from, 0.05 and 0.00;
 * those two have no black and white skeleton now. The other two keep one and
 * differ enough in their colours.
 *
 * Everything random here is drawn from the top of the generator's word rather
 * than its bottom byte. Each is driven from the seed the caller set, so an
 * algorithm and a seed give the same palette every time -- which is all a
 * saved position records of its colours. */

/* A number from nought to n - 1 out of the top sixteen bits of the generator,
 * which have a period of their own of millions where the bottom eight have one
 * of 256. */
static int prand(int n)
{
    return (int)(((long)(XaoS_random() >> 8) * n) >> 16);
}

/* 4: smog -- dark, melancholy, and grey the way a city's air is grey.
 *
 * The first attempt here was 1's mechanism on another ring of 256 colours,
 * and however far the ring was chosen from the classic one it looked like 1:
 * measured as the Jensen-Shannon divergence between the colours the two show,
 * by name, it was nearer 1 than 1 is to 2 (0.05 against 0.11). The mechanism
 * was the likeness -- black half the time, a colour from anywhere in the cube
 * the other half.
 *
 * So the colours carry the contrast themselves, and they are the other end of
 * 1's: where 1 puts vivid colours on black, this is dark and spent. Deep stops
 * are a colour barely above black; light ones an overcast grey with a tint,
 * never white; every fourth stop a plain grey. The hues are walked from a ring
 * of 256 as 1 and 2 walk theirs, so neighbouring hues recur the way theirs do,
 * over the hues of smog -- ochre and olive, 30 to 100 degrees, and steel, slate
 * and a dusty violet, 190 to 290 -- and, for a quarter of the ring, the reds
 * that go with it, crimson to scarlet: oxblood when deep and a dimmed scarlet
 * when light, both saturated, since a dark red that is not is read as brown.
 * One stop in five reads as red, against one in thirty in 1. Grey, a tinted
 * near-black, dark grey, burgundy, red, slate and dark brown are its commonest
 * colours by name, and against 1, 2 and 3 it stands 0.43 to 0.56 off. */
#define RING4_A 209
#define RING4_C 111
static int ring4(int x)
{
    return (RING4_A * x + RING4_C) & 255;
}

/* how many of the ring's 256 bytes go to crimson and scarlet */
#define SMOG_RED 64

static void smog_segments(int whitemode, int nsegments)
{
    int x = prand(256);
    int i;
    for (i = 0; i < nsegments; i++) {
        if (i % 4 == 3) {
            /* the plain grey of the air */
            int v = 77 + prand(64);
            int s = prand(8);
            hsv_to_rgb(0, s, v, colors[i], colors[i] + 1, colors[i] + 2);
            continue;
        }
        int u = x, h, s, v;
        int deep = (i & 1) == whitemode;
        x = ring4(ring4(ring4(x)));
        if (u < SMOG_RED) {
            /* 350 to 372 degrees, in hsv_to_rgb's 256ths of a turn */
            h = (350 + 22 * u / SMOG_RED) * 256 / 360;
            s = deep ? 191 + prand(65) : 178 + prand(65);
            v = deep ? 71 + prand(44) : 153 + prand(65);
        } else {
            int deg = (u - SMOG_RED) * 170 / (256 - SMOG_RED);
            deg = deg < 70 ? 30 + deg : 120 + deg;
            h = deg * 256 / 360;
            s = deep ? 89 + prand(77) : 15 + prand(62);
            v = deep ? 18 + prand(39) : 140 + prand(65);
        }
        hsv_to_rgb(h, s, v, colors[i], colors[i] + 1, colors[i] + 2);
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

/* 5: the two windows of 3, with other colours in them.
 *
 * 3 puts warm bright colours -- violet through red to orange -- over dark
 * greens and teals, and a bright colour on a dark green is what its commonest
 * pairs are. This turns it round: warm colours of every brightness, from rose
 * through red, burgundy, brown and orange to ochre and yellow, over night
 * colours kept dark, from teal through petrol and blue to indigo. Both windows
 * are 125 degrees wide, near the 135 of 3's, in the skeleton 3 has: a black or
 * a white stop every third, the two between taken from the two windows in 3's
 * order. They were a third as wide at first and read as two colours. Navy is
 * its commonest colour, and navy with olive, burgundy or brown -- brown against
 * blue was asked for by name -- its commonest pairs. */
static void earth_segments(int /*whitemode*/, int nsegments)
{
    int i;
    for (i = 0; i < nsegments; i++) {
        if (!(i % 3)) {
            int v = (i % 6) ? 255 : 0;
            colors[i][0] = colors[i][1] = colors[i][2] = v;
        } else if (((i) % 6 > 3) ^ ((i) % 3 == 1)) {
            /* night: 160 to 285 degrees, dark */
            int h = 114 + prand(89);
            int s = 128 + prand(128);
            int v = 51 + prand(103);
            hsv_to_rgb(h, s, v, colors[i], colors[i] + 1, colors[i] + 2);
        } else {
            /* warm: 320 round to 85 degrees, any brightness */
            int h = 228 + prand(89);
            int s = 89 + prand(154);
            int v = 89 + prand(167);
            hsv_to_rgb(h, s, v, colors[i], colors[i] + 1, colors[i] + 2);
        }
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

/* 6: the pairs the first three favour, side by side.
 *
 * Every pair of colours that one of the first three puts one after the other at
 * least three times as often as chance would, measured over twenty thousand
 * palettes each, with the two colours as they actually came out on average.
 * There they are one pair in a hundred or so and always an anchor apart; here
 * each pair is set down whole between a black stop and a white one, the two
 * colours touching, so what makes the pair is what the picture shows. Which
 * way round, and a nudge of up to ten in each channel so that one pair is not
 * always one colour, are left to the dice. */
static const unsigned char favoured_pairs[][2][3] = {
    {{205, 130, 147}, {243, 249, 113}}, /* pink + yellow, x23.8 in 2 */
    {{82, 35, 32}, {217, 158, 127}}, /* dark brown + salmon, x10.7 in 2 */
    {{87, 68, 45}, {152, 241, 214}}, /* dark brown + mint, x10.7 in 2 */
    {{228, 77, 2}, {18, 99, 96}}, /* red + petrol, x8.0 in 2 */
    {{205, 130, 147}, {27, 184, 145}}, /* pink + turquoise, x7.8 in 1 */
    {{91, 120, 81}, {227, 224, 153}}, /* grey-green + sand, x6.7 in 1 */
    {{145, 41, 68}, {160, 144, 35}}, /* burgundy + olive, x6.5 in 2 */
    {{12, 85, 106}, {153, 94, 63}}, /* petrol + brown, x6.4 in 2 */
    {{135, 180, 221}, {82, 35, 32}}, /* sky blue + dark brown, x6.4 in 2 */
    {{149, 170, 155}, {158, 127, 76}}, /* grey-green + taupe, x6.1 in 2 */
    {{203, 22, 188}, {210, 90, 86}}, /* magenta + salmon, x5.4 in 1 */
    {{216, 49, 22}, {151, 132, 109}}, /* red + taupe, x5.3 in 2 */
    {{129, 231, 42}, {159, 219, 163}}, /* lime + mint, x5.2 in 2 */
    {{126, 223, 44}, {152, 113, 86}}, /* lime + taupe, x5.2 in 1 */
    {{58, 200, 217}, {91, 128, 91}}, /* cyan + grey-green, x5.1 in 1 */
    {{31, 87, 29}, {169, 155, 56}}, /* dark green + olive, x5.0 in 1 */
    {{56, 17, 118}, {158, 127, 76}}, /* navy + taupe, x4.9 in 1 */
    {{155, 56, 17}, {217, 158, 127}}, /* brown + salmon, x4.9 in 1 */
    {{77, 2, 19}, {155, 56, 17}}, /* burgundy + brown, x4.7 in 1 */
    {{203, 34, 57}, {78, 62, 244}}, /* red + blue, x4.4 in 2 */
    {{127, 76, 149}, {17, 118, 119}}, /* plum + petrol, x4.4 in 1 */
    {{215, 68, 45}, {95, 44, 117}}, /* red + plum, x4.4 in 1 */
    {{5, 90, 139}, {19, 80, 73}}, /* slate + petrol, x4.3 in 1 */
    {{24, 113, 86}, {215, 196, 173}}, /* dark green + sand, x4.3 in 2 */
    {{139, 104, 129}, {38, 103, 20}}, /* mauve + dark green, x4.3 in 2 */
    {{202, 59, 88}, {31, 108, 53}}, /* burgundy + dark green, x4.3 in 2 */
    {{227, 224, 153}, {94, 63, 12}}, /* sand + brown, x4.3 in 2 */
    {{123, 24, 113}, {3, 128, 185}}, /* plum + slate, x4.3 in 1 */
    {{226, 220, 55}, {157, 224, 176}}, /* yellow + mint, x4.1 in 1 */
    {{163, 26, 71}, {29, 146, 213}}, /* burgundy + sky blue, x4.1 in 1 */
    {{115, 14, 120}, {78, 62, 39}}, /* plum + dark brown, x4.0 in 1 */
    {{206, 42, 159}, {141, 139, 53}}, /* fuchsia + olive, x4.0 in 1 */
    {{118, 249, 245}, {150, 206, 186}}, /* cyan + mint, x3.9 in 1 */
    {{35, 32, 217}, {158, 127, 76}}, /* blue + taupe, x3.9 in 2 */
    {{108, 126, 22}, {224, 208, 146}}, /* olive + sand, x3.9 in 2 */
    {{97, 24, 76}, {6, 16, 62}}, /* plum + navy, x3.9 in 1 */
    {{143, 28, 37}, {161, 198, 135}}, /* burgundy + mint, x3.9 in 1 */
    {{141, 66, 83}, {219, 120, 81}}, /* burgundy + salmon, x3.8 in 1 */
    {{98, 243, 176}, {87, 68, 45}}, /* turquoise + dark brown, x3.8 in 2 */
    {{11, 104, 129}, {93, 210, 163}}, /* petrol + turquoise, x3.8 in 2 */
    {{208, 67, 156}, {69, 118, 168}}, /* fuchsia + slate, x3.7 in 2 */
    {{85, 106, 91}, {94, 63, 12}}, /* grey-green + brown, x3.7 in 2 */
    {{140, 213, 234}, {146, 99, 96}}, /* cyan + taupe, x3.6 in 1 */
    {{26, 75, 40}, {212, 125, 114}}, /* dark green + salmon, x3.6 in 1 */
    {{166, 103, 148}, {46, 207, 220}}, /* mauve + cyan, x3.6 in 1 */
    {{10, 123, 152}, {241, 214, 87}}, /* petrol + yellow, x3.6 in 2 */
    {{39, 30, 73}, {218, 213, 52}}, /* navy + yellow, x3.6 in 2 */
    {{200, 97, 134}, {229, 186, 107}}, /* mauve + sand, x3.6 in 2 */
    {{193, 102, 167}, {154, 203, 168}}, /* mauve + mint, x3.6 in 2 */
    {{71, 116, 157}, {153, 94, 63}}, /* slate + brown, x3.6 in 1 */
    {{68, 45, 98}, {10, 123, 152}}, /* navy + petrol, x3.5 in 1 */
    {{127, 114, 196}, {30, 105, 51}}, /* lilac + dark green, x3.5 in 1 */
    {{229, 186, 107}, {243, 176, 41}}, /* sand + orange, x3.4 in 1 */
    {{93, 137, 88}, {113, 140, 11}}, /* grey-green + olive, x3.3 in 2 */
    {{75, 143, 167}, {72, 34, 27}}, /* slate + dark brown, x3.3 in 1 */
    {{60, 64, 221}, {35, 205, 227}}, /* blue + cyan, x3.3 in 1 */
    {{134, 71, 116}, {224, 153, 94}}, /* mauve + salmon, x3.3 in 1 */
    {{43, 65, 146}, {87, 222, 18}}, /* navy + lime, x3.3 in 1 */
    {{118, 119, 228}, {155, 56, 17}}, /* lilac + brown, x3.2 in 2 */
    {{122, 43, 136}, {33, 70, 7}}, /* plum + dark green, x3.2 in 2 */
    {{23, 111, 84}, {226, 116, 48}}, /* dark green + orange, x3.2 in 1 */
    {{143, 28, 37}, {250, 171, 8}}, /* burgundy + orange, x3.1 in 2 */
    {{250, 171, 8}, {161, 198, 135}}, /* orange + mint, x3.1 in 2 */
    {{245, 138, 251}, {215, 196, 173}}, /* lilac + sand, x3.0 in 1 */
    {{139, 104, 129}, {73, 78, 111}}, /* mauve + slate, x3.0 in 1 */
    {{131, 0, 57}, {245, 138, 251}}, /* burgundy + lilac, x3.0 in 1 */
};
#define NFAVOURED ((int)(sizeof(favoured_pairs) / sizeof(favoured_pairs[0])))

static int nudge(int c)
{
    c += prand(21) - 10;
    return c < 0 ? 0 : (c > 255 ? 255 : c);
}

static void pair_segments(int whitemode, int nsegments)
{
    int i, p = 0, flip = 0;
    for (i = 0; i < nsegments; i++) {
        int slot = i % 4;
        if (slot == 0 || slot == 3) {
            int v = (slot == 0) != !!whitemode ? 0 : 255;
            colors[i][0] = colors[i][1] = colors[i][2] = v;
            continue;
        }
        if (slot == 1) {
            p = prand(NFAVOURED);
            flip = prand(2);
        }
        const unsigned char *c = favoured_pairs[p][(slot == 2) ^ flip];
        colors[i][0] = nudge(c[0]);
        colors[i][1] = nudge(c[1]);
        colors[i][2] = nudge(c[2]);
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

/* 7: every stop a colour at random, and nothing else.
 *
 * It was 2 with the colours made truly random, black, a colour, white, and it
 * showed as 2 does: the black and the white every third stop are what the eye
 * reads, and against 2 it measured 0.00 -- the same palette as far as the
 * colours it shows go. So there is no skeleton at all: every channel of every
 * stop is drawn on its own from the top of the generator, any of sixteen million
 * colours, the light and the dark wherever the dice put them. It has less range
 * from dark to light than the anchored ones, which is what random means. */
static void random_segments(int /*whitemode*/, int nsegments)
{
    int i;
    for (i = 0; i < nsegments; i++) {
        colors[i][0] = prand(256);
        colors[i][1] = prand(256);
        colors[i][2] = prand(256);
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

/* What each is called in the palette dialog, counting from one as the dialog
 * does. */
const char *palette_algorithm_name(int algorithm)
{
    static const char *const names[PALGORITHMS] = {
        "Dark and colour", "Black, colour, white", "Warm over dark green",
        "Smog", "Warm over night", "Favoured pairs", "Random colours"};
    if (algorithm < 1 || algorithm > PALGORITHMS)
        return "";
    return names[algorithm - 1];
}

#define MYLONG_MAX 0xffffff
#define rrandom(i) ((int)(((int)XaoS_random() / (double)MYLONG_MAX) * (i)))
/*Do not use modulo type random since it should bring very different results
 *for slightly different sizes
 */

int mkpalette(struct palette *c, int seed, int algorithm)
{
    int i, ncolors = c->size;
    int whitemode;
    int i1;

    context = c;
    needupdate = 0;

    if (c->flags & DONOTCHANGE)
        return 0;
    XaoS_srandom(seed);
    seed = (int)XaoS_random();
    whitemode = (int)XaoS_random() % 2;

    if ((c->flags & UNKNOWNENTRIES) || !c->size) {
        maxentries = context->maxentries - nprecells;
        segmentsize = (rrandom(maxentries / 2)) & (~3);
        if (segmentsize < 1)
            segmentsize = 1;
    } else {
        if (maxentries > 8) {
            int qq = 255;

            maxentries = context->maxentries - nprecells;
            segmentsize = rrandom(qq / 3 + 4);
            segmentsize += rrandom(qq / 3 + 4);
            segmentsize += rrandom(qq / 3 + 4);
            segmentsize += rrandom(
                qq / 3 + 4); /*Make smaller segments with higher probability */

            segmentsize = abs(segmentsize / 2 - qq / 3 + 3);
            if (segmentsize < 8)
                segmentsize = 8;
            if (segmentsize > maxentries / 3)
                segmentsize = maxentries / 3;
        }
    }

    if (c->flags & UNKNOWNENTRIES)
        i = rrandom(maxentries);
    else
        i = (maxentries + segmentsize - 5) / segmentsize;

    if (i < 0)
        i = 1;
    if (i > MAXNSEGMENTS)
        i1 = MAXNSEGMENTS;
    else
        i1 = i;

    XaoS_srandom(seed);

    switch (algorithm) {
        case 6:
            random_segments(whitemode, i1);
            break;
        case 5:
            pair_segments(whitemode, i1);
            break;
        case 4:
            earth_segments(whitemode, i1);
            break;
        case 3:
            smog_segments(whitemode, i1);
            break;
        case 2:
            randomize_segments3(whitemode, i1);
            break;
        case 1:
            randomize_segments2(whitemode, i1);
            break;
        default:
            randomize_segments(whitemode, i1);
    }
    mksmooth(i1, i);

    if (!(c->flags & FINISHLATER)) {
        if (c->allocfinished != NULL)
            c->allocfinished(c);
    } else
        c->flags |= UNFINISHED;
    if (context->size != ncolors || needupdate) {
        context->version++;
        return 1;
    }

    return 0;
}

int mkstereogrampalette(struct palette *c)
{
    int i, ncolors = c->size;
    context = c;
    needupdate = 0;
    for (i = 0; i < 16; i++)
        allocate(i * 4, i * 4, i * 16, i == 0);
    if (!(c->flags & FINISHLATER)) {
        if (c->allocfinished != NULL)
            c->allocfinished(c);
    } else
        c->flags |= UNFINISHED;
    if (context->size != ncolors || needupdate) {
        context->version++;
        return 1;
    }
    return 0;
}

int mkstarfieldpalette(struct palette *c)
{
    int i, ncolors = c->size;
    context = c;
    needupdate = 0;
    for (i = 0; i < 16; i++)
        if (i % 2)
            allocate(i * 4, i * 4, i * 16, i == 0);
        else
            allocate(i * 16, i * 16, i * 16, i == 0);
    if (!(c->flags & FINISHLATER)) {
        if (c->allocfinished != NULL)
            c->allocfinished(c);
    } else
        c->flags |= UNFINISHED;
    if (context->size != ncolors || needupdate) {
        context->version++;
        return 1;
    }
    return 0;
}

int mkblurpalette(struct palette *c)
{
    int i, ncolors = c->size;
    context = c;
    needupdate = 0;
    for (i = 0; i < 63; i++)
        allocate(i * 2, i * 2, i * 4, i == 0);
    allocate(i * 2, i * 2, i * 4 - 1, 0);
    if (!(c->flags & FINISHLATER)) {
        if (c->allocfinished != NULL)
            c->allocfinished(c);
    } else
        c->flags |= UNFINISHED;
    if (context->size != ncolors || needupdate) {
        context->version++;
        return 1;
    }
    return 0;
}

int mkgraypalette(struct palette *c)
{
    int i, ncolors = c->size;
    context = c;
    needupdate = 0;
    for (i = 0; i < 64; i++)
        allocate(i * 4, i * 4, i * 4, i == 0);
    for (i = 0; i < 16; i++)
        allocate(255, 255 - i * 16, 255 - i * 16, 0);
    for (i = 0; i < 16; i++)
        allocate(255 - i * 16, 0, 0, 0);
    if (!(c->flags & FINISHLATER)) {
        if (c->allocfinished != NULL)
            c->allocfinished(c);
    } else
        c->flags |= UNFINISHED;
    if (context->size != ncolors || needupdate) {
        context->version++;
        return 1;
    }
    return 0;
}

int mkdefaultpalette(struct palette *c)
{
    int i, ncolors = c->size;
    context = c;
    needupdate = 0;
    segmentsize = 8;

    if (c->flags & DONOTCHANGE)
        return 0;
    memcpy(colors, colors1, sizeof(colors1));
    maxentries = context->maxentries - nprecells;
    if (c->flags & UNKNOWNENTRIES)
        i = 128 / 8;
    else
        i = (maxentries + 3) / 8;
    if (i < 0)
        i = 1;
    mksmooth(255 / 8, i);
    if (!(c->flags & FINISHLATER)) {
        if (c->allocfinished != NULL)
            c->allocfinished(c);
    } else
        c->flags |= UNFINISHED;
    if (context->size != ncolors || needupdate) {
        context->version++;
        return 1;
    }
    return 0;
}

int shiftpalette(struct palette *c, int shift)
{
    if (!c->size)
        return 0;

    while (shift < 0)
        shift += c->size - 1;
    shift = shift % (c->size - 1);

    if (!shift)
        return 0;

    if (c->cyclecolors != NULL) {
        if (c->flags & UNFINISHED) {
            cycle_entries(c, shift);
        } else {
            c->cyclecolors(c, shift);
        }
        return 0;
    }

    if (c->type & (TRUECOLOR | TRUECOLOR24 | TRUECOLOR16)) {
        int i;
        int i3;
        int *co;

        if (shift > 0)
            shift %= c->size - 1;
        else
            shift = -((-shift) % (c->size - 1));
        if (!shift)
            return 0;
        co = (int *)malloc(c->size * sizeof(*co));
        memcpy(co, c->pixels, sizeof(*co) * c->size);
        i3 = (c->size - 1 + shift) % (c->size - 1) + 1;
        for (i = 1; i < c->size; i++) {
            c->pixels[i] = co[i3];
            i3++;
            if (i3 >= c->size)
                i3 = 1;
        }
        c->version++;
        free(co);
    }
    return 1;
}

static int allocrgb(struct palette */*c*/, int r1, int g1, int b1)
{
    int r, g, b;
    int f = 1;
    for (g = 0; g < g1; g++)
        for (b = 0; b < b1; b++) {
            for (r = 0; r < r1; r++) {
                if (!allocate(r * 255 / (r1 - 1), g * 255 / (g1 - 1),
                              b * 255 / (b1 - 1), f))
                    return 0;
                f = 0;
            }
        }
    return 1;
}

int mkrgb(struct palette *c)
{
    int ncolors = c->size;
    int red = 8, green = 8, blue = 4;

    context = c;
    needupdate = 0;

    if (c->flags & UNKNOWNENTRIES) {
        while (blue > 0) {
            if (allocrgb(c, red, green, blue))
                break;
            red--;
            if (allocrgb(c, red, green, blue))
                break;
            green--;
            if (allocrgb(c, red, green, blue))
                break;
            red--;
            if (allocrgb(c, red, green, blue))
                break;
            green--;
            if (allocrgb(c, red, green, blue))
                break;
            blue--;
        }
    } else {
        number_t n =
            pow((c->maxentries - nprecells) / (0.5 * 0.2 * 0.3), 1.0 / 3);
        green = (int)(n * 0.5);
        blue = (int)(n * 0.2);
        red = (int)(n * 0.3);
        while ((blue + 1) * red * green < (c->maxentries - nprecells))
            blue++;
        while ((red + 1) * blue * green < (c->maxentries - nprecells))
            red++;
        while ((green + 1) * blue * red < (c->maxentries - nprecells))
            green++;
        allocrgb(c, red, green, blue);
    }

    if (!(c->flags & FINISHLATER)) {
        if (c->allocfinished != NULL)
            c->allocfinished(c);
    } else
        c->flags |= UNFINISHED;
    if (context->size != ncolors || needupdate) {
        context->version++;
    }

    return red * 256 * 256 + green * 256 + blue;
}

void getPaletteColor(struct palette *src, int seed, int algorithm, int shift, int newColors[101][3]) {
    mkpalette(src, seed, algorithm);
    shiftpalette(src, shift);
    for (int i = 0; i < 101; i++) {
        int r = 0, g = 0, b = 0;
        switch (src->type) {
            case SMALLITER:
                r = g = b = i;
                break;
            case LARGEITER:
                r = g = b = i / 256;
                break;
            case GRAYSCALE:
                r = g = b = src->pixels[i];
                break;
            case C256:
            case FIXEDCOLOR:
            case MBITMAP:
            case LBITMAP:
            case MIBITMAP:
            case LIBITMAP:
                r = src->rgb[src->pixels[i]][0];
                g = src->rgb[src->pixels[i]][1];
                b = src->rgb[src->pixels[i]][2];
                break;
            case TRUECOLOR:
            case TRUECOLOR16:
            case TRUECOLOR24:
                r = (((src->pixels[i] & src->info.truec.rmask) >>
                      src->info.truec.rshift))
                    << src->info.truec.rprec;
                g = (((src->pixels[i] & src->info.truec.gmask) >>
                      src->info.truec.gshift))
                    << src->info.truec.gprec;
                b = (((src->pixels[i] & src->info.truec.bmask) >>
                      src->info.truec.bshift))
                    << src->info.truec.bprec;
                break;
        }
        newColors[i][0]=r;
        newColors[i][1]=g;
        newColors[i][2]=b;
    }
}

void getDEFSEGMENTColor(unsigned char newColors[][3]) {
    memcpy(newColors, colors, sizeof(colors1));
}

int mkcustompalette(struct palette *c, unsigned char newColors[31][3])
{
    int i, ncolors = c->size;
    context = c;
    needupdate = 0;
    segmentsize = 8;

    if (c->flags & DONOTCHANGE)
        return 0;
    memcpy(colors, newColors, sizeof(colors1));
    maxentries = context->maxentries - nprecells;
    if (c->flags & UNKNOWNENTRIES)
        i = 128 / 8;
    else
        i = (maxentries + 3) / 8;
    if (i < 0)
        i = 1;
    mksmooth(255 / 8, i);
    if (!(c->flags & FINISHLATER)) {
        if (c->allocfinished != NULL)
            c->allocfinished(c);
    } else
        c->flags |= UNFINISHED;
    if (context->size != ncolors || needupdate) {
        context->version++;
        return 1;
    }
    return 0;
}
