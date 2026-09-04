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


/* --- palettes that are a colour scheme rather than a scatter ---------------
 *
 * One rule runs through all four and is worth stating before them: every
 * segment keeps its colour. The three above get their character by pinning
 * segments to black and white, and a segment pinned to either has no hue at
 * all; do the same in a palette whose whole point is which hues it holds and
 * the hues stop being visible. Two of these four went that way -- alternating
 * near-black with a washed-out near-white -- and came out grey with a tint,
 * which is what "monochrome" means when it is a complaint. So contrast here is
 * between a deep colour and a bright one, never between nothing and nothing.
 *
 * The three above all work the same way: every second or third segment is
 * pinned to black or white and the rest are picked at random, which is what
 * gives XaoS its banded, high-contrast look. It is a good look and it is the
 * only one it had. These four choose their colours in relation to each other
 * instead, so that a palette holds together as one thing.
 *
 * They cost what the others cost. A palette is made once, when it is asked
 * for, and is only some four to nine segments long -- the two hundred odd
 * entries it ends up with are interpolated between those -- so nothing here is
 * per pixel and none of it does more per segment than one conversion out of
 * hue, saturation and value. That shortness is also why none of them tries to
 * lay a pattern along the palette: there is no length to lay one along.
 *
 * Each is driven by XaoS_random from the seed the caller set, so the same
 * algorithm and seed give the same palette, which is what lets a saved
 * position bring its colours back with it.
 */

/* A cycle right round the hue circle, at full strength: the smooth spectrum a
 * fractal is usually shown in elsewhere and could not be shown in here. One
 * turn, or two or three, so that a palette of the same length can band the
 * iterations coarsely or finely. */
static void spectrum_segments(int whitemode, int nsegments)
{
    int start = (int)XaoS_random() % 256;
    /* Round it one way or the other: warm to cool reads differently from cool
     * to warm, and with four to nine stops a single turn is all a cycle can
     * be -- two turns across five stops is not a spectrum, it is two hues
     * alternating, which is what algorithm 7 is for. */
    int direction = ((int)XaoS_random() & 1) ? 1 : -1;
    int bright = whitemode ? 255 : 220 + (int)XaoS_random() % 36;
    /* deep rather than dark: at a fifth of full value a hue is already too
     * near black to read, and half the circle went missing from the picture */
    int dark = 90 + (int)XaoS_random() % 35;
    /* Saturation runs the other way from brightness, which is how a gradient
     * gets its range without either end going colourless: deep and rich at one
     * end, light and tinted at the other. Holding it steady instead, the
     * channel sum could not span more than twice the value and the palette was
     * a fifth of what it can be. */
    int deepsat = 215 + (int)XaoS_random() % 41;
    int palesat = 70 + (int)XaoS_random() % 45;
    int span = nsegments > 1 ? nsegments - 1 : 1;
    int i;

    for (i = 0; i < nsegments; i++) {
        int h = start + direction * 256 * i / span;
        /* Round the circle in hue, and light and deep by turns in brightness.
         *
         * Swelling once across the whole palette instead put every band a
         * quarter of the way from its neighbour, and the colour then changes
         * over sixty entries where the three older ways change over thirty:
         * the gradient reads as slow. Alternating gives every band an edge
         * against both its neighbours, which is what makes a fractal's rings
         * visible, and costs the spectrum nothing -- the hue goes round all
         * the same. */
        int deep = !(i & 1);
        int val = deep ? dark : bright;
        int sat = deep ? deepsat : palesat;
        /* the first and last segment are one segment, so an odd count would
         * put two deep ones side by side at the seam */
        if ((nsegments & 1) && i == nsegments - 2) {
            val = (dark + bright) / 2;
            sat = (deepsat + palesat) / 2;
        }
        hsv_to_rgb(h, sat, val, colors[i], colors[i] + 1, colors[i] + 2);
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

/* Two hues far enough apart to be told apart, one owning the shadows and the
 * other the highlights, with the palette running from the first through the
 * second: a duotone, as a press does it with two inks.
 *
 * It was one hue with a few degrees of drift, and one hue is one colour: a
 * quarter of the ways to make a palette came out monochrome, which is a thing
 * to have one of and not two. Seventy degrees or more between the ends is what
 * makes the two ends read as two colours. */
static void duotone_segments(int whitemode, int nsegments)
{
    int shadow = (int)XaoS_random() % 256;
    /* between a fifth and two fifths of the way round, either way about */
    int apart = 50 + (int)XaoS_random() % 55;
    if ((int)XaoS_random() & 1)
        apart = -apart;
    int span = nsegments > 1 ? nsegments - 1 : 1;
    int i;

    for (i = 0; i < nsegments; i++) {
        /* The hue walks from the first ink to the second and back, so the
         * palette is still the two of them mixing; the weight alternates band
         * by band, so the picture has an edge every band instead of one long
         * fade from end to end. */
        int rise = i * 2 <= span ? i * 2 : (span - i) * 2;
        int top = (span / 2) * 2;
        int t = rise * 255 / (top ? top : 1);
        int deep = !(i & 1);
        int val = deep ? (whitemode ? 75 : 60) + t * 30 / 255
                       : 215 + t * 40 / 255;
        int sat = deep ? 210 + t * 45 / 255 : 95 + t * 45 / 255;
        if ((nsegments & 1) && i == nsegments - 2) {
            val = (val + 255) / 2;
            sat = sat * 2 / 3;
        }
        hsv_to_rgb(shadow + apart * t / 255, sat, val, colors[i],
                   colors[i] + 1, colors[i] + 2);
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

/* Three hues spread round the circle, taken in turn, dark and light
 * alternating: a scheme rather than a scatter, and one that has three colours
 * in it however few segments there are to hold them.
 *
 * This was hues from one narrow arc, which is a thing a painter would call
 * analogous and anybody else would call one colour. A narrow arc and four
 * segments is four shades of the same thing; a wide one is three colours that
 * still agree, which is what was wanted. */
static void triad_segments(int whitemode, int nsegments)
{
    int first = (int)XaoS_random() % 256;
    /* a third of the circle apart, give or take, so the three are spread but
     * not evenly enough to look mechanical */
    int step = 70 + (int)XaoS_random() % 32;
    int i;

    for (i = 0; i < nsegments; i++) {
        int h = first + step * (i % 3) + (int)XaoS_random() % 13 - 6;
        /* the bright turn is also the tinted one, for the same reason the
         * spectrum's is: two ends that differ only in value cannot span much */
        int s = (i & 1) ? 80 + (int)XaoS_random() % 55
                        : 200 + (int)XaoS_random() % 56;
        /* Deep and bright take turns by position, not by the dice: left to the
         * dice a short palette comes out all one weight. Deep, not black, and
         * bright, not white -- pinning the turns to those left thirty entries
         * of two hundred and fifty with any hue in them at all. */
        int v = (i & 1) ? 215 + (int)XaoS_random() % 41
                        : 85 + (int)XaoS_random() % 45;
        if (whitemode && (i & 1))
            s = s * 2 / 3;
        hsv_to_rgb(h, s, v, colors[i], colors[i] + 1, colors[i] + 2);
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
}

/* Two hues from opposite sides of the circle, taking turns: the loudest
 * agreement two colours can make, and the one that shows a boundary between
 * two bands most plainly. */
static void complementary_segments(int whitemode, int nsegments)
{
    int hue = (int)XaoS_random() % 256;
    int spread = (int)XaoS_random() % 24 - 12; /* not exactly opposite */
    int i;

    for (i = 0; i < nsegments; i++) {
        int h = (i & 1) ? hue + 128 + spread : hue;
        int s = 170 + (int)XaoS_random() % 86;
        int v = 110 + (int)XaoS_random() % 146;
        /* one segment in four goes to an extreme, which gives the pair
         * somewhere to come from and somewhere to go */
        if (i % 4 == 2) {
            s = 60 + (int)XaoS_random() % 60;
            v = whitemode ? 225 + (int)XaoS_random() % 31
                          : 45 + (int)XaoS_random() % 35;
        }
        hsv_to_rgb(h, s, v, colors[i], colors[i] + 1, colors[i] + 2);
    }
    colors[i - 1][0] = colors[0][0];
    colors[i - 1][1] = colors[0][1];
    colors[i - 1][2] = colors[0][2];
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
            complementary_segments(whitemode, i1);
            break;
        case 5:
            triad_segments(whitemode, i1);
            break;
        case 4:
            duotone_segments(whitemode, i1);
            break;
        case 3:
            spectrum_segments(whitemode, i1);
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
