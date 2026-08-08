/* The overlay save/restore protocol in wstack.cpp.
 *
 * An overlay -- the status bar, a message, the menu -- is drawn into the
 * fractal image itself, so the pixels it covers have to be kept somewhere and
 * put back when it goes away. Where "somewhere" is turns out to matter a great
 * deal: the code used to park them in image->oldlines, which is not scratch
 * space but the previous frame, the one moveoldpoints() in zoom.cpp reads to
 * reuse rows while zooming. Overlay pixels landing there come back out as
 * fractal data and get stretched across the image.
 *
 * These checks pin the properties that failure violated, so that a future
 * change back to a shared scratch area is caught here rather than as smeared
 * text on someone's screen.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "filter.h"
#include "ui_helper.h"

#define W 64
#define H 32
#define OVX 8
#define OVY 4
#define OVW 40
#define OVH 6
#define OVERLAY_PIXEL 0xFEEDFACE

static int failures = 0;

static void check(int ok, const char *what)
{
    printf("%-6s %s\n", ok ? "ok" : "FAIL", what);
    if (!ok)
        failures++;
}

/* A recognisable value per pixel, so that a restore putting back the wrong
 * row or the wrong column is visible rather than merely plausible. */
static unsigned int pattern(int x, int y) { return 0x01000000u + y * W + x; }

static void fill_pattern(struct image *img)
{
    for (int y = 0; y < H; y++) {
        unsigned int *row = (unsigned int *)img->currlines[y];
        for (int x = 0; x < W; x++)
            row[x] = pattern(x, y);
    }
}

static int differs_from_pattern(struct image *img)
{
    for (int y = 0; y < H; y++) {
        unsigned int *row = (unsigned int *)img->currlines[y];
        for (int x = 0; x < W; x++)
            if (row[x] != pattern(x, y))
                return 1;
    }
    return 0;
}

static void overlay_pos(struct uih_context *, int *x, int *y, int *w, int *h,
                        void *)
{
    *x = OVX;
    *y = OVY;
    *w = OVW;
    *h = OVH;
}

static void overlay_draw(struct uih_context *c, void *)
{
    for (int y = OVY; y < OVY + OVH; y++) {
        unsigned int *row = (unsigned int *)c->image->currlines[y];
        for (int x = OVX; x < OVX + OVW; x++)
            row[x] = OVERLAY_PIXEL;
    }
}

static int overlay_pixels_in(pixel_t *buf, int count)
{
    unsigned int *p = (unsigned int *)buf;
    int n = 0;
    for (int i = 0; i < count; i++)
        if (p[i] == OVERLAY_PIXEL)
            n++;
    return n;
}

int main(void)
{
    const int pixels = W * H;
    pixel_t *buf1 = (pixel_t *)calloc(pixels, sizeof(unsigned int));
    pixel_t *buf2 = (pixel_t *)calloc(pixels, sizeof(unsigned int));
    union paletteinfo info;
    info.truec.rmask = 0xff0000;
    info.truec.gmask = 0x00ff00;
    info.truec.bmask = 0x0000ff;
    struct palette *pal =
        createpalette(0, 0, TRUECOLOR, 0, 0, NULL, NULL, NULL, NULL, &info);
    struct image *img =
        create_image_cont(W, H, W * sizeof(unsigned int), 2, buf1, buf2, pal,
                          flipgeneric, 0, 1, 1);
    if (!buf1 || !buf2 || !pal || !img) {
        printf("FAIL  cannot set up the image\n");
        return 1;
    }

    struct uih_context ctx;
    memset(&ctx, 0, sizeof ctx);
    ctx.image = img;

    struct uih_window *w =
        uih_registerw(&ctx, overlay_pos, overlay_draw, NULL, 0);
    if (!w) {
        printf("FAIL  cannot register the overlay\n");
        return 1;
    }

    /* Two frames, because the buffer holding the saved pixels is now kept
     * between them: a second pass catches a stale restore that a single one
     * would not. */
    for (int frame = 0; frame < 2; frame++) {
        fill_pattern(img);
        memset(img->oldlines[0], 0, 0); /* oldlines rows are not contiguous */
        for (int y = 0; y < H; y++)
            memset(img->oldlines[y], 0xA5, W * sizeof(unsigned int));

        uih_drawwindows(&ctx);
        check(overlay_pixels_in((pixel_t *)((unsigned int *)img->currlines[OVY] +
                                            OVX),
                                OVW) == OVW,
              frame ? "frame 2: the overlay reaches the image"
                    : "frame 1: the overlay reaches the image");

        /* The point of the whole exercise: whatever the overlay needed to
         * save, none of it may go through the other frame buffer. */
        int leaked = 0;
        for (int y = 0; y < H; y++) {
            unsigned char *row = (unsigned char *)img->oldlines[y];
            for (unsigned int i = 0; i < W * sizeof(unsigned int); i++)
                if (row[i] != 0xA5)
                    leaked = 1;
        }
        check(!leaked, frame ? "frame 2: the previous frame is left alone"
                             : "frame 1: the previous frame is left alone");

        uih_clearwindows(&ctx);
        check(!differs_from_pattern(img),
              frame ? "frame 2: the image is restored exactly"
                    : "frame 1: the image is restored exactly");
    }

    /* A clear with nothing drawn must not put anything back: drawwindows
     * skips a window whose area is already covered by another, and a restore
     * that ignored that would paint a stale frame over a fresh one. */
    fill_pattern(img);
    uih_clearwindows(&ctx);
    check(!differs_from_pattern(img), "a clear with nothing drawn is a no-op");

    uih_removew(&ctx, w);

    if (failures)
        printf("\n%d check(s) failed\n", failures);
    return failures != 0;
}
