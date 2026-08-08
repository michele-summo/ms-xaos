#include <cstdlib>
#include <cstring>
#include <cassert>

#include "config.h"
#include "filter.h"
#include "ui_helper.h"
#include "grlib.h"

/* The "windows" implemented here are just regions of the image, not user
 * interface windows.  They are used to save the pixels that are going to
 * be drawn over (e.g., with a text message) so they can be restored when the
 * overlay should disappear. Windows are redrawn every frame and removed before
 * calculation.
 */

struct uih_window *uih_registerw(struct uih_context *uih, uih_getposfunc getpos,
                                 uih_drawfunc draw, void *data, int flags)
{
    struct uih_window *w =
        (struct uih_window *)calloc(1, sizeof(struct uih_window));
    struct uih_window *w1;
    assert(uih != NULL && getpos != NULL && draw != NULL && flags >= 0);
    if (w == NULL)
        return NULL;
    uih_clearwindows(uih);
    w1 = uih->wtop;
    w->getpos = getpos;
    w->draw = draw;
    w->data = data;
    w->flags = flags;
    w->savedline = -1;
    w->saveddata = NULL;
    w->next = NULL;
    if (w1 == NULL) {
        uih->wtop = w;
    } else {
        while (w1->next != NULL)
            w1 = w1->next;
        w1->next = w;
    }
    w->previous = w1;
    w->x = -65536;
    return w;
}

void uih_setline(struct uih_context *uih, struct uih_window *w, int color,
                 int x1, int y1, int x2, int y2)
{
    if (w->savedline != color || w->x != x1 || w->y != y1 ||
        w->width != x2 - x1 || w->height != y2 - y1) {
        uih_clearwindows(uih);
        uih->display = 1;
        w->savedline = color;
        w->x = x1;
        w->y = y1;
        w->width = x2 - x1;
        w->height = y2 - y1;
    }
}

struct uih_window *uih_registerline(struct uih_context *uih, int color, int x1,
                                    int y1, int x2, int y2)
{
    struct uih_window *w =
        (struct uih_window *)calloc(1, sizeof(struct uih_window));
    struct uih_window *w1;
    if (w == NULL)
        return NULL;
    uih_clearwindows(uih);
    w1 = uih->wtop;
    uih->display = 1;
    w->getpos = NULL;
    w->savedline = color;
    w->flags = 0;
    w->x = x1;
    w->y = y1;
    w->width = x2 - x1;
    w->height = y2 - y1;
    w->saveddata = NULL;
    w->next = NULL;
    if (w1 == NULL) {
        uih->wtop = w;
    } else {
        while (w1->next != NULL)
            w1 = w1->next;
        w1->next = w;
    }
    w->previous = w1;
    return w;
}

void uih_removew(struct uih_context *uih, struct uih_window *w)
{
    uih_clearwindows(uih);
    assert(uih->wtop != NULL);
    assert(w != NULL);
    uih->display = 1;

    if (w->previous == NULL) {
        assert(uih->wtop == w);
        uih->wtop = w->next;
    } else {
        w->previous->next = w->next;
    }
    if (w->next != NULL) {
        w->next->previous = w->previous;
    }
    /* The saved area outlives the frame it was taken in, so releasing the
     * window is what releases it. */
    free(w->saveddata);
    free(w);
}

/*Remove all drawn windows from screen */
void uih_clearwindows(struct uih_context *uih)
{
    struct uih_window *w = uih->wtop;
    if (!uih->wdisplayed)
        return;
    uih->wdisplayed = 0;
    if (uih->wflipped)
        uih->image->flip(uih->image), uih->wflipped = 0;
    while (w) {
        if (w->getpos == NULL) {
            if (w->saveddata != NULL) {
                xrestoreline(uih->image, w->saveddata, w->x, w->y,
                             w->width + w->x, w->height + w->y);
                free(w->saveddata);
                w->saveddata = NULL;
            }
        } else {
            /* savedvalid is the marker for "this window was saved during
             * the current frame": drawwindows skips a window whose area
             * another one already covers, and without the marker the stale
             * contents of its buffer would be painted back over the image. */
            if (w->savedvalid) {
                int i;
                int xskip = w->x * uih->image->bytesperpixel;
                int width = w->width * uih->image->bytesperpixel;
                if (!uih->image->bytesperpixel) {
                    xskip = w->x / 8;
                    width = (w->x + w->width + 7) / 8 - xskip;
                }
                assert(w->width);
                assert(w->height);
                assert(w->x >= 0);
                assert(w->y >= 0);
                assert(w->x + w->width <= uih->image->width);
                assert(w->y + w->height <= uih->image->height);
                assert(w->saveddata);
                for (i = w->y; i < w->y + w->height; i++) {
                    unsigned char *data = uih->image->currlines[i] + xskip;
                    memcpy(data, w->saveddata + (i - w->y) * width, width);
                }
                /* The buffer stays allocated for the next frame; only the
                 * marker is cleared. */
                w->savedvalid = 0;
            }
        }
        w = w->next;
    }
}

void uih_drawwindows(struct uih_context *uih)
{
    struct uih_window *w = uih->wtop;
    struct image *img = uih->image;
    if (uih->wdisplayed)
        return;
    uih->wdisplayed = 1;
    while (w) {
        if (w->getpos != NULL) {
            w->getpos(uih, &w->x, &w->y, &w->width, &w->height, w->data);
            if (w->x < 0)
                w->width -= w->x, w->x = 0;
            if (w->y < 0)
                w->height -= w->y, w->y = 0;
            if (w->x > img->width)
                w->width = 0, w->height = 0, w->x = 0;
            if (w->y > img->height)
                w->width = 0, w->height = 0, w->y = 0;
            if (w->x + w->width > img->width)
                w->width = img->width - w->x;
            if (w->y + w->height > img->height)
                w->height = img->height - w->y;
            if (w->width < 0)
                w->width = 0;
            if (w->height < 0)
                w->height = 0;
            assert(w->width >= 0);
            assert(w->height >= 0);
            assert(w->x >= 0);
            assert(w->y >= 0);
            assert(w->x + w->width <= uih->image->width);
            assert(w->y + w->height <= uih->image->height);
        }
        w = w->next;
    }
    /* Windows covering more than half the image used to be handled by copying
     * the whole image into oldlines and flipping the two buffers, so that the
     * restore was a pointer swap. That has the same flaw as parking the small
     * saves there: oldlines is the previous frame, which the zoom engine
     * still needs. Every window now saves what it covers into its own buffer,
     * which costs one more copy back on the largest windows and nothing at
     * all on a status bar. */
    {
        int savedminx = -1;
        int savedmaxx = -1;
        int savedminy = -1;
        int savedmaxy = -1;
        uih->wflipped = 0;
        w = uih->wtop;
        while (w) {
            int i;
            if (w->getpos == NULL) {
                if ((w->x < savedminx || w->y < savedminy ||
                     w->x + w->width > savedmaxx ||
                     w->x + w->height > savedmaxy ||
                     w->x + w->width < savedminx ||
                     w->y + w->height < savedminy || w->x > savedmaxx ||
                     w->y > savedmaxy)) {
                    w->saveddata = xsaveline(uih->image, w->x, w->y,
                                             w->width + w->x, w->height + w->y);
                }
            } else {
                if (w->width && w->height &&
                    (w->x < savedminx || w->y < savedminy ||
                     w->x + w->width > savedmaxx ||
                     w->y + w->height > savedmaxy)) {
                    int xskip = w->x * uih->image->bytesperpixel;
                    int width = w->width * uih->image->bytesperpixel;
                    savedminx = w->x;
                    savedminy = w->y;
                    savedmaxx = w->x + w->width;
                    savedmaxy = w->y + w->height;
                    if (!uih->image->bytesperpixel) {
                        xskip = w->x / 8;
                        width = (w->x + w->width + 7) / 8 - xskip;
                    }
                    /* The area under the window goes into a buffer of its
                     * own. It used to be parked in image->oldlines unless
                     * PROTECTBUFFERS said otherwise, on the assumption that
                     * the second buffer is free scratch -- and it is not: it
                     * holds the previous frame, which is what moveoldpoints()
                     * in zoom.cpp reads to reuse rows while zooming. Since
                     * the progress callback draws the windows in the middle
                     * of a calculation, the two ran over each other: the
                     * engine picked the saved pixels up as if they were
                     * fractal data and stretched them across the image, and
                     * the restore afterwards read back whatever the engine
                     * had written over the save. That is the smeared status
                     * text, and the block of noise left behind it. */
                    int needed = width * w->height + 1;
                    if (w->savedsize < needed) {
                        free(w->saveddata);
                        w->saveddata = (char *)malloc(needed);
                        w->savedsize = w->saveddata ? needed : 0;
                    }
                    if (w->saveddata != NULL) {
                        for (i = w->y; i < w->y + w->height; i++) {
                            unsigned char *data = img->currlines[i] + xskip;
                            memcpy(w->saveddata + (i - w->y) * width, data,
                                   width);
                        }
                        w->savedvalid = 1;
                    }
                }
            }
            w = w->next;
        }
    }
    w = uih->wtop;
    while (w) {
        if (w->getpos == NULL) {
#define lwi 0
#define lwj 0
            xline(uih->image, w->x + lwi, w->y + lwj, w->width + w->x + lwi,
                  w->height + w->y + lwj, w->savedline);
        } else if (w->width && w->height) {
            w->draw(uih, w->data);
        }
        w = w->next;
    }
}
