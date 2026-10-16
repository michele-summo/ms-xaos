#ifndef _TEXT_H
#define _TEXT_H
#define TEXT_PRESSED 1
#include "ui_helper.h"
int xprint(struct image *image, void *font, int x, int y,
           const char *text, int fgcolor, int bgcolor, int mode);
int xtextwidth(struct image *image, void *font, const char *text);
int xtextheight(struct image *image, void *font);
int xtextcharw(struct image *image, void *font, const char c);
char *xsaveline(struct image *img, int x1, int y1, int x2, int y2);
void xrestoreline(struct image *img, char *data, int x1, int y1, int x2,
                  int y2);
void xline(struct image *img, int x1, int y1, int x2, int y2, int color);
void xprepareimage(struct image *img);
void overlayGrid(uih_context *c, int fgcolor);

#endif
