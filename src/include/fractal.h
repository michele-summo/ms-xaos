/*
 *     XaoS, a fast portable realtime fractal zoomer
 *                  Copyright (C) 1996,1997 by
 *
 *      Jan Hubicka          (hubicka@paru.cas.cz)
 *      Thomas Marsh         (tmarsh@austin.ibm.com)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef FRACTAL1_H
#define FRACTAL1_H

#include "config.h"
#ifdef USE_SFFE
#include "sffe.h"
#endif

#define INCOLORING 11
#define OUTCOLORING (OutColormodeClass::ColOut_MAXMode + 1)
#define TCOLOR 15
#define COLORFUN 10

class OutColormodeClass {
public:
enum OutColormode {
                    ColOut_iter,
                    ColOut_iter_plus_real,
                    ColOut_iter_plus_imag,
                    ColOut_iter_plus_real_div_imag,
                    ColOut_iter_plus_real_plus_imag_plus_real_div_imag,
                    ColOut_binary_decomposition,
                    ColOut_biomorphs,
                    ColOut_potential,
                    ColOut_color_decomposition,
                    ColOut_smooth,
                    ColOut_smooth_log,
                    ColOut_True_color,
                    ColOut_MAXMode = ColOut_True_color
                  };
   enum OutColormode ColorMode;
   int AsInt() {return ColorMode; }
   OutColormodeClass &operator=(enum OutColormode Color) {ColorMode = Color; return *this;}

   bool operator==(enum OutColormode Color) {return ColorMode == Color;}
   bool operator!=(enum OutColormode Color) {return ColorMode != Color;}
   bool operator!=(const OutColormodeClass &Color) {return ColorMode == Color.ColorMode;}
};

typedef OutColormodeClass OutColormodeType;

typedef struct {
    number_t y0, k;
} symmetrytype;

struct symmetryinfo {
    number_t xsym, ysym;
    int nsymmetries;
    const symmetrytype *symmetry;
};
typedef struct {
    number_t mc, nc;
    number_t mi, ni;
} vrect;
typedef struct {
    number_t cr, ci;
    number_t rr, ri;
} vinfo;
typedef unsigned int (*iterationfunc)(number_t, number_t, number_t, number_t);
struct formula {
    int magic;
    iterationfunc calculate, calculate_periodicity, smooth_calculate,
        smooth_calculate_periodicity;
    void (*calculate_julia)(struct image *img, number_t pre, number_t pim);
    const char *name[2];
    const char *shortname;
    vinfo v;
    int hasperiodicity;
    int mandelbrot;
    number_t pre, pim;
    struct symmetryinfo out[OUTCOLORING + 1];
    struct symmetryinfo in[INCOLORING + 1];
    int flags;
};



/* The region the iteration has to leave for a point to count as escaped. It
 * has always been a circle; the shape is visible in the bands outside the set,
 * so the others are as much a drawing tool as a numerical one. */
#define BAILOUT_CIRCLE 0  /* |re|^2 + |im|^2  -- what every version did */
#define BAILOUT_SQUARE 1  /* either component alone */
#define BAILOUT_DIAMOND 2 /* |re| + |im| */
#define BAILOUT_REAL 3    /* the real component only */
#define BAILOUT_IMAG 4    /* the imaginary component only */
#define BAILOUT_BOTH 5    /* both components at once */
/* Regular polygons, measured by their apothem -- the distance from the centre
 * to the middle of a side, which is what the bailout value gives, so that a
 * polygon and the circle of the same bailout touch at the sides rather than
 * at the corners. The orientation is part of the shape: a triangle turned by
 * 120 degrees is the same triangle, so 0, 90 and -90 are three different
 * ones, while a hexagon repeats every 60 and an octagon every 45. */
#define BAILOUT_TRIANGLE0 6
#define BAILOUT_TRIANGLE90 7
#define BAILOUT_TRIANGLEM90 8
#define BAILOUT_HEXAGON0 9
#define BAILOUT_HEXAGON90 10
#define BAILOUT_OCTAGON 11
#define BAILOUTMODES 12
#define BAILOUT_MAXSIDES 8

struct fractal_context {
    number_t pre, pim;
    number_t bre, bim;
    const struct formula *currentformula;
#ifdef USE_SFFE
    sffe *userformula;
    sffe *userinitial;
#endif
    number_t angle;
    int periodicity;
    unsigned int maxiter;
    number_t bailout;
    /* Which shape the bailout tests against; see the BAILOUT_ constants
     * below. Zero is the circle every version before this one used, so a
     * context never told otherwise behaves as it always did. */
    int bailoutmode;
    /* Derived from bailoutmode and bailout by set_fractalc, never saved and
     * set by nothing else: the outward normals of the bailout polygon and its
     * apothem. Precomputed because the alternative is a sine and a cosine per
     * side on every iteration. bailoutsides is zero for the shapes that need
     * none of this. */
    number_t bailoutnx[BAILOUT_MAXSIDES], bailoutny[BAILOUT_MAXSIDES];
    number_t bailoutapothem;
    int bailoutsides;
    OutColormodeType coloringmode;
    int incoloringmode;
    int intcolor, outtcolor;
    //MSUMMO BEGIN HACK
    number_t incolorspeed, outcolorspeed;
    int incolorfun, outcolorfun;
    int incolorshift, outcolorshift;
    //MSUMMO HACK 20220409
    int pndefault;
    int newtonmodesffe;
    number_t newtonconvergence;
    //MSUMMO END HACK
    int mandelbrot;
    int plane;
    int version;
    int range;
    float windowwidth, windowheight;
    vinfo s;
    vrect rs;
    number_t sin, cos;
    int slowmode; /* 1 in case we want to be exact, not fast */
    /*values temporary filled by set_fractal_context */
    iterationfunc calculate[2];
    number_t periodicity_limit;
    struct palette *palette; /*fractal's palette */
};
typedef struct fractal_context fractal_context;
typedef struct {
    double y0, k, kk, y0k;
} symmetry2;

struct symmetryinfo2 {
    number_t xsym, ysym;
    int nsymmetries;
    symmetry2 *symmetry;
    number_t xmul, ymul, xdist, ydist;
};
#define STARTZERO 1
#define JULIA_BTRACE 2
#define MANDEL_BTRACE 4
#ifdef USE_SFFE
#define SFFE_FRACTAL 8
#endif

#ifdef USE_SFFE
/* Set when the user formula calls randsc or randscq. Boundary tracing walks
 * the edge of a region, finds one colour all the way round and fills the
 * inside without computing it -- true of a fractal, false of a noise field,
 * where the inside is whatever the noise says. Left on, some pixels are
 * filled rather than computed and the noise is wrong there; and since the
 * flag that enables tracing differs between mandelbrot and julia mode, the
 * two modes disagreed about which pixels those were. */
extern int sffe_formula_noise;
#define BTRACE_NOISE_OK (!sffe_formula_noise)
#else
#define BTRACE_NOISE_OK 1
#endif

#define BTRACEOK                                                               \
    ((cformula.flags & (2 << cfractalc.mandelbrot)) &&                         \
     !cfractalc.incoloringmode && BTRACE_NOISE_OK &&                           \
     cfractalc.coloringmode != OutColormodeType::ColOut_potential)
#define my_rotate(f, x, y)                                                        \
    {                                                                          \
        number_t tmp;                                                          \
        tmp = (x) * (f).cos - (y) * (f).sin;                                   \
        y = (x) * (f).sin + (y) * (f).cos;                                     \
        x = tmp;                                                               \
    }
#define rotateback(f, x, y)                                                    \
    {                                                                          \
        number_t tmp;                                                          \
        tmp = (x) * (f).cos + (y) * (f).sin;                                   \
        y = -(x) * (f).sin + (y) * (f).cos;                                    \
        x = tmp;                                                               \
    }

#ifdef USE_SFFE
void sffe_setlocal(fractal_context *c);
#endif

extern struct symmetryinfo2 cursymmetry;
extern struct fractal_context cfractalc;
extern struct formula cformula;
extern struct palette cpalette;
extern struct image cimage;

#ifdef STATISTICS
/*This is an statistics variables printed from various parts
 *of XaoS.
 */
extern int nadded2, nsymmetry2, nskipped2;
extern int tocalculate2, avoided2, frames2;
extern int ncalculated2, ninside2;
extern int niter2, niter1;
extern int nperi;

extern int iters2, guessed2, unguessed2, total2;

#endif

void set_formula(fractal_context *, int);
void set_fractalc(fractal_context *, struct image *img);
void fractalc_resize_to(fractal_context *, float, float);
void update_view(fractal_context *);
void free_fractalc(fractal_context *);
fractal_context *make_fractalc(const int, float, float);
void speed_test(fractal_context *, struct image *img);
unsigned int calculateswitch(number_t x1, number_t y1, number_t x2, number_t y2,
                             int periodicity);

/* needs struct formula */
#include "formulas.h"

#endif /* FRACTAL_H */
