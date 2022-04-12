﻿/*
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
#define OUTCOLORING 11
#define TCOLOR 15
#define COLORFUN 10

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
    int coloringmode, incoloringmode;
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

#define BTRACEOK                                                               \
    ((cformula.flags & (2 << cfractalc.mandelbrot)) &&                         \
     !cfractalc.incoloringmode && cfractalc.coloringmode != 7)
#define rotate(f, x, y)                                                        \
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
