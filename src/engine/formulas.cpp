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
/* Hello reader!
 * code you are looking at is dangerous for both you and your hardware! PLEASE
 * CLOSE THIS FILE UNLESS YOU REALLY KNOW WHAT YOU ARE DOING.
 *
 * Main purpose of this file is to generate optimal caluclation loops for
 * various formulas/algorithms. It heavily includes docalc.c - set of
 * caluclation loops, that then uses macros instead of formulas. This lets me
 * to change calculation loops easily. At the other hand it looks very ugly.
 * You have been warned :)
 */

// Some help can be read below about line 700. :-)

#include <climits>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cstdio>

#include "config.h"
#include "number_math.h"
#include "phist.h"
#include "cmplx.h"
#include "filter.h"
#include "fractal.h"
#include "julia.h"
#include "ui_helper.h"
#include "xthread.h"
#ifndef M_PI
#define M_PI 3.1415
#endif

#ifdef USE_SFFE
#include "sffe.h"
#endif

const char *const incolorname[] = {"0",
                                   "zmag",
                                   "Decomposition-like",
                                   "real/imag",
                                   "abs(abs(c)-abs(r))",
                                   "cos(mag)",
                                   "mag*cos(real^2)",
                                   "sin(real^2-imag^2)",
                                   "atan(real*imag*creal*cimag)",
                                   "squares",
                                   "real/mag",
                                   "max(|real|,|imag|)",
                                   "|real|+|imag|",
                                   "min(|real|,|imag|)",
                                   "|z-c|",
                                   "|z*c|",
                                   "angle(z)-angle(c)",
                                   "real*imag",
                                   "sin(real)*sin(imag)",
                                   "sign(imag)",
                                   "frac(mag)",
                                   "log(mag)",
                                   "True-color",
                                   NULL};


const char *const outcolorname[] = {"iter",
                                    "iter+real",
                                    "iter+imag",
                                    "iter+real/imag",
                                    "iter+real+imag+real/imag",
                                    "binary decomposition",
                                    "biomorphs",
                                    "potential",
                                    "color decomposition",
                                    "smooth",
                                    "smooth log",
                                    "iter+angle",
                                    "iter+log(mag)",
                                    "iter+real*imag",
                                    "max(|real|,|imag|)",
                                    "iter banded",
                                    "|real|-|imag|",
                                    "True-color",
                                    NULL};

const char *const tcolorname[] = {
    "black",
    "re*im sin(re^2) angle",
    "sin(re) sin(im) sin(square)",
    "hsv",
    "hsv2",
    "cos(re^c) cos(im^2) cos(square)",
    "abs(re^2) abs(im^2) abs(square)",
    "re*im re*re im*im",
    "abs(im*cim) abs(re*cre) abs(re*cim)",
    "abs(re*im-csqr) abs(re^2-csqr) abs(im^2-csqr)",
    "angle angle2 angle",
    "Disable truecolor colouring",
    "simple red (for education purposes)",
    "simple blue (for education purposes)",
    NULL};

const char *const colorfun[] = {
    "ident",
    "log",
    "atan",
    "sin",
    "cos",
    "square",
    "cube",
    "square root",
    "cube root",
    NULL
};

#define SHIFT 8
#define SMUL 256

/* Bailout shapes.
 *
 * The iteration stops when z leaves a region, and which region that is has
 * always been a circle: |z|^2 > bailout. The shape shows in the picture --
 * the bands outside the set follow it -- so the others here are as much a
 * drawing tool as a numerical one.
 *
 * rp and ip are the squared components, so every test is in squares and only
 * the diamond needs a root. The mode is one value for a whole image, so the
 * branch below predicts perfectly and costs about what the load of
 * cfractalc.bailout beside it costs. */
/* True while the point is still inside, which is the sense BTEST is used in.
 *
 * Takes the components as well as their squares: the squares are what the
 * first few shapes want and are already to hand, while a polygon needs the
 * signs, since a side faces one way only. */
static inline int bailout_inside(number_t zre, number_t zim, number_t rp,
                                 number_t ip)
{
    const number_t b = cfractalc.bailout;
    switch (cfractalc.bailoutmode) {
        case BAILOUT_SQUARE:
            return rp < b && ip < b;
        case BAILOUT_DIAMOND:
            return rp + ip + 2 * nsqrt(rp * ip) < b;
        case BAILOUT_REAL:
            return rp < b;
        case BAILOUT_IMAG:
            return ip < b;
        case BAILOUT_BOTH:
            return rp < b || ip < b;
        case BAILOUT_CIRCLE:
            return rp + ip < b;
        default: {
            /* A regular polygon: still inside while every side's outward
             * normal projects the point short of the apothem. Three to eight
             * multiply-adds, with the normals worked out in set_fractalc. */
            const number_t a = cfractalc.bailoutapothem;
            for (int k = 0; k < cfractalc.bailoutsides; k++)
                if (zre * cfractalc.bailoutnx[k] +
                        zim * cfractalc.bailoutny[k] >=
                    a)
                    return 0;
            return 1;
        }
    }
}


/* The quantity a bailout shape actually tests, and what it tests it against.
 *
 * bailout_inside answers a yes or no; smooth colouring needs the number behind
 * it. It interpolates between one pass and the next by asking where the
 * quantity crossed the threshold, and that only means anything if it is the
 * same quantity the escape was decided on. Interpolating on the modulus while
 * escaping on a hexagon is what left seams along the sides of the shape.
 *
 * Every mode compares one scalar against one threshold, so both are here.
 * For the circle the scalar is the squared modulus and the threshold the
 * bailout, which is what smooth colouring always used -- so nothing moves for
 * the shape every version had. */
static inline number_t bailout_measure(number_t zre, number_t zim, number_t rp,
                                       number_t ip)
{
    switch (cfractalc.bailoutmode) {
        case BAILOUT_SQUARE:
            return rp > ip ? rp : ip;
        case BAILOUT_DIAMOND:
            return rp + ip + 2 * nsqrt(rp * ip);
        case BAILOUT_REAL:
            return rp;
        case BAILOUT_IMAG:
            return ip;
        case BAILOUT_BOTH:
            return rp < ip ? rp : ip;
        case BAILOUT_CIRCLE:
            return rp + ip;
        default: {
            number_t furthest = 0;
            for (int k = 0; k < cfractalc.bailoutsides; k++) {
                number_t side = zre * cfractalc.bailoutnx[k] +
                                zim * cfractalc.bailoutny[k];
                if (side > furthest)
                    furthest = side;
            }
            return furthest;
        }
    }
}

static inline number_t bailout_threshold(void)
{
    /* The polygons are measured by their apothem and against a projection,
     * which is a length; the rest compare squares. */
    return cfractalc.bailoutmode > BAILOUT_BOTH ? cfractalc.bailoutapothem
                                                : cfractalc.bailout;
}

#ifndef less_than_4
#define less_than_0(x) ((x) < 0)
#define less_than_4(x) ((x) < cfractalc.bailout)
/* Newton mode stops when successive iterates stop moving. The limit was
 * written into the macro, so the Newton convergence the menu offers was
 * set, saved and reloaded while nothing ever read it. */
#define not_converged(n) ((n) > cfractalc.newtonconvergence)
#define abs_less_than(x, y) (myabs(x) < y)
#define greater_than(x, y) ((x) > (y))
#endif

#define PERIINOUTPUT()                                                         \
    STAT(nperi++; ninside2++);                                                 \
    return (cpalette.pixels[0])

#define OUTOUTPUT()                                                            \
    STAT(niter2 += iter);                                                      \
    return (color_output(zre, zim, iter))
    //return (cfractalc.coloringmode == OutColormodeType::ColOut_iter          \                                 \
                ? cpalette.pixels[(iter % (cpalette.size - 1)) + 1]            \
                : color_output(zre, zim, iter))
#define INOUTPUT()                                                             \
    STAT(niter1 += iter; ninside2++);                                          \
    return (cfractalc.incoloringmode                                           \
                ? incolor_output(zre, zim, pre, pim, iter)                     \
                : cpalette.pixels[0])

#define OUTPUT()                                                               \
    if (iter >= (unsigned int)cfractalc.maxiter) {                             \
        if (cfractalc.incoloringmode == INCOLORING_TRUECOLOR)                                    \
            return (                                                           \
                truecolor_output(zre, zim, pre, pim, cfractalc.intcolor, 1));  \
        INOUTPUT();                                                            \
    } else {                                                                   \
        if (cfractalc.coloringmode == OutColormodeClass::ColOut_True_color)                                      \
            return (                                                           \
                truecolor_output(zre, zim, pre, pim, cfractalc.outtcolor, 0)); \
        OUTOUTPUT();                                                           \
    }

#define SMOOTHOUTPUT()                                                         \
    {                                                                          \
        PRESMOOTH;                                                             \
        {                                                                      \
            /* Enough to keep a logarithm away from zero and no more. It was   \
             * a millionth, which is nothing beside a bailout of four and      \
             * everything beside a Newton convergence of a millionth; below    \
             * one it now follows the threshold down. At the bailouts anyone   \
             * uses it is the millionth it always was. */                      \
            number_t nudge = SMOOTHTHRESHOLD < 1                               \
                                 ? (number_t)SMOOTHTHRESHOLD * (number_t)1e-6  \
                                 : (number_t)0.000001;                         \
            zre += nudge;                                                      \
            szmag += nudge;                                                    \
        }                                                                      \
        iter = (int)(((cfractalc.maxiter - iter) * 256 +                       \
                      log((double)(SMOOTHTHRESHOLD / (szmag))) /               \
                          log((double)((zre) / (szmag))) * 256));              \
        if (cfractalc.coloringmode == OutColormodeClass::ColOut_smooth_log) { \
           iter = log(iter) * ((cpalette.size - 1))/log(cfractalc.maxiter * 256) + 1;  \
        }\
        /* The colouring speed, the shift and the colouring function. Every    \
         * other outside mode meets color_precalc on its way out of            \
         * color_output; this one returns a pixel of its own and so skipped    \
         * it, which is why neither speed nor shift did anything at all while  \
         * smooth or smooth log was chosen. What is held here is the same      \
         * fixed point the others use -- the count times 256, plus a fraction  \
         * -- so it goes through unchanged.                                    \
         *                                                                     \
         * The negative case is handled as color_output handles it: a shift    \
         * can carry the value below zero, and the modulus of a negative is    \
         * not the wrap that is wanted. */                                     \
        {                                                                      \
            number_t smoothed = (number_t)iter;                                \
            int wrapped;                                                       \
            int cycle = (int)(((unsigned int)(cpalette.size - 1)) << 8);       \
            color_precalc(smoothed, 0);                                        \
            wrapped = (int)smoothed;                                           \
            if (wrapped < 0) {                                                 \
                wrapped = cycle - ((-wrapped) % cycle) - 1;                    \
                if (wrapped < 0)                                               \
                    wrapped = 0;                                               \
            }                                                                  \
            iter = (unsigned int)wrapped;                                      \
        }                                                                      \
        iter %= ((unsigned int)(cpalette.size - 1)) << 8;                      \
                                                                               \
        if ((cpalette.type & (C256 | SMALLITER)) || !(iter & 255))             \
            return (cpalette.pixels[1 + (iter >> 8)]);                         \
        {                                                                      \
            unsigned int i1, i2;                                               \
            i1 = cpalette.pixels[1 + (iter >> 8)];                             \
            if ((iter >> 8) == (unsigned int)(cpalette.size - 2))              \
                i2 = cpalette.pixels[1];                                       \
            else                                                               \
                i2 = cpalette.pixels[2 + (iter >> 8)];                         \
            iter &= 255;                                                       \
            return (interpoltype(cpalette, i2, i1, iter));                     \
        }                                                                      \
    }


/* 2021-02-09 MSUMMO calculate color functions */
const number_t V_MAX = DBL_MAX / SMUL;
static void color_f_normalize(number_t &v) {
    if (v < -V_MAX) {
        v = -V_MAX;
    } else if (v > V_MAX) {
        v = V_MAX;
    }
}

static float color_f_func(number_t v, const float speed, const int func, const int shift) {
    color_f_normalize(v);
    switch (func) {
    case 1:
        // LOG
        if (v > 0)
            v = logl(v + 1);
        else
            v = -logl(-v + 1);
        break;
    case 2:
        // ATAN
        v = atanl(v);
        break;
    case 3:
        // SIN
        v = sinl(v);
        break;
    case 4:
        // COS
        v = cosl(v);
        break;
    case 5:
        // SQR
        if (v > 0)
            v = v * v;
        else
            v = -(v * v);
        color_f_normalize(v);
        break;
    case 6:
        // CUBE
        v = v * v * v;
        color_f_normalize(v);
        break;
    case 7:
        // SQRT
        if (v > 0)
            v = sqrtl(v);
        else
            v = -sqrtl(-v);
        break;
    case 8:
        // CBRT
        if (v > 0)
            v = cbrtl(v);
        else
            v = -cbrtl(-v);
        break;
    }
    v = v * speed;
    v += shift;
    color_f_normalize(v);
    return v;
}

static void color_precalc(number_t &iter, int inset) {
    const float speed = inset ? cfractalc.incolorspeed : cfractalc.outcolorspeed;
    const int func = inset ? cfractalc.incolorfun : cfractalc.outcolorfun;
    const int shift = inset ? cfractalc.incolorshift : cfractalc.outcolorshift;
    iter = (color_f_func(iter / SMUL, speed, func, shift) * SMUL);
}
static void color_precalc(number_t &zre, number_t &zim, number_t &pre,
                          number_t &pim, int inset) {
    const float speed = inset ? cfractalc.incolorspeed : cfractalc.outcolorspeed;
    const int func = inset ? cfractalc.incolorfun : cfractalc.outcolorfun;
    const int shift = inset ? cfractalc.incolorshift : cfractalc.outcolorshift;
    zre = color_f_func(zre, speed, func, shift);
    zim = color_f_func(zim, speed, func, shift);
    pre = color_f_func(pre, speed, func, shift);
    pim = color_f_func(pim, speed, func, shift);
}


/* 2009-07-30 JB Langston:
 * Fixing bug #3: HSV modes are completely black when compiled with GCC 4...
 * Removed  qualifier from hsv_to_rgb declaration.  macro is
 * defined to __attribute__((__const__)), on which I found some more details
 * here: http://unixwiz.net/techtips/gnu-c-attributes.html#const.  Apparently
 * this should never be used with a function that takes a pointer or relies on
 * side-effects, and hsv_to_rgb does both.  Therefore, it should never have
 * been declared this way in the first place.
 */

static inline void hsv_to_rgb(int h, int s, int v, int *red, int *green,
                              int *blue)
    /**/;
static inline void hsv_to_rgb(int h, int s, int v, int *red, int *green,
                              int *blue)
{
    int hue;
    int f, p, q, t;

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

static unsigned int truecolor_output(number_t zre, number_t zim, number_t pre,
                                     number_t pim, int mode, int inset);

static unsigned int truecolor_output(number_t zre, number_t zim, number_t pre,
                                     number_t pim, int mode, int inset)
{
    /* WARNING: r and b fields are swapped for HISTORICAL REASONS (BUG :),
     * in other words: use r for blue and b for red. */
    /*
     * MSUMMO: Fixed bug!
     */
    int b = 0, g = 0, r = 0, w = 0;
    color_precalc(zre, zim, pre, pim, inset);
    switch (mode) {
    case 0:
        break;
    case 1:
        r = (int)((sinl(atan2l(zre, zim) * 20) + 1) *
                   127);
        w = (int)((sinl(zim / zre)) * 127);
        b = (int)((int)(zre * zim));
        g = (int)((sinl((zre * zre) / 2) + 1) * 127);
        break;
    case 2:
        if (!inset) {
            b = (int)((sinl(zre * 2) + 1) * 127);
            g = (int)((sinl(zim * 2) + 1) * 127);
            r = (int)((sinl((zim * zim + zre * zre) / 2) + 1) * 127);
        } else {
            b = (int)((sinl(zre * 50) + 1) * 127);
            g = (int)((sinl(zim * 50) + 1) * 127);
            r = (int)((sinl((zim * zim + zre * zre) * 50) + 1) *
                       127);
        }
        w = (int)((sinl(zim / zre)) * 127);
        break;
    case 3:
        if (inset)
            hsv_to_rgb((int)(atan2l(zre, zim) * 256 / M_PI),
                       (int)((sinl((zre * 50)) + 1) * 128),
                       (int)((sinl((zim * 50)) + 1) * 128), &b, &g,
                       &r);
        else
            hsv_to_rgb((int)(atan2l(zre, zim) * 256 / M_PI),
                       (int)((sinl(zre) + 1) * 128),
                       (int)((sinl(zim) + 1) * 128), &b, &g, &r);
        break;
    case 4:
        if (inset)
            hsv_to_rgb(
                (int)(sinl((zre * zre + zim * zim) * 0.1) * 256),
                (int)(sinl(atan2l(zre, zim) * 10) * 128 +
                       128),
                (int)((sinl((zre + zim) * 10)) * 65 + 128), &b, &g,
                &r);
        else
            hsv_to_rgb(
                (int)(sinl((zre * zre + zim * zim) * 0.01) * 256),
                (int)(sinl(atan2l(zre, zim) * 10) * 128 +
                       128),
                (int)((sinl((zre + zim) * 0.3)) * 65 + 128), &b, &g,
                &r);
        break;
    case 5: {
        if (!inset) {
            b = (int)(cosl(myabs(zre * zre)) * 128) + 128;
            g = (int)(cosl(myabs(zre * zim)) * 128) + 128;
            r = (int)(cosl(myabs(zim * zim + zre * zre)) * 128) +
                128;
        } else {
            b = (int)(cosl(myabs(zre * zre) * 10) * 128) + 128;
            g = (int)(cosl(myabs(zre * zim) * 10) * 128) + 128;
            r = (int)(cosl(myabs(zim * zim + zre * zre) * 10) *
                       128) +
                128;
        }
    } break;
    case 6: {
        if (!inset) {
            b = (int)(zre * zim * 64);
            g = (int)(zre * zre * 64);
            r = (int)(zim * zim * 64);
        } else
            b = (int)(zre * zim * 256);
        g = (int)(zre * zre * 256);
        r = (int)(zim * zim * 256);
    } break;
    case 7: {
        if (!inset) {
            b = (int)((zre * zre + zim * zim - pre * pre - pim * pim) * 16);
            g = (int)((zre * zre * 2 - pre * pre - pim * pim) * 16);
            r = (int)((zim * zim * 2 - pre * pre - pim * pim) * 16);
        } else {
            b = (int)((zre * zre + zim * zim - pre * pre - pim * pim) *
                       256);
            g = (int)((zre * zre * 2 - pre * pre - pim * pim) * 256);
            r = (int)((zim * zim * 2 - pre * pre - pim * pim) * 256);
        }
    } break;
    case 8: {
        if (!inset) {
            b = (int)((myabs(zim * pim)) * 64);
            g = (int)((myabs(zre * pre)) * 64);
            r = (int)((myabs(zre * pim)) * 64);
        } else {
            b = (int)((myabs(zim * pim)) * 256);
            g = (int)((myabs(zre * pre)) * 256);
            r = (int)((myabs(zre * pim)) * 256);
        }
    } break;
    case 9: {
        if (!inset) {
            b = (int)((myabs(zre * zim - pre * pre - pim * pim)) * 64);
            g = (int)((myabs(zre * zre - pre * pre - pim * pim)) * 64);
            r = (int)((myabs(zim * zim - pre * pre - pim * pim)) * 64);
        } else {
            b = (int)((myabs(zre * zim - pre * pre - pim * pim)) * 256);
            g = (int)((myabs(zre * zre - pre * pre - pim * pim)) * 256);
            r = (int)((myabs(zim * zim - pre * pre - pim * pim)) * 256);
        }
    } break;
    case 10: {
        b = (int)(atan2l(zre, zim) * 128 / M_PI) + 128;
        g = (int)(atan2l(zre, zim) * 128 / M_PI) + 128;
        r = (int)(atan2l(zim, zre) * 128 / M_PI) + 128;
    } break;
        // case 11 is for disabling truecolor mode
    case 12: {
        r = 255;
        g = 0;
        b = 0;
        w = 50;
    } break;
    case 13: {
        b = 255;
        g = 0;
        r = 0;
        w = 0;
    } break;
    }

    b += w;
    g += w;
    r += w;
    if (b < 0)
        b = 0;
    else if (b > 255)
        b = 255;
    if (g < 0)
        g = 0;
    else if (g > 255)
        g = 255;
    if (r < 0)
        r = 0;
    else if (r > 255)
        r = 255;

    switch (cpalette.type) {
    case GRAYSCALE:
        return ((unsigned int)(b * 76 + g * 151 + r * 29) *
                    (cpalette.end - cpalette.start) >>
                16) +
               cpalette.start;
    case TRUECOLOR:
    case TRUECOLOR24:
    case TRUECOLOR16:
        b >>= cpalette.info.truec.bprec;
        g >>= cpalette.info.truec.gprec;
        r >>= cpalette.info.truec.rprec;
        return ((b << cpalette.info.truec.bshift) +
                (g << cpalette.info.truec.gshift) +
                (r << cpalette.info.truec.rshift) +
                cpalette.info.truec.alpha);
    }

    return cpalette.pixels[inset];
}

static unsigned int color_output(number_t zre, number_t zim, unsigned int iter);
static unsigned int color_output(number_t zre, number_t zim, unsigned int iter)
{
    int i;
    number_t i_f;
    iter <<= SHIFT;
    i_f = iter;
    if (cfractalc.coloringmode.ColorMode != OutColormodeType::ColOut_iter) {
        color_f_normalize(zre);
        color_f_normalize(zim);
    }

    switch (cfractalc.coloringmode.ColorMode) {
        case OutColormodeType::ColOut_smooth:
        case OutColormodeType::ColOut_smooth_log:
            break;
        case OutColormodeType::ColOut_iter_plus_real: /* real */
            i_f = (iter + zre * SMUL);
            break;
        case OutColormodeType::ColOut_iter_plus_imag: /* imag */
            i_f = (iter + zim * SMUL);
            break;
        case OutColormodeType::ColOut_iter_plus_real_div_imag: /* real / imag */
            i_f = (iter + (zre / zim) * SMUL);
            break;
        case OutColormodeType::ColOut_iter_plus_real_plus_imag_plus_real_div_imag: /* all of the above */
            i_f = (iter + (zre + zim + zre / zim) * SMUL);
            break;
        case OutColormodeType::ColOut_binary_decomposition:
            if (zim > 0)
                i_f = ((cfractalc.maxiter << SHIFT) - iter);
            break;
        case OutColormodeType::ColOut_biomorphs:
            if (myabs(zim) < 2.0 || myabs(zre) < 2.0)
                i_f = ((cfractalc.maxiter << SHIFT) - iter);
            break;
        case OutColormodeType::ColOut_potential:
            zre = zre * zre + zim * zim;
            i_f = (sqrtl(logl(zre) / i_f) * 256 * 256);
            break;
        case OutColormodeType::ColOut_color_decomposition:
            i_f = ((atan2l(zre, zim) / (M_PI + M_PI) + 0.75) *
                   20000);
            break;
        case OutColormodeType::ColOut_iter_plus_angle:
            /* the angle it left by, blended into the count: the bands acquire
             * a twist and the escape directions show as spokes */
            i_f = (iter + (natan2(zim, zre) / (N_PI + N_PI) + (number_t)0.5) *
                              SMUL * 8);
            break;
        case OutColormodeType::ColOut_iter_plus_log_mag:
            /* how far past the bailout it went, which is a cheaper smoothing
             * than the log of a log and bands differently */
            i_f = (iter + nlog(zre * zre + zim * zim + 1) * SMUL * 4);
            break;
        case OutColormodeType::ColOut_iter_plus_real_times_imag:
            /* biomorphs made continuous: the same quantity that the biomorph
             * test thresholds, added to the count instead */
            i_f = (iter + zre * zim * SMUL);
            break;
        case OutColormodeType::ColOut_max_real_imag:
            /* which side it left by, and how far, with no count at all */
            i_f = ((myabs(zre) > myabs(zim) ? myabs(zre) : myabs(zim)) * SMUL *
                   200);
            break;
        case OutColormodeType::ColOut_iter_banded:
            /* the count folded into eight bands, which shows the shape of the
             * level sets rather than their number */
            i_f = ((number_t)(((unsigned int)i_f >> SHIFT) % 8) * SMUL * 8);
            break;
        case OutColormodeType::ColOut_abs_real_minus_abs_imag:
            /* how lopsided the escape point is, with no count: the level sets
             * of the difference rather than of the modulus */
            i_f = ((myabs(zre) - myabs(zim)) * SMUL * 200);
            break;
        default:
            break;
    }

    color_precalc(i_f, 0);
    i = (int)i_f;

    if (i < 0) {
        i = (((unsigned int)(cpalette.size - 1)) << 8) -
            ((-i) % (((unsigned int)(cpalette.size - 1) << 8))) - 1;
        if (i < 0)
            i = 0;
    }
    iter = ((unsigned int)i) % ((cpalette.size - 1) << 8);
    if ((cpalette.type & (C256 | SMALLITER)) || !(iter & 255))
        return (cpalette.pixels[1 + (iter >> 8)]);
    {
        unsigned int i1, i2;

        i1 = cpalette.pixels[1 + (iter >> 8)];

        if ((int)(iter >> 8) == cpalette.size - 2)
            i2 = cpalette.pixels[1];
        else
            i2 = cpalette.pixels[2 + (iter >> 8)];

        iter &= 255;
        return (interpoltype(cpalette, i2, i1, iter));
    }
}

static unsigned int incolor_output(number_t zre, number_t zim, number_t pre,
                                   number_t pim, unsigned int iter);

static unsigned int incolor_output(number_t zre, number_t zim, number_t pre,
                                   number_t pim, unsigned int iter)
{
    int i = iter;
    number_t i_f = iter;
    if (cfractalc.coloringmode.ColorMode) {
        color_f_normalize(zre);
        color_f_normalize(zim);
        color_f_normalize(pre);
        color_f_normalize(pim);
    }

    switch (cfractalc.incoloringmode) {
    case 1: /* zmag */
        i_f = (((zre * zre + zim * zim) *
                    (number_t)(cfractalc.maxiter >> 1) * SMUL +
                SMUL));
        break;
    case 2: /* real */
        i_f = ((
            (atan2l(zre, zim) / (M_PI + M_PI) + 0.75) *
            20000));
        break;
    default:
        break;
    case 3: /* real / imag */
        i_f = (100 + (zre / zim) * SMUL * 10);
        break;
    case 4:
        zre = myabs(zre);
        zim = myabs(zim);
        pre = myabs(pre);
        pre = myabs(pim);
        i_f += (myabs(pre - zre) * 256 * 64);
        i_f += (myabs(pim - zim) * 256 * 64);
        break;
    case 5:
        if (((int)((zre * zre + zim * zim) * 10)) % 2)
            i_f = (cosl((zre * zim * pre * pim)) * 256 * 256);
        else
            i_f = (sinl((zre * zim * pre * pim)) * 256 * 256);
        break;
    case 6:
        i_f = ((zre * zre + zim * zim) * cosl((zre * zre)) * 256 *
               256);
        break;
    case 7:
        i_f = (sinl((zre * zre - zim * zim)) * 256 * 256);
        break;
    case 8:
        i_f = (atanl((zre * zim * pre * pim)) * 256 * 64);
        break;
    case 9:
        if ((abs((int)(zre * 40)) % 2) ^ (abs((int)(zim * 40)) % 2))
            i_f = ((
                (atan2l(zre, zim) / (M_PI + M_PI) + 0.75) *
                20000));
        else
            i_f = ((
                (atan2l(zim, zre) / (M_PI + M_PI) + 0.75) *
                20000));
        break;
    case 10: /* real/mag: the cosine of the angle. Where atan2 has a seam
              * wherever the angle wraps, this closes on itself. */
        i_f = ((zre / nsqrt(zre * zre + zim * zim) + 1) * 20000);
        break;
    case 11: /* max(|real|,|imag|): the square norm, so the bands are squares
              * where zmag's are circles */
        i_f = ((myabs(zre) > myabs(zim) ? myabs(zre) : myabs(zim)) *
                   (number_t)(cfractalc.maxiter >> 1) * SMUL +
               SMUL);
        break;
    case 12: /* |real|+|imag|: the same bands turned by an eighth of a turn */
        i_f = ((myabs(zre) + myabs(zim)) * (number_t)(cfractalc.maxiter >> 1) *
                   SMUL +
               SMUL);
        break;
    case 13: /* min(|real|,|imag|): near zero along either axis, so rays */
        i_f = ((myabs(zre) < myabs(zim) ? myabs(zre) : myabs(zim)) *
                   (number_t)cfractalc.maxiter * SMUL +
               SMUL);
        break;
    case 14: /* |z-c|: how far the orbit settled from the parameter, which is
              * one thing that tells the bulbs apart */
        i_f = (nsqrt((zre - pre) * (zre - pre) + (zim - pim) * (zim - pim)) *
                   (number_t)(cfractalc.maxiter >> 1) * SMUL +
               SMUL);
        break;
    case 15: /* |z*c|: the two moduli against each other */
        i_f = (nsqrt((zre * zre + zim * zim) * (pre * pre + pim * pim)) *
                   (number_t)(cfractalc.maxiter >> 1) * SMUL +
               SMUL);
        break;
    case 16: /* angle(z) - angle(c) */
        i_f = ((natan2(zim, zre) - natan2(pim, pre)) / (N_PI + N_PI) +
               (number_t)0.75) *
              20000;
        break;
    case 17: /* real*imag: a saddle, four lobes about the axes */
        i_f = (zre * zim * (number_t)cfractalc.maxiter * SMUL + SMUL);
        break;
    case 18: /* sin(real)*sin(imag): a grid laid over where the orbit settled */
        i_f = ((nsin(zre) * nsin(zim) + 1) * 20000);
        break;
    case 19: /* sign(imag): two flat tones. A decomposition of the inside,
              * where binary decomposition does the same to the outside. */
        i_f = (zim > 0 ? 20000 : 45000);
        break;
    case 20: /* frac(mag): contours at every eighth of the modulus */
        {
            number_t bands = nsqrt(zre * zre + zim * zim) * 8;
            i_f = ((bands - nfloor(bands)) * 60000);
        }
        break;
    case 21: /* log(mag): zmag with the contrast moved onto the small values,
              * which is where an attracting orbit spends its time */
        i_f = (nlog(zre * zre + zim * zim + (number_t)1 / 1000000) * 8000 +
               30000);
        break;
    };
    color_precalc(i_f, 1);
    i = (int) i_f;

    if (i < 0) {
        i = (((unsigned int)(cpalette.size - 1)) << 8) -
            ((-i) % (((unsigned int)(cpalette.size - 1) << 8))) - 1;
        if (i < 0)
            i = 0;
    }
    iter = ((unsigned int)i) % ((cpalette.size - 1) << 8);

    if ((cpalette.type & (C256 | SMALLITER)) || !(iter & 255))
        return (cpalette.pixels[1 + ((unsigned int)iter >> 8)]);
    {
        unsigned int i1, i2;
        i1 = cpalette.pixels[1 + ((unsigned int)iter >> 8)];
        if (((unsigned int)iter >> 8) == (unsigned int)(cpalette.size - 2))
            i2 = cpalette.pixels[1];
        else
            i2 = cpalette.pixels[2 + ((unsigned int)iter >> 8)];
        iter &= 255;
        return (interpoltype(cpalette, i2, i1, iter));
    }
}

#define VARIABLES
#define INIT
#define UNCOMPRESS
#define PRETEST 0
#define FORMULA                                                                \
    zim = (zim * zre) * 2 + pim;                                               \
    zre = rp - ip + pre;                                                       \
    ip = zim * zim;                                                            \
    rp = zre * zre;

/* Some help for the brave ones. :-)
 *
 * Mandelbrot's original formula is z=z^2+c which means
 * z[next]=z[previous]^2+c.
 * Here c is the pixel coordinates from the screen and z[0] is usually 0
 * (if not perturbation was added.)
 * In the following code z[previous] is described by (zre;zim)
 * and z[next] will also be zre and zim.
 * c is described by (pre;pim).
 * Finally rp and ip are helper variables, mostly for checking the bailout
 * (which usually means abs(z)>=4, see BTEST).
 *
 * Both basic operations and some other functions (c_mul, c_pow3, ...) can
 * be used. For a "detailed" description refer to ../include/complex.h.
 *
 * If you add/modify fractals, please note that struct formula_formulas
 * (at line cca. 1300) should be also edited for proper initialization
 * and for menu entries. However it is not self-explanatory, just copy-paste
 * existing tables and give it a try.
 *
 * Finally, please also edit the calculateswitch function and
 * the nmformulas constant (at the end of this file).
 *
 * -- Zoltan, 2009-07-30
 */

#define BTEST bailout_inside(zre, zim, rp, ip)
#define SMOOTH
#define SCALC smand_calc
#define SPERI smand_peri
#define CALC mand_calc
#define PERI mand_peri
#define JULIA mand_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define PRETEST 0
#define FORMULA                                                                \
    rp = zre * (rp - 3 * ip);                                                  \
    zim = zim * (3 * zre * zre - ip) + pim;                                    \
    zre = rp + pre;                                                            \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define BTEST bailout_inside(zre, zim, rp, ip)
#define SMOOTH
#define SCALC smand3_calc
#define SPERI smand3_peri
#define CALC mand3_calc
#define PERI mand3_peri
#define JULIA mand3_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define UNCOMPRESS
#define VARIABLES number_t br, tmp;
#define FORMULA                                                                \
    br = zre + zre + pre - 2;                                                  \
    tmp = zre * zim;                                                           \
    zre = rp - ip + pre - 1;                                                   \
    ip = zim + zim + pim;                                                      \
    zim = tmp + tmp + pim;                                                     \
    tmp = 1 / (br * br + ip * ip);                                             \
    rp = (zre * br + zim * ip) * tmp;                                          \
    ip = (zim * br - zre * ip) * tmp;                                          \
    zre = (rp + ip) * (rp - ip);                                               \
    zim = rp * ip;                                                             \
    zim += zim;                                                                \
    rp = zre - 1;                                                              \
    ip = zim * zim;                                                            \
    rp = zre * zre;
#define BTEST                                                                  \
    (rp + ip < (number_t)100 * 100 &&                                          \
     (rp - 2 * zre + ip) > 0.04 / cfractalc.bailout - 1)
#define POSTCALC                                                               \
    if (rp - 2 * zre + ip > 0.04 / cfractalc.bailout - 1) {                    \
        zre *= 0.08 / cfractalc.bailout, zim *= 0.08 / cfractalc.bailout;      \
        if (iter)                                                              \
            iter = cfractalc.maxiter - iter + 1;                               \
    }
#define CALC magnet_calc
#define PERI magnet_peri
#define SCALC smagnet_calc
#define SPERI smagnet_peri
#define SMOOTH
#define PRESMOOTH                                                              \
    szmag /= 100 * 100 / 4;                                                    \
    zre = (rp + ip) / (100 * 100 * 4);
#define JULIA magnet_julia
#define RANGE 4
#define RPIP
#include "docalc.h"

#define UNCOMPRESS
#define VARIABLES number_t inre, inim, tmp1, tmp2, dnre, nmre, dnim;
#define INIT                                                                   \
    inre = pre * pre - pim * pim - pre - pre - pre;                            \
    inim = pre * pim;                                                          \
    inim = inim + inim - pim - pim - pim;
#define FORMULA                                                                \
    tmp1 = rp - ip;                                                            \
    tmp2 = zre * pre - zim * pim - zre;                                        \
    dnre =                                                                     \
        tmp1 + tmp1 + tmp1 + tmp2 + tmp2 + tmp2 - zre - zre - zre + inre + 3;  \
    tmp1 = zre * ip;                                                           \
    nmre = zre * rp - tmp1 - tmp1 - tmp1 + tmp2 + tmp2 + tmp2 + inre + 2;      \
    tmp1 = zre * zim;                                                          \
    tmp2 = zre * pim + zim * pre - zim;                                        \
    dnim = tmp1 + tmp1 + tmp1 + tmp1 + tmp1 + tmp1 + tmp2 + tmp2 + tmp2 -      \
           zim - zim - zim + inim;                                             \
    tmp1 = zim * rp;                                                           \
    zim = tmp1 + tmp1 + tmp1 - zim * ip + tmp2 + tmp2 + tmp2 + inim;           \
    zre = nmre;                                                                \
    ip = dnim;                                                                 \
    tmp1 = 1 / (dnre * dnre + ip * ip);                                        \
    rp = (zre * dnre + zim * ip) * tmp1;                                       \
    ip = (zim * dnre - zre * ip) * tmp1;                                       \
    zre = (rp + ip) * (rp - ip);                                               \
    zim = rp * ip;                                                             \
    zim += zim;                                                                \
    ip = zim * zim;                                                            \
    rp = zre * zre;
#define BTEST                                                                  \
    (rp + ip < (number_t)100 * 100 &&                                          \
     (rp - 2 * zre + ip) > 0.04 / cfractalc.bailout - 1)
#define POSTCALC                                                               \
    if (rp - 2 * zre + ip > 0.04 / cfractalc.bailout - 1) {                    \
        zre *= 0.08 / cfractalc.bailout, zim *= 0.08 / cfractalc.bailout;      \
        if (iter)                                                              \
            iter = cfractalc.maxiter - iter + 1;                               \
    }
#define CALC magnet2_calc
#define PERI magnet2_peri
#define SCALC smagnet2_calc
#define SPERI smagnet2_peri
#define SMOOTH
#define PRESMOOTH                                                              \
    szmag /= 100 * 100 / 4;                                                    \
    zre = (rp + ip) / (100 * 100 * 4);
#define JULIA magnet2_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    rp = rp * rp - 6 * rp * ip + ip * ip + pre;                                \
    zim = 4 * zre * zre * zre * zim - 4 * zre * ip * zim + pim;                \
    zre = rp;                                                                  \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC smand4_calc
#define SPERI smand4_peri
#define CALC mand4_calc
#define PERI mand4_peri
#define JULIA mand4_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES number_t t;
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    c_pow4(zre, zim, rp, ip);                                                  \
    c_mul(zre, zim, rp, ip, t, zim);                                           \
    zre = t + pre;                                                             \
    zim += pim;                                                                \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC smand5_calc
#define SPERI smand5_peri
#define CALC mand5_calc
#define PERI mand5_peri
#define JULIA mand5_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES number_t t;
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    c_pow3(zre, zim, rp, ip);                                                  \
    c_mul(rp, ip, rp, ip, t, zim);                                             \
    zre = t + pre;                                                             \
    zim += pim;                                                                \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC smand6_calc
#define SPERI smand6_peri
#define CALC mand6_calc
#define PERI mand6_peri
#define JULIA mand6_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES number_t t;
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    c_pow3(zre, zim, rp, ip);                                                  \
    c_pow3(rp, ip, t, zim);                                                    \
    zre = t + pre;                                                             \
    zim += pim;                                                                \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC smand9_calc
#define SPERI smand9_peri
#define CALC mand9_calc
#define PERI mand9_peri
#define JULIA mand9_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    zim = zre * zim + zim / 2 + pim;                                           \
    zre = (rp - ip + zre) / 2 + pre;                                           \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC strice_calc
#define SPERI strice_peri
#define CALC trice_calc
#define PERI trice_peri
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES number_t zor, zoi;
/* For some reason Cat's Eye renders as an empty circle unless the bailout
 * is slightly more than 4.  It was first observed in 2009 on Mac OS X but
 * more recently started happening on other operating systems. I suspect
 * a compiler bug, but I haven't been able to figure out exactly what's
 * happening.  I can work around it by subtracting a small amount from the
 * magnitude before performing the bailout test.
 */
/* Its own bailout test, so its own quantity to smooth on: the shape chosen
 * for the others means nothing here. */
#define CUSTOMSAVEZMAG szmag = rp + ip
#define PRESMOOTH zre = rp + ip
#define SMOOTHTHRESHOLD cfractalc.bailout
#define BTEST less_than_4(rp + ip - 0.00000001)
#define FORMULA                                                                \
    c_div(pre, pim, zre, zim, rp, ip);                                         \
    c_div(zre, zim, pre, pim, zor, zoi);                                       \
    zre = zor + rp;                                                            \
    zim = zoi + ip;                                                            \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC scatseye_calc
#define SPERI scatseye_peri
#define CALC catseye_calc
#define PERI catseye_peri
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    zim = (zim * zre) * (-2.0) + pim;                                          \
    zre = rp - ip + pre;                                                       \
    ip = zim * zim;                                                            \
    rp = zre * zre;
#define SMOOTH
#define SCALC smbar_calc
#define SPERI smbar_peri
#define CALC mbar_calc
#define PERI mbar_peri
#define JULIA mbar_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES
#define INIT                                                                   \
    rp = zre;                                                                  \
    zre = pre;                                                                 \
    pre = rp;                                                                  \
    ip = zim;                                                                  \
    zim = pim;                                                                 \
    pim = ip;
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    rp = ip - rp + zre;                                                        \
    ip = zim - 2 * zre * zim;                                                  \
    c_mul(rp, ip, pre, pim, zre, zim);                                         \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC smlambda_calc
#define SPERI smlambda_peri
#define CALC mlambda_calc
#define PERI mlambda_peri
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES number_t zre1, zim1, zre2, zim2;
#define INIT                                                                   \
    zre1 = zre;                                                                \
    zim1 = zim;
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    zre2 = zre;                                                                \
    zim2 = zim;                                                                \
    zim = (zim * zre) * 2 + pim + zim1;                                        \
    zre = rp - ip + pre + zre1;                                                \
    zre1 = zre2;                                                               \
    zim1 = zim2;                                                               \
    ip = zim * zim;                                                            \
    rp = zre * zre;
#define SMOOTH
#define SCALC smanowar_calc
#define CALC manowar_calc
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES number_t zre1, zim1;
#define INIT                                                                   \
    zre1 = pre;                                                                \
    zim1 = pim;
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    zim = (zim * zre) * 2 + zim1;                                              \
    zre = rp - ip + zre1;                                                      \
    zre1 = zre1 / 2 + zre;                                                     \
    zim1 = zim1 / 2 + zim;                                                     \
    ip = zim * zim;                                                            \
    rp = zre * zre;
#define SMOOTH
#define SCALC sspider_calc
#define CALC spider_calc
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES
#define INIT                                                                   \
    if ((zre == pre) && (zim == pim)) {                                        \
        pre = 0.5;                                                             \
        pim = 0.8660254;                                                       \
    }                                                                          \
    if (pim < 0)                                                               \
        pim = (-pim);                                                          \
    if (((pim * zre - pre * zim) < 0) || (zim < 0)) {                          \
        zre = 2 * pre + 2;                                                     \
        zim = 2 * pim;                                                         \
    }
#define BTEST ((pim * zre + (1 - pre) * zim) < pim)
#define FORMULA                                                                \
    zre = 2 * zre;                                                             \
    zim = 2 * zim;                                                             \
    if ((pim * zre - pre * zim) > pim)                                         \
        zre = zre - 1;                                                         \
    if (zim > pim) {                                                           \
        zim = zim - pim;                                                       \
        zre = zre - pre;                                                       \
    }
#define CALC sier_calc
#define RANGE 2
#include "docalc.h"

#define VARIABLES
#define INIT                                                                   \
    if ((zre == pre) && (zim == pim)) {                                        \
        pre = 0.5;                                                             \
        pim = 0.8660254;                                                       \
    }                                                                          \
    if (pim < 0)                                                               \
        pim = (-pim);                                                          \
    if (((pim * zre - pre * zim) < 0) || (zim < 0)) {                          \
        zre = 2 * pre + 2;                                                     \
        zim = 2 * pim;                                                         \
    }
#define BTEST ((pim * zre + (1 - pre) * zim) < pim)
#define FORMULA                                                                \
    zre = 1.6180339 * zre;                                                     \
    zim = 1.6180339 * zim;                                                     \
    if ((pim * zre - pre * zim) > pim * 0.6180339)                             \
        zre = zre - 0.6180339;                                                 \
    if (zim > pim * 0.6180339) {                                               \
        zim = zim - pim * 0.6180339;                                           \
        zre = zre - pre * 0.6180339;                                           \
    }
#define CALC goldsier_calc
#define RANGE 2
#include "docalc.h"

#define VARIABLES
#define INIT
#define BTEST (zre * zre + zim * zim < 1)
#define FORMULA                                                                \
    zre = 3 * zre;                                                             \
    zim = 3 * zim;                                                             \
    if ((zim - 2) * (zim - 2) + zre * zre < 1)                                 \
        zim = zim - 2;                                                         \
    if ((zim + 2) * (zim + 2) + zre * zre < 1)                                 \
        zim = zim + 2;                                                         \
    if ((zim - 1) * (zim - 1) + (zre - 1.7320508) * (zre - 1.7320508) < 1) {   \
        zim = zim - 1;                                                         \
        zre = zre - 1.7320508;                                                 \
    }                                                                          \
    if ((zim + 1) * (zim + 1) + (zre - 1.7320508) * (zre - 1.7320508) < 1) {   \
        zim = zim + 1;                                                         \
        zre = zre - 1.7320508;                                                 \
    }                                                                          \
    if ((zim - 1) * (zim - 1) + (zre + 1.7320508) * (zre + 1.7320508) < 1) {   \
        zim = zim - 1;                                                         \
        zre = zre + 1.7320508;                                                 \
    }                                                                          \
    if ((zim + 1) * (zim + 1) + (zre + 1.7320508) * (zre + 1.7320508) < 1) {   \
        zim = zim + 1;                                                         \
        zre = zre + 1.7320508;                                                 \
    }
#define CALC circle7_calc
#define RANGE 2
#include "docalc.h"

#define VARIABLES
#define INIT
#define BTEST (zre * zre + zim * zim < 0.6944444444444445)
#define FORMULA                                                                \
  const double rd = 0.25;                                       \
  zre = 3 * zre;                                                               \
  zim = 3 * zim;                                                               \
  if ((zim + 2) * (zim + 2) + zre * zre < rd) {                                \
    zim = zim + 2;                                                             \
  }                                                                            \
  if ((zim + 1.7320508075688772) * (zim + 1.7320508075688772) +                \
          (zre + 1) * (zre + 1) <                                              \
      rd) {                                                                    \
    zim = zim + 1.7320508075688772;                                            \
    zre = zre + 1;                                                             \
  }                                                                            \
  if ((zim + 1) * (zim + 1) +                                                  \
          (zre + 1.7320508075688772) * (zre + 1.7320508075688772) <            \
      rd) {                                                                    \
    zim = zim + 1;                                                             \
    zre = zre + 1.7320508075688772;                                            \
  }                                                                            \
  if ((zre + 2) * (zre + 2) + zim * zim < rd) {                                \
    zre = zre + 2;                                                             \
  }                                                                            \
  if ((zim - 1) * (zim - 1) +                                                  \
          (zre + 1.7320508075688772) * (zre + 1.7320508075688772) <            \
      rd) {                                                                    \
    zim = zim - 1;                                                             \
    zre = zre + 1.7320508075688772;                                            \
  }                                                                            \
  if ((zim - 1.7320508075688772) * (zim - 1.7320508075688772) +                \
          (zre + 1) * (zre + 1) <                                              \
      rd) {                                                                    \
    zim = zim - 1.7320508075688772;                                            \
    zre = zre + 1;                                                             \
  }                                                                            \
  if ((zim - 2) * (zim - 2) + zre * zre < rd) {                                \
    zim = zim - 2;                                                             \
  }                                                                            \
  if ((zim - 1.7320508075688772) * (zim - 1.7320508075688772) +                \
          (zre - 1) * (zre - 1) <                                              \
      rd) {                                                                    \
    zim = zim - 1.7320508075688772;                                            \
    zre = zre - 1;                                                             \
  }                                                                            \
  if ((zim - 1) * (zim - 1) +                                                  \
          (zre - 1.7320508075688772) * (zre - 1.7320508075688772) <            \
      rd) {                                                                    \
    zim = zim - 1;                                                             \
    zre = zre - 1.7320508075688772;                                            \
  }                                                                            \
  if ((zre - 2) * (zre - 2) + zim * zim < rd) {                                \
    zre = zre - 2;                                                             \
  }                                                                            \
  if ((zim + 1) * (zim + 1) +                                                  \
          (zre - 1.7320508075688772) * (zre - 1.7320508075688772) <            \
      rd) {                                                                    \
    zim = zim + 1;                                                             \
    zre = zre - 1.7320508075688772;                                            \
  }                                                                            \
  if ((zim + 1.7320508075688772) * (zim + 1.7320508075688772) +                \
          (zre - 1) * (zre - 1) <                                              \
      rd) {                                                                    \
    zim = zim + 1.7320508075688772;                                            \
    zre = zre - 1;                                                             \
  }
#define CALC clock_calc
#define RANGE 2
#include "docalc.h"

#define VARIABLES
#define INIT
#define BTEST less_than_4((rp + ip) / 4.0)
#define FORMULA                                                                \
    if (less_than_0(zre)) {                                                    \
        rp = zre + 1;                                                          \
    } else {                                                                   \
        rp = zre - 1;                                                          \
    }                                                                          \
    if (less_than_0(zim)) {                                                    \
        ip = zim + 1;                                                          \
    } else {                                                                   \
        ip = zim - 1;                                                          \
    }                                                                          \
    c_mul(rp, ip, pre, pim, zre, zim);                                         \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define CALC symbarn_calc
#define JULIA symbarn_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES
#define INIT                                                                   \
    if ((zre == pre) && (zim == pim)) {                                        \
        pre = 1;                                                               \
        pim = 1;                                                               \
    }                                                                          \
    if (pre < 0)                                                               \
        pre = (-pre);                                                          \
    if (pim < 0)                                                               \
        pim = (-pim);                                                          \
    if ((zre < 0) || (zre > pre)) {                                            \
        zre = pre / 2;                                                         \
        zim = pim / 2;                                                         \
    }                                                                          \
    if ((zim < 0) || (zim > pim)) {                                            \
        zre = pre / 2;                                                         \
        zim = pim / 2;                                                         \
    }
#define BTEST                                                                  \
    ((zre < pre / 3) || (zre > 2 * pre / 3) || (zim < pim / 3) ||              \
     (zim > 2 * pim / 3))
#define FORMULA                                                                \
    zre = 3 * zre;                                                             \
    zim = 3 * zim;                                                             \
    if (zre > 2 * pre)                                                         \
        zre = zre - 2 * pre;                                                   \
    else if (zre > pre)                                                        \
        zre = zre - pre;                                                       \
    if (zim > 2 * pim)                                                         \
        zim = zim - 2 * pim;                                                   \
    else if (zim > pim)                                                        \
        zim = zim - pim;
#define CALC carpet_calc
#define RANGE 2
#include "docalc.h"

#define VARIABLES
#define BTEST                                                                  \
    ((((1.5 * zre + 0.8660254038 * zim) > 0.8660254038) ||                     \
      ((0.8660254038 * zim - 1.5 * zre) > 0.8660254038) || (zim < (-0.5))) &&  \
     (((1.5 * zre + 0.8660254038 * zim) < -0.8660254038) ||                    \
      ((0.8660254038 * zim - 1.5 * zre) < -0.8660254038) || (zim > 0.5)))
#define FORMULA                                                                \
    zre = 3 * zre;                                                             \
    zim = 3 * zim;                                                             \
    if ((0.2886751346 * zim - 0.5 * zre) > 0.0) {                              \
        if ((0.2886751346 * zim + 0.5 * zre) > 0.0) {                          \
            zim = zim - 2.0;                                                   \
        } else {                                                               \
            if (zim > 0) {                                                     \
                zre = zre + 1.732050808;                                       \
                zim = zim - 1.0;                                               \
            } else {                                                           \
                zre = zre + 1.732050808;                                       \
                zim = zim + 1.0;                                               \
            }                                                                  \
        }                                                                      \
    } else {                                                                   \
        if ((0.2886751346 * zim + 0.5 * zre) < 0.0) {                          \
            zim = zim + 2.0;                                                   \
        } else {                                                               \
            if (zim > 0) {                                                     \
                zre = zre - 1.732050808;                                       \
                zim = zim - 1.0;                                               \
            } else {                                                           \
                zre = zre - 1.732050808;                                       \
                zim = zim + 1.0;                                               \
            }                                                                  \
        }                                                                      \
    }
#define CALC koch_calc
#define RANGE 2
#include "docalc.h"

#define VARIABLES number_t zre1, zim1;
#define INIT                                                                   \
    pim = fabs(pim);                                                           \
    zre = pre;                                                                 \
    zim = pim;
#define BTEST                                                                  \
    (!((zre < 0) && (zim > 0) &&                                               \
       (-1.0 * zre + 1.732050808 * zim < 1.732050808)))
#define FORMULA                                                                \
    zre1 = 1.5 * zre - 0.866 + 0.866 * zim;                                    \
    zim1 = -1.5 + 1.5 * zim - 0.866 * zre;                                     \
    zre = zre1;                                                                \
    zim = zim1;
#define CALC hornflake_calc
#define RANGE 2
#include "docalc.h"

#define VARIABLES
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    if (less_than_0(zre)) {                                                    \
        rp = zre + 1;                                                          \
    } else {                                                                   \
        rp = zre - 1;                                                          \
    }                                                                          \
    c_mul(rp, zim, pre, pim, zre, zim);                                        \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC sbarnsley1_calc
#define CALC barnsley1_calc
#define JULIA barnsley1_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    if (less_than_0(zre * pim + zim * pre)) {                                  \
        rp = zre + 1;                                                          \
    } else {                                                                   \
        rp = zre - 1;                                                          \
    }                                                                          \
    c_mul(rp, zim, pre, pim, zre, zim);                                        \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC sbarnsley2_calc
#define CALC barnsley2_calc
#define JULIA barnsley2_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    if (!less_than_0(-zre)) {                                                  \
        zim = 2 * zre * zim + pim * zre;                                       \
        zre = rp - ip - 1 + pre * zre;                                         \
    } else {                                                                   \
        zim = 2 * zre * zim;                                                   \
        zre = rp - ip - 1;                                                     \
    }                                                                          \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define SMOOTH
#define SCALC sbarnsley3_calc
#define CALC barnsley3_calc
#define JULIA barnsley3_julia
#define RANGE 2
#define RPIP
#include "docalc.h"

#define VARIABLES number_t n, sqrr, sqri, zre1, zim1;
#define INIT                                                                   \
    sqri = zim * zim, n = zre, zre = pre, pre = n, n = zim, zim = pim,         \
    pim = n, n = (number_t)1;
#define BTEST not_converged(n)
#define FORMULA                                                                \
    zre1 = zre;                                                                \
    zim1 = zim;                                                                \
    n = zim * zim;                                                             \
    sqri = zre * zre;                                                          \
    sqrr = sqri - n;                                                           \
    sqri = n + sqri;                                                           \
    n = 0.3333333333 / ((sqri * sqri));                                        \
    zim = (0.66666666) * zim - (zre + zre) * zim * n + pim;                    \
    zre = (0.66666666) * zre + (sqrr)*n + pre;                                 \
    zre1 -= zre;                                                               \
    zim1 -= zim;                                                               \
    n = zre1 * zre1 + zim1 * zim1;
/* Newton converges rather than escaping, so what smooth colouring interpolates
 * on is the step length falling past the convergence limit rather than the
 * modulus climbing past the bailout. The formula is the same either way: where
 * between the two passes did the quantity cross. */
#define CUSTOMSAVEZMAG szmag = n
#define PRESMOOTH zre = n
#define SMOOTHTHRESHOLD cfractalc.newtonconvergence
#define SMOOTH
#define SCALC snewton_calc
#define CALC newton_calc
#include "docalc.h"

#define VARIABLES number_t n, sqrr, sqri, zre1, zim1;
#define INIT                                                                   \
    sqri = zim * zim, n = zre, zre = pre, pre = n, n = zim, zim = pim,         \
    pim = n, n = (number_t)1;
#define BTEST not_converged(n)
#define FORMULA                                                                \
    zre1 = zre;                                                                \
    zim1 = zim;                                                                \
    sqrr = zre * zre;                                                          \
    sqri = zim * zim;                                                          \
    n = sqri + sqrr;                                                           \
    n = 1 / ((n * n * n));                                                     \
    zim = (0.25) * zim * (3 + (sqri - 3 * sqrr) * n) + pim;                    \
    zre = (0.25) * zre * (3 + (sqrr - 3 * sqri) * n) + pre;                    \
    zre1 -= zre;                                                               \
    zim1 -= zim;                                                               \
    n = zre1 * zre1 + zim1 * zim1;
/* Newton converges rather than escaping, so what smooth colouring interpolates
 * on is the step length falling past the convergence limit rather than the
 * modulus climbing past the bailout. The formula is the same either way: where
 * between the two passes did the quantity cross. */
#define CUSTOMSAVEZMAG szmag = n
#define PRESMOOTH zre = n
#define SMOOTHTHRESHOLD cfractalc.newtonconvergence
#define SMOOTH
#define SCALC snewton4_calc
#define CALC newton4_calc
#include "docalc.h"

#define VARIABLES number_t zpr, zip;
#define SAVEVARIABLES number_t szpr, szip;
#define SAVE szpr = zpr, szip = zip;
#define RESTORE zpr = szpr, zip = szip;
#define INIT zpr = zip = (number_t)0;
#define BTEST bailout_inside(zre, zim, rp, ip)
#define FORMULA                                                                \
    rp = rp - ip + pre + pim * zpr;                                            \
    ip = 2 * zre * zim + pim * zip;                                            \
    zpr = zre, zip = zim;                                                      \
    zre = rp;                                                                  \
    zim = ip;                                                                  \
    rp = zre * zre, ip = zim * zim;
#define SMOOTH
#define SCALC sphoenix_calc
#define SPERI sphoenix_peri
#define CALC phoenix_calc
#define PERI phoenix_peri
#define RPIP
#include "docalc.h"

#define VARIABLES number_t tr, ti, zpr, zpm, rp1, ip1;
#define INIT                                                                   \
    zpr = zpm = 0, tr = zre, zre = pre, pre = tr, tr = zim, zim = pim,         \
    pim = tr, tr = 1;
#define BTEST less_than_4(zpr *zpr + zpm * zpm)
#define FORMULA                                                                \
    rp1 = zre;                                                                 \
    ip1 = zim;                                                                 \
    c_pow3(zre, zim, tr, ti);                                                  \
    c_add(tr, ti, zpr, zpm, zre, zim);                                         \
    zpr = rp1 + pre;                                                           \
    zpm = ip1 + pim;
#define CALC octo_calc
#define SCALC socto_calc
#define SMOOTH
#define CUSTOMSAVEZMAG szmag = zpr * zpr + zpm * zpm
#define PRESMOOTH zre = zpr * zpr + zpm * zpm
#include "docalc.h"

#define VARIABLES number_t yre, yim, re1tmp, re2tmp, im1tmp;
#define BTEST (rp + ip < 9 || (yre * yre + yim * yim) < 4 * (rp + ip))
#define INIT                                                                   \
    yre = pre;                                                                 \
    yim = pim;
#define FORMULA                                                                \
    re1tmp = zre;                                                              \
    re2tmp = yre;                                                              \
    im1tmp = zim;                                                              \
    zre = re1tmp + yre;                                                        \
    zim = im1tmp + yim;                                                        \
    yre = (re1tmp * re2tmp - im1tmp * yim);                                    \
    yim = (re1tmp * yim + re2tmp * im1tmp);                                    \
    rp = zre * zre;                                                            \
    ip = zim * zim;
#define CALC beryl_calc
#define PERI beryl_peri
#define RANGE 2
#define RPIP
#include "docalc.h"

#ifdef USE_SFFE

// Parser is not thread safe so each thread needs its own instance
thread_local bool sffe_formula_valid = false;
thread_local sffe *sffe_formula_local = NULL;
thread_local bool sffe_initial_valid = false;
thread_local sffe *sffe_initial_local = NULL;
/* sffe_z is defined beside sffe_position in sffe_cmplx_gsl.cpp: the
 * figures follow the value from one pass to the next and so have to read
 * it, the way the noise reads the position. */
thread_local cmplx sffe_c, sffe_n, sffe_x;
/* Set when the user formula calls randsc or randscq; read by BTRACEOK. */
int sffe_formula_noise = 0;

// Copy the formula from the main parser to this thread's local parser
// Possibly initializing the parser if this is the first time
void sffe_setmine(void *data, struct taskinfo * /*task*/, int /*r1*/,
                  int /*r2*/)
{
    fractal_context *c = (fractal_context *)data;
    if (!sffe_formula_local) {
        sffe_formula_local = sffe_alloc();
        /* p, p1, p2 ... are asked for as they are met; see sffe_resolve_p */
        sffe_formula_local->resolve = sffe_resolve_p;
        sffe_regvar(&sffe_formula_local, &sffe_z, "z");
        sffe_regvar(&sffe_formula_local, &sffe_c, "c");
        sffe_regvar(&sffe_formula_local, &sffe_n, "n");
        sffe_regvar(&sffe_formula_local, &sffe_x, "x");
    }
    if (c->userformula->expression) {
        if (sffe_parse(&sffe_formula_local, c->userformula->expression) == 0)
            sffe_formula_valid = true;
        else
            sffe_formula_valid = false;
        /* Boundary tracing fills a region it found one colour around without
         * computing it. True of a fractal, false of a noise field, so a
         * formula that calls randsc or randscq turns it off. */
        sffe_formula_noise = sffe_formula_valid &&
                             sffe_uses_noise(sffe_formula_local);
    }

    if (!sffe_initial_local) {
        sffe_initial_local = sffe_alloc();
        sffe_initial_local->resolve = sffe_resolve_p;
        sffe_regvar(&sffe_initial_local, &sffe_c, "c");
        sffe_regvar(&sffe_initial_local, &sffe_n, "n");
        /* On the initial parser, which is what "usrformInit" is parsed by.
         * x was registered on the formula parser instead, and z on neither, so
         * an initialization naming either was refused as an unknown variable
         * by the thread that had to compute it -- while the dialog, whose
         * parser has both, accepted the formula and gave no sign that it was
         * then thrown away. In an initialization z is where z would have
         * started had there been no initialization, which is the same thing x
         * is; "z" alone therefore says exactly what saying nothing says. */
        sffe_regvar(&sffe_initial_local, &sffe_x, "x");
        sffe_regvar(&sffe_initial_local, &sffe_z, "z");
    }
    if (c->userinitial->expression) {
        if (sffe_parse(&sffe_initial_local, c->userinitial->expression) == 0)
            sffe_initial_valid = true;
        else
            sffe_initial_valid = false;
    }

    /* Both are parsed: work out which places they read and how deep a history
     * that needs. Once per formula, not once per pass. */
    sffe_ptaps_build(sffe_formula_valid ? sffe_formula_local : NULL,
                     sffe_initial_valid ? sffe_initial_local : NULL);
}

// Tell all threads copy the formula into their local parser
void sffe_setlocal(fractal_context *c)
{
    xth_function(sffe_setmine, c, nthreads);
    xth_sync();
}

/* z and n are set whether or not there is an initialization to read them: an
 * initialization that named either used to read whatever the pixel before it
 * had left behind, which is to say it read the order the pixels happened to
 * be computed in. Before the first pass z is where z starts, which is what x
 * is as well, so an initialization of "z" says what saying nothing says. */
#define INIT                                                                   \
    if (pndef) {                                                               \
            sffe_phist_reset(pre, pim);                                        \
    } else {                                                                   \
            sffe_phist_reset(0, 0);                                            \
    }                                                                          \
    cmplxset(sffe_c, pre, pim);                                                \
    cmplxset(sffe_x, zre, zim);                                                \
    cmplxset(sffe_z, zre, zim);                                                \
    cmplxset(sffe_n, 1, 0);                                                    \
    sffe_iteration = 0;                                                        \
    sffe_maxiter = maxit;                                                      \
    if (sffe_initial_valid)                                                    \
        sffe_z = sffe_eval(sffe_initial_local);                                \
    n = INFINITY;
//#define SAVE cmplxset(pZ,real(Z),imag(Z));
//#define PRETEST 0
#define VARIABLES number_t n; bool newtok = cfractalc.newtonmodesffe;          \
bool pndef = cfractalc.pndefault; unsigned int maxit                           \
    = (unsigned int)cfractalc.maxiter;
/* The newton emulation measures the step from one pass to the next, so it
 * takes where z was before the two lines below put where it now is in its
 * place. */
#define FORMULA                                                                \
    if (sffe_formula_valid)                                                    \
        sffe_z = sffe_eval(sffe_formula_local);                                \
    sffe_phist_push(zre, zim);                                                 \
    if (newtok) {                                                              \
        pre = zre;                                                             \
        pim = zim;                                                             \
    }                                                                          \
    zre = real(sffe_z);                                                        \
    zim = imag(sffe_z);                                                        \
    if (newtok) {                                                              \
        n = iter > 1 ? zre*zre - 2*zre*pre + zim*zim - 2*pim*zim + pre*pre + pim*pim : 0; \
    }                                                                          \
    cmplxset(sffe_n, maxit - iter + 1, 0);                                     \
    sffe_iteration += 1;
/* The escape half goes through the shape selector like every other
 * escape-time formula; the convergence half above it is a different test
 * and a shape means nothing there. This one was missed when the other
 * fifteen were converted, because it spells the sum out rather than
 * using rp and ip -- so choosing a bailout mode did nothing at all for a
 * user formula, which is the one place it is most wanted. */
#define BTEST newtok ?                                                         \
not_converged(n) \
    : bailout_inside(zre, zim, zre * zre, zim * zim)
#define CALC sffe_calc
#define JULIA sffe_julia
#define SCALC ssffe_calc
#define SMOOTH
/* The user formula does not keep rp and ip -- its bailout test spells the two
 * squares out -- so the places smooth colouring reads its quantity have to be
 * told where to get it. And it has two bailout tests, one for each mode: the
 * shape when it escapes, the step length when it converges, so the quantity
 * and the threshold follow whichever is in force.
 *
 * n starts at infinity to force the first pass, and an orbit that converges on
 * that very pass would leave a logarithm with nothing to work on; one stands
 * in for it, as the built-in Newtons start theirs. */
#define CUSTOMSAVEZMAG                                                         \
    szmag = newtok ? (n < (number_t)1e30 ? n : (number_t)1)                    \
                   : bailout_measure(zre, zim, zre * zre, zim * zim)
#define PRESMOOTH                                                              \
    zre = newtok ? n : bailout_measure(zre, zim, zre * zre, zim * zim)
#define SMOOTHTHRESHOLD                                                        \
    (newtok ? cfractalc.newtonconvergence : bailout_threshold())
#include "docalc.h"
#endif

static const symmetrytype sym6[] = {{0, 1.73205080758}, {0, -1.73205080758}};

static const symmetrytype sym8[] = {{0, 1}, {0, -1}};

static const symmetrytype sym16[] = {{0, 1},        {0, -1},
                                     {0, 0.414214}, {0, -0.414214},
                                     {0, 2.414214}, {0, -2.414214}};

const struct formula formulas[] = {
    {                           /* 0 */
     FORMULAMAGIC,
     mand_calc,
     mand_peri,
     smand_calc,
     smand_peri,
     mand_julia,
     {"Mandelbrot", "Julia"},
     "mandel",
     /*{0.5, -2.0, 1.25, -1.25}, */
     {-0.75, 0.0, 2.5, 2.5},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* 1 */
     FORMULAMAGIC,
     mand3_calc,
     mand3_peri,
     smand3_calc,
     smand3_peri,
     mand3_julia,
     {"Mandelbrot^3", "Julia^3"},
     "mandel3",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     1, 1, 0.0, 0.0,
     {
      {0, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {0, 0, 0, NULL},
      {0, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, 0, 0, NULL},
      {0, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {0, 0, 0, NULL},
      {0, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* 2 */
     FORMULAMAGIC,
     mand4_calc,
     mand4_peri,
     smand4_calc,
     smand4_peri,
     mand4_julia,
     {"Mandelbrot^4", "Julia^4"},
     "mandel4",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* 3 */
     FORMULAMAGIC,
     mand5_calc,
     mand5_peri,
     smand5_calc,
     smand5_peri,
     mand5_julia,
     {"Mandelbrot^5", "Julia^5"},
     "mandel5",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     1, 1, 0.0, 0.0,
     {
      {0, 0, 2, sym8},
      {INT_MAX, 0, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {0, 0, 2, sym8},
      {0, 0, 2, sym8},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, 0, 2, sym8},
      {0, 0, 2, sym8},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {0, 0, 2, sym8},
      {0, 0, 2, sym8},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* 4 */
     FORMULAMAGIC,
     mand6_calc,
     mand6_peri,
     smand6_calc,
     smand6_peri,
     mand6_julia,
     {"Mandelbrot^6", "Julia^6"},
     "mandel6",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* 5 */
     FORMULAMAGIC,
     newton_calc,
     NULL,
     snewton_calc,
     NULL,
     NULL,
     {"Newton", "Newton julia?"},
     "newton",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     0, 1, 1.0199502202048319698, 0.0,
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     STARTZERO,
     },
    {                           /* formula added by Andreas Madritsch *//* 6 */
     FORMULAMAGIC,
     newton4_calc,
     NULL,
     snewton4_calc,
     NULL,
     NULL,
     {"Newton^4", "Newton^4 julia?"},
     "newton4",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     0, 1, 1.0199502202048319698, 0.0,
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     STARTZERO,
     },
    {                           /* 7 */
     FORMULAMAGIC,
     barnsley1_calc,
     NULL,
     sbarnsley1_calc,
     NULL,
     barnsley1_julia,
     {"Barnsley1 Mandelbrot", "Barnsley1"},
     "barnsley",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     0, 0, -0.6, 1.1,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     STARTZERO | MANDEL_BTRACE,
     },
    {                           /* formula added by Andreas Madritsch *//* 8 */
     FORMULAMAGIC,
     barnsley2_calc,
     NULL,
     sbarnsley2_calc,
     NULL,
     barnsley2_julia,
     {"Barnsley2 Mandelbrot", "Barnsley2"},
     "barnsley2",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 5.5},
     0, 0, -0.6, 1.1,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     STARTZERO | MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete *//* 9 */
     FORMULAMAGIC,
     barnsley3_calc,
     NULL,
     sbarnsley3_calc,
     NULL,
     barnsley3_julia,
     {"Barnsley3 Mandelbrot", "Barnsley3"},
     "barnsley3",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 3.5},
     0, 0, 0.0, 0.4,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     STARTZERO | MANDEL_BTRACE,
     },
    {                           /* 10 */
     FORMULAMAGIC,
     octo_calc,
     /*octo_peri, */ NULL,
     socto_calc,
     /*socto_peri, */ NULL,
     NULL,
     {"Octo", "Octo"},
     "octo",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     0, 1, 0.0, 0.0,
     {
      {0, 0, 6, sym16},
      {INT_MAX, 0, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {0, 0, 0, NULL},
      {0, 0, 6, sym16},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, 0, 6, sym16},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {0, 0, 6, sym16},
      {0, 0, 6, sym16},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE | STARTZERO,
     },
    {                           /* 11 */
     FORMULAMAGIC,
     phoenix_calc,
     phoenix_peri,
     sphoenix_calc,
     sphoenix_peri,
     NULL,
     {"MandPhoenix", "Phoenix"},
     "phoenix",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     1, 0, 0.56667000000000001, -0.5,
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* 12 */
     FORMULAMAGIC,
     magnet_calc,
     magnet_peri,
     smagnet_calc,
     smagnet_peri,
     magnet_julia,
     {"Magnet", "Magnet"},
     "magnet",
     /*{3, 0, 2.2, -2.2}, */
     {1.5, 0.0, 3.0, 4.4},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     STARTZERO,
     },
    {                           /* formula added by Andreas Madritsch *//* 13 */
     FORMULAMAGIC,
     magnet2_calc,
     magnet2_peri,
     smagnet2_calc,
     smagnet2_peri,
     magnet2_julia,
     {"Magnet2", "Magnet2"},
     "magnet2",
     /*{3, 0, 2.2, -2.2}, */
     {1.0, 0.0, 3.0, 3.2},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     STARTZERO,
     },
    {                           /* formula added by Arpad Fekete *//* 14 */
     FORMULAMAGIC,
     trice_calc,
     trice_peri,
     strice_calc,
     strice_peri,
     NULL,
     {"Triceratops", "Triceratops Julia"},
     "trice",
     {0.0, 0.0, 2.5, 4.5},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete *//* 15 */
     FORMULAMAGIC,
     catseye_calc,
     catseye_peri,
     scatseye_calc,
     scatseye_peri,
     NULL,
     {"Catseye", "Catseye Julia"},
     "catseye",
     {0.0, 0.0, 2.5, 4.5},
     1, 1, 0.0, 0.0,
     {
      {0, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, 0, 0, NULL},
      {0, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {0, 0, 0, NULL},
      {0, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, 0, 0, NULL},
      {0, 0, 0, NULL},
      {0, 0, 0, NULL},
      {0, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /*formula added by Arpad Fekete *//* 16 */
     /*in Gnofract4d from mathworld.wolfram.com */
     FORMULAMAGIC,
     mbar_calc,
     mbar_peri,
     smbar_calc,
     smbar_peri,
     mbar_julia,
     {"Mandelbar", "Mandelbar Julia"},
     "mbar",
     {0.0, 0.0, 2.5, 3.5},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, 0, 2, sym6},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete (from fractint) *//* 17 */
     FORMULAMAGIC,
     mlambda_calc,
     mlambda_peri,
     smlambda_calc,
     smlambda_peri,
     NULL,
     {"Lambda Mandelbrot", "Lambda"},
     "mlambda",
     {0.0, 0.0, 2.5, 5.5},
     0, 0, 0.5, 0.0,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete (from fractint) *//* 18 */
     FORMULAMAGIC,
     manowar_calc,
     NULL,
     smanowar_calc,
     NULL,
     NULL,
     {"Manowar", "Manowar Julia"},
     "manowar",
     {0.0, 0.0, 2.5, 2.5},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete (from fractint) *//* 19 */
     FORMULAMAGIC,
     spider_calc,
     NULL,
     sspider_calc,
     NULL,
     NULL,
     {"Spider", "Spider Julia"},
     "spider",
     {0.0, 0.0, 2.5, 4.5},
     1, 1, 0.0, 0.0,
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, 0, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete, method from fractint *//* 20 */
     FORMULAMAGIC,
     sier_calc,
     NULL,
     NULL,
     NULL,
     NULL,
     {"Sierpinski", "Sierpinski"},
     "sier",
     {0.5, 0.43, 1.5, 1.0},
     0, 0, 0.5, 0.8660254,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete, method from fractint *//* 21 */
     FORMULAMAGIC,
     carpet_calc,
     NULL,
     NULL,
     NULL,
     NULL,
     {"Sierpinski Carpet", "Sierpinski Carpet"},
     "carpet",
     {0.5, 0.5, 1.5, 1.5},
     0, 0, 1, 1,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete, method from fractint *//* 22 */
     FORMULAMAGIC,
     koch_calc,
     NULL,
     NULL,
     NULL,
     NULL,
     {"Koch Snowflake", "Koch Snowflake"},
     "koch",
     {0.0, 0.0, 2.5, 2.5},
     0, 1, 0, 0,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Z. Kovacs *//* 23 */
     FORMULAMAGIC,
     hornflake_calc,
     NULL,
     NULL,
     NULL,
     NULL,
     {"Spidron Hornflake", "Spidron Hornflake"},
     "hornflake",
     {-0.75, 0, 3.8756, 3.8756},
     0, 1, 0, 0,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     0,
     },
    {                           /* formula added by Z. Kovacs, originally mand6 but it was mand9 by accident *//* 24 */
     FORMULAMAGIC,
     mand9_calc,
     mand9_peri,
     smand9_calc,
     smand9_peri,
     mand9_julia,
     {"Mandelbrot^9", "Julia^9"},
     "mandel9",
     /*{1.25, -1.25, 1.25, -1.25}, */
     {0.0, 0.0, 2.5, 2.5},
     1, 1, 0.0, 0.0,
     {
      {0, 0, 6, sym16},
      {INT_MAX, 0, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, INT_MAX, 0, NULL},
      {0, 0, 0, NULL},
      {0, 0, 6, sym16},
      {INT_MAX, INT_MAX, 0, NULL},
      {0, 0, 6, sym16},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {0, 0, 6, sym16},
      {0, 0, 6, sym16},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },

    {                           /* formula added by S. Bizien *//* 25 */
     FORMULAMAGIC,
     beryl_calc,
     beryl_peri,
     NULL,
     NULL,
     NULL,
     {"Beryl", "Beryl"},
     "beryl",
     {-0.6, 0, 2, 2},
     0, 0, 1.0, 0.0,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete *//* 26 */
     FORMULAMAGIC,
     goldsier_calc,
     NULL,
     NULL,
     NULL,
     NULL,
     {"Golden Sierpinski", "Golden Sierpinski"},
     "goldsier",
     {0.5, 0.43, 1.5, 1.0},
     0, 0, 0.5, 0.8660254,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete *//* 27 */
     FORMULAMAGIC,
     circle7_calc,
     NULL,
     NULL,
     NULL,
     NULL,
     {"Circle 7", "Circle 7"},
     "circle7",
     {0.0, 0.0, 2.5, 2.5},
     0, 0, 0.0, 0.0,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {
     FORMULAMAGIC,
     clock_calc,
     NULL,
     NULL,
     NULL,
     NULL,
     {"Clock", "Clock"},
     "clock",
     {0.0, 0.0, 2.5, 2.5},
     0, 0, 0.0, 0.0,
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     },
    {                           /* formula added by Arpad Fekete *//* 28 */
     FORMULAMAGIC,
     symbarn_calc,
     NULL,
     NULL,
     NULL,
     symbarn_julia,
     {"Symmetric Barnsley M.", "Symmetric Barnsley"},
     "symbarn",
     {0.0, 0.0, 8.0, 1.0},
     0, 0, 1.3, 1.3,
     /* Arpad hasn't created the symmetry properties, */
     /* because he doesn't considered it to be important */
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     {
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      {INT_MAX, INT_MAX, 0, NULL},
      },
     MANDEL_BTRACE,
     }

#ifdef USE_SFFE
    , {                         /* formula added by M. Malczak - SFFE *//* 29 */
       FORMULAMAGIC,
       sffe_calc,
       NULL,
       ssffe_calc,
       NULL,
#endif
       sffe_julia,
       {"User defined", "User defined"},
       "user",
       /*{0.5, -2.0, 1.25, -1.25}, */
       /*{-0.75, 0.0, 1, 1}, */
       /* 2009-08-01 JB Langston
        * Changed default zoom level to match Mandelbrot
        */
       //{-0.75, 0.0, 2.5, 2.5},
       //0, 1, 0.0, 0.0,
       // Changed default view to match Burning Ship
       {0.0, 0.0, 5, 5},
       0, 1, 0.0, 0.0,
       {
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        },
       {
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        {INT_MAX, INT_MAX, 0, NULL},
        },
       MANDEL_BTRACE | SFFE_FRACTAL,
       }
};

const struct formula *currentformula;
const int nformulas = sizeof(formulas) / sizeof(struct formula);
const int nmformulas = 16; // Is this correct here? -- Zoltan, 2009-07-30
