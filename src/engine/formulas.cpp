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
#include "fbm_noise.h"
#include "randsctile.h"
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
                                   "fbm",
                                   "fbm + zmag",
                                   "fbm + decomposition",
                                   "cos(|z|)",
                                   "hexagons",
                                   "spiral",
                                   "cos(|z-c|)",
                                   "tessellation",
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
                                    "fbm + smooth",
                                    "fbm + iter",
                                    "fbm",
                                    "tessellation + smooth",
                                    "tessellation + iter",
                                    "tessellation",
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
    QT_TRANSLATE_NOOP("Menu", "Disable truecolor colouring"),
    QT_TRANSLATE_NOOP("Menu", "simple red (for education purposes)"),
    QT_TRANSLATE_NOOP("Menu", "simple blue (for education purposes)"),
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
    /* on the end, so that every number a saved position holds still names
     * what it named */
    "log log",
    "exp",
    "triangle",
    NULL
};

/* The speed change after ITER, by the number fractal_context holds: off, and
 * the two ways of joining a stretch to the one before it. See color_change. */
const char *const colorfun_modes[] = {"off", "jump", "smooth", NULL};

/* Whose angle an inside decomposition reads, by the number indecompof holds.
 * Formulas, and so the same in every language, as the mode names are. */
const char *const decomp_of[] = {"z", "z-c", "z*c", NULL};

/* Which way a colour decomposition walks the palette round a turn, by the
 * number decompway holds; see indecompway. */
const char *const decomp_ways[] = {QT_TRANSLATE_NOOP("Dialog", "right"),
                                   QT_TRANSLATE_NOOP("Dialog", "left"),
                                   QT_TRANSLATE_NOOP("Dialog", "both"), NULL};

/* How the outside's fbm marks are sized, by the number outfbmscaling
 * holds: the names the dialog offers and a position saves -- one word each,
 * since a keyword is, and "smooth" for the smooth count. */
const char *const fbm_scalings[] = {QT_TRANSLATE_NOOP("Dialog", "fixed"),
                                    QT_TRANSLATE_NOOP("Dialog", "iteration"),
                                    QT_TRANSLATE_NOOP("Dialog", "smooth"),
                                    NULL};

/* What the inside's tessellation mode colours by and how its tiles stand
 * apart, by the numbers intessellcolor and intessellfill hold: the names the
 * dialog offers and a position saves. */
const char *const tessell_colors[] = {QT_TRANSLATE_NOOP("Dialog", "angle"),
                                      QT_TRANSLATE_NOOP("Dialog", "tile-angle"),
                                      QT_TRANSLATE_NOOP("Dialog", "zmag"),
                                      QT_TRANSLATE_NOOP("Dialog", "none"),
                                      NULL};
const char *const tessell_fills[] = {QT_TRANSLATE_NOOP("Dialog", "turned"),
                                     QT_TRANSLATE_NOOP("Dialog", "shaded"),
                                     QT_TRANSLATE_NOOP("Dialog", "edged"),
                                     QT_TRANSLATE_NOOP("Dialog", "flat"), NULL};
/* And its tilings, the first numbered nought: each as the Values tab of the
 * formula reference names it, after the number randsctile and a position
 * know it by. Only the dialog shows them; the tiling is saved as its number,
 * by intessellation. */
const char *const tessell_tilings[] = {
    QT_TRANSLATE_NOOP("Dialog", "1: squares"),
    QT_TRANSLATE_NOOP("Dialog", "2: equilateral triangles"),
    QT_TRANSLATE_NOOP("Dialog", "3: hexagons"),
    QT_TRANSLATE_NOOP("Dialog", "4: octagons and squares (4.8.8)"),
    QT_TRANSLATE_NOOP("Dialog", "5: hexagons and triangles, the trihexagonal (3.6.3.6)"),
    QT_TRANSLATE_NOOP("Dialog", "6: hexagons, squares and triangles (3.4.6.4)"),
    QT_TRANSLATE_NOOP("Dialog", "7: dodecagons and triangles (3.12.12)"),
    QT_TRANSLATE_NOOP("Dialog", "8: dodecagons, hexagons and squares (4.6.12)"),
    QT_TRANSLATE_NOOP("Dialog", "9: rows of squares and of triangles (3.3.3.4.4)"),
    QT_TRANSLATE_NOOP("Dialog", "10: the snub square (3.3.4.3.4)"),
    QT_TRANSLATE_NOOP("Dialog", "11: the snub hexagonal, hexagons in a sea of triangles (3.3.3.3.6)"),
    QT_TRANSLATE_NOOP("Dialog", "12: tetrakis square: right triangles, eight round a point"),
    QT_TRANSLATE_NOOP("Dialog", "13: rhombille: rhombi that read as cubes"),
    QT_TRANSLATE_NOOP("Dialog", "14: deltoidal trihexagonal: kites"),
    QT_TRANSLATE_NOOP("Dialog", "15: triakis triangular: obtuse triangles"),
    QT_TRANSLATE_NOOP("Dialog", "16: kisrhombille: half-equilateral triangles, twelve round a point"),
    QT_TRANSLATE_NOOP("Dialog", "17: prismatic pentagons"),
    QT_TRANSLATE_NOOP("Dialog", "18: Cairo pentagons"),
    QT_TRANSLATE_NOOP("Dialog", "19: floret pentagons, six to a flower"),
    QT_TRANSLATE_NOOP("Dialog", "20: bricks, running bond"),
    QT_TRANSLATE_NOOP("Dialog", "21: Flemish bond: a brick and a half brick by turns"),
    QT_TRANSLATE_NOOP("Dialog", "22: herringbone"),
    QT_TRANSLATE_NOOP("Dialog", "23: basketweave"),
    QT_TRANSLATE_NOOP("Dialog", "24: Pythagorean: squares of two sizes"),
    QT_TRANSLATE_NOOP("Dialog", "25: chevrons"),
    QT_TRANSLATE_NOOP("Dialog", "26: squares and rhombi in alternate rows"),
    QT_TRANSLATE_NOOP("Dialog", "27: houses: a pentagon, up and hanging by turns"),
    QT_TRANSLATE_NOOP("Dialog", "28: rows: two of squares, one of triangles"),
    QT_TRANSLATE_NOOP("Dialog", "29: rows: one of squares, two of triangles"),
    QT_TRANSLATE_NOOP("Dialog", "30: hexagons among triangles"),
    QT_TRANSLATE_NOOP("Dialog", "31: Greek crosses"),
    QT_TRANSLATE_NOOP("Dialog", "32: T tetrominoes"),
    QT_TRANSLATE_NOOP("Dialog", "33: eight-pointed stars and crosses"),
    QT_TRANSLATE_NOOP("Dialog", "34: six-pointed stars and hexagons"),
    QT_TRANSLATE_NOOP("Dialog", "35: eight-pointed stars with four-pointed ones and rhombi"),
    QT_TRANSLATE_NOOP("Dialog", "36: twelve-pointed stars with three-pointed ones and rhombi"),
    QT_TRANSLATE_NOOP("Dialog", "37: twelve-pointed stars with notched hexagons and squares"),
    QT_TRANSLATE_NOOP("Dialog", "38: Voronoi cells, irregular -- randscp's"),
    QT_TRANSLATE_NOOP("Dialog", "39: Penrose rhombs: never repeats, fivefold"),
    QT_TRANSLATE_NOOP("Dialog", "40: Penrose kites and darts"),
    QT_TRANSLATE_NOOP("Dialog", "41: Ammann-Beenker: squares and rhombi, eightfold, never repeats"),
    QT_TRANSLATE_NOOP("Dialog", "42: rhombi in twelve directions, never repeats"),
    QT_TRANSLATE_NOOP("Dialog", "43: rhombi in seven directions, never repeats"),
    QT_TRANSLATE_NOOP("Dialog", "44: the pinwheel: triangles facing every way there is"),
    QT_TRANSLATE_NOOP("Dialog", "45: the chair: Ls cut into Ls"),
    NULL};

/* And the four 1.8.0 had, which a position from then names: read, and turned
 * into the two by uih_setcolorfunlegacy, but no longer offered. */
const char *const colorfun_after_modes[] = {"off", "jump", "continuous",
                                            "smooth", NULL};

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

/* --- the smooth count --------------------------------------------------------
 *
 * How far into its last iteration a point escaped, which smooth colouring adds
 * to the count of whole ones; in 256ths of an iteration.
 *
 * At the bailout -- smooth_count_near, what XaoS always did -- it is where the
 * threshold falls between the quantity before the escaping pass and after it,
 * linearly in its logarithm. The count runs on unbroken from one iteration to
 * the next, but not its pace: the quantity is all but squared from one pass to
 * the next, so the count's pace changes across every iteration and turns at
 * once where the next begins -- by four times, at the edge of an iteration in
 * a view of SPEED_PROBLEM_4. At a colouring speed of one that moves a colour
 * by a few hundredths of an entry; ten thousand times faster it is the bands
 * changing density at a stroke along every edge of an iteration.
 *
 * Past it -- smooth_count_far -- the orbit is carried on until the quantity
 * passes SMOOTH_FAR, where the squaring is all there is to it, and the
 * fraction is taken in the logarithm of the logarithm, where a squaring is a
 * step of one length:
 *
 *     done + ln(ln T / ln s) / ln(ln m / ln s)
 *
 * with s and m the quantity either side of SMOOTH_FAR, T the bailout and done
 * the passes made. The passes past the bailout take themselves back: each adds
 * one to done and doubles ln m, which takes one from the fraction. The degree
 * is the orbit's own, ln m / ln s, so a cube is counted as a cube. The count
 * is the one at the bailout to within a tenth of an iteration, and its pace
 * runs on unbroken from one iteration to the next: the turn of four at the
 * edge of that view comes out at 1.0000, and along a ray from outside the set
 * the turns of 1.3 to 1.8 at every edge at 1.0025 -- 1.0003 at 10^12, which
 * costs a pass more for nothing anyone could see.
 *
 * With the circle only. Every other shape measures the orbit by something that
 * turns with it while the passes double its angle, and carried on, the degree
 * came out unsteady enough to make the turns worse, and the square's even turn
 * back. Those, and the formulas that converge or test their escape their own
 * way, are counted at the bailout.
 *
 * And with a power only. The Barnsleys multiply the orbit by c, so that its
 * logarithm grows by a step of its own rather than doubling; the count at the
 * bailout, linear in that logarithm, is already the right one for them, and in
 * the logarithm of the logarithm their steps shrink to nothing and the count
 * runs off. They stay at the bailout, and so does an orbit of a user formula
 * that is well out -- past SMOOTH_POWER -- and does not grow as a power of one
 * and a half at least, which is given up on there. */
static const number_t SMOOTH_FAR = 1e6;
static const number_t SMOOTH_POWER = 1e3;

/* At the bailout: done passes, the last from s to m either side of t. */
static inline number_t smooth_count_near(unsigned int done, number_t s,
                                         number_t m, number_t t)
{
    /* Enough to keep a logarithm away from zero and no more. It was a
     * millionth, which is nothing beside a bailout of four and everything
     * beside a Newton convergence of a millionth; below one it now follows the
     * threshold down. At the bailouts anyone uses it is the millionth it
     * always was. */
    const number_t nudge =
        t < 1 ? (number_t)t * (number_t)1e-6 : (number_t)0.000001;
    m += nudge;
    s += nudge;
    return (number_t)(done * 256) +
           nlog((number_t)(t / s)) / nlog((number_t)(m / s)) * 256;
}

/* Past it: done passes, the last from s to m either side of SMOOTH_FAR, the
 * bailout t. Says whether there was a count to be had: not where the
 * logarithms have nothing to work on, a quantity that never went past one or
 * one that ran to infinity, nor where the last pass was no power of one and a
 * half. */
static inline int smooth_count_far(unsigned int done, number_t s, number_t m,
                                   number_t t, number_t &smoothed)
{
    const number_t ls = nlog(s), lm = nlog(m);
    if (!(ls > 0) || !(lm >= ls * (number_t)1.5) || !(t > 1) ||
        !(lm - lm == 0))
        return 0;
    /* ln ln t, worked out once a picture by set_fractalc */
    const number_t llt = t == cfractalc.bailoutloglogof ? cfractalc.bailoutloglog
                                                        : nlog(nlog(t));
    const number_t lls = nlog(ls);
    smoothed = ((number_t)done + (llt - lls) / (nlog(lm) - lls)) * 256;
    return 1;
}

#define SMOOTHOUTPUT()                                                         \
    {                                                                          \
        /* The smooth count, in 256ths of an iteration as the others are, and  \
         * whole, in number_t as the rest of the colouring is: carried past    \
         * the bailout where the orbit can be (SMOOTHBEYOND, see docalc.h),    \
         * at the bailout where not. It was made an int here, before the       \
         * colouring speed, the function and the speed change had it: a        \
         * 256th of an iteration is nothing at a speed of one, and a stripe    \
         * of one colour wherever they run the count faster -- ten entries a   \
         * step at 2500 an iteration, where the bands should have been finer   \
         * than the pixels. */                                                 \
        number_t smoothed = 0;                                                 \
        int beyond = 0;                                                        \
        SMOOTHBEYOND(beyond, smoothed);                                        \
        if (!beyond) {                                                         \
            PRESMOOTH;                                                         \
            smoothed = smooth_count_near(cfractalc.maxiter - iter, szmag,      \
                                         zre, SMOOTHTHRESHOLD);                \
        }                                                                      \
        /* fbm + iter and fbm, here only because their marks are sized by    \
         * the smooth count: that handed on, and the rest drawn by            \
         * color_output as ever. They read neither zre nor zim, which the     \
         * passes past the bailout have moved. The tessellation modes, here   \
         * for the smooth count as well, read them only for the angle of z,   \
         * which is then the angle a few passes on. */                        \
        if (cfractalc.coloringmode == OutColormodeClass::ColOut_fbm_iter ||    \
            cfractalc.coloringmode == OutColormodeClass::ColOut_fbm ||         \
            cfractalc.coloringmode.ColorMode >=                                \
                OutColormodeClass::ColOut_tessell_smooth) {                    \
            color_smooth = smoothed;                                           \
            color_smooth_held = 1;                                             \
            return (color_output(zre, zim, cfractalc.maxiter - iter));         \
        }                                                                      \
        if (cfractalc.coloringmode == OutColormodeClass::ColOut_smooth_log) {  \
            smoothed = nlog(smoothed) * ((cpalette.size - 1)) /                \
                           nlog((number_t)(cfractalc.maxiter * 256)) +         \
                       1;                                                      \
        }                                                                      \
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
            /* The smooth count worn.                                          \
             *                                                                 \
             * Here, and not in color_output: this macro works the count out    \
             * and returns a pixel of its own without ever passing through      \
             * color_output, so a case added there is never reached by a mode   \
             * that uses the smooth iteration function.                         \
             *                                                                 \
             * And here rather than four lines up, where it was: iter is an     \
             * unsigned int and the motion moves the count both ways, so a      \
             * negative addition wrapped it to four thousand million and the    \
             * cast below overflowed. What that drew was flat patches wherever  \
             * the motion pulled the count under nought. smoothed is the same   \
             * count in a number_t, which is what the negative case below is    \
             * already written for. */                                          \
            if (cfractalc.coloringmode ==                                      \
                OutColormodeClass::ColOut_fbm_smooth)                          \
                smoothed = fbm_fold(smoothed +                                 \
                                    fbm_at_pixel(0, fbm_count(                 \
                                        cfractalc.maxiter - iter,              \
                                        smoothed)) * 256);                     \
            int cycle = (int)(((unsigned int)(cpalette.size - 1)) << 8);       \
            color_precalc(smoothed, 0);                                        \
            iter = color_wrap(smoothed, (unsigned int)cycle);                  \
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


/* Where the pixel being coloured stands; see fractal.h. Written once a pixel
 * by calculate(), read by the fbm modes. */
thread_local number_t color_px = 0, color_py = 0;

/* The smooth count, in 256ths, that SMOOTHOUTPUT hands color_output for
 * fbm + iter and fbm when their marks are sized by it; held says it is there
 * for the pixel being coloured, and color_output takes it down again. */
static thread_local number_t color_smooth = 0;
static thread_local int color_smooth_held = 0;

/* --- a fractional Brownian motion, for the colouring modes that want the
 * picture to look used rather than clean ----------------------------------
 *
 * The motion itself is fbm_noise, which fbm() in the parser draws as well;
 * what it is and why it is gradient noise is written there.
 *
 * The roughness is the share of height each octave keeps of the one before.
 * At a half -- the plain motion -- the eighth octave carries a two hundred and
 * fiftieth of the whole, so octaves past four or five change nothing one can
 * see; measured, four and eighteen came out the same picture. Raising it gives
 * the fine octaves something to spend and the count begins to matter, which is
 * why both are offered.
 *
 * The seed is a number the user sets and nothing else: no clock, no pass, no
 * order the pixels happened to be computed in, so a picture comes back the
 * same tomorrow.
 *
 * Read at the pixel and not at the escape point. The escape point is not a
 * continuous function of the pixel -- it jumps wherever the escape iteration
 * does -- so noise read there comes out as grain, which was tried and thrown
 * away.
 */

/* The motion at the pixel, in bands of colour: nought to the intensity asked
 * for, and never below.
 *
 * It used to be centred, moving the value both ways and leaving it where it was
 * on average, and that is what took it under nought wherever the value it was
 * added to was smaller than the motion. What came of that had to be repaired --
 * wrapped round the palette, which drew a hard edge, or mirrored, which creased
 * the field -- and neither repair is as good as not needing one.
 *
 * A motion that only ever adds cannot take a value under nought, and it draws
 * the same picture: adding a half of the intensity to every point of a centred
 * motion is a constant, and a constant moves the whole picture along the
 * palette without changing a mark of it. What it costs is that shift, which
 * Coloring shift undoes for anyone who wants the old place back. */
/* A value under nought, mirrored.
 *
 * The motion cannot put one there any more -- it only ever adds -- so this is a
 * guard and not the mechanism: what it catches is a quantity that was negative
 * before the motion was added to it, which the smooth count can be when its
 * logarithm goes the other way. The colouring would wrap such a value round to
 * the far end of the palette and draw a hard edge along the crossing; a mirror
 * meets itself there. */
static inline number_t fbm_fold(number_t v) { return v < 0 ? -v : v; }

/* The count the outside's marks are sized by, in iterations: the whole one
 * or the smooth one, given in 256ths, as outfbmscaling says. */
static inline number_t fbm_count(unsigned int whole, number_t smooth)
{
    return cfractalc.outfbmscaling == FBM_SCALE_SMOOTH ? smooth / 256
                                                       : (number_t)whole;
}

/* The motion at the pixel, its marks sized by count when the outside asks
 * for it.
 *
 * Fixed, or at a factor of nought, it is the motion as it always was, and
 * the code that drew it runs as it did. Otherwise the frequency is doubled
 * factor times every ten iterations of the count. Laid over the pixel's
 * place as it stands, a frequency that changes from one pixel to the next
 * stretches the field -- the further from the origin, the more -- so it is
 * not: the motion is read at the two powers of two either side, each a
 * field of its own at a fixed frequency, and the two are blended by where
 * between them the count stands, as a texture is between its mipmaps. Two
 * fields blended half and half swing only seven tenths as far as either, so
 * the blend is taken back to the swing of one, and held within nought to
 * one, which it leaves only at the rarest extremes.
 *
 * With the whole count the place between is the same for a whole band, and
 * the marks change size at its edges; with the smooth count they change as
 * it does, without a seam. */
static number_t fbm_at_pixel(int inset, number_t count = 0)
{
    number_t freq = inset ? cfractalc.infbmfrequency : cfractalc.outfbmfrequency;
    number_t rough = inset ? cfractalc.infbmroughness : cfractalc.outfbmroughness;
    number_t much = inset ? cfractalc.infbmintensity : cfractalc.outfbmintensity;
    int octaves = inset ? cfractalc.infbmoctaves : cfractalc.outfbmoctaves;
    int seed = inset ? cfractalc.infbmseed : cfractalc.outfbmseed;

    if (inset || cfractalc.outfbmscaling == FBM_SCALE_FIXED ||
        !(cfractalc.outfbmfactor != 0))
        return fbm_noise(color_px * freq, color_py * freq,
                         (uint64_t)(int64_t)seed, octaves, rough) *
               much;

    number_t doublings = cfractalc.outfbmfactor * count / 10;
    /* past a thousand doublings either way the marks are long gone below a
     * pixel or above the picture, and the power of two has to fit an int */
    if (!(doublings > -1000))
        doublings = -1000;
    if (doublings > 1000)
        doublings = 1000;
    number_t whole = nfloor(doublings);
    double between = (double)(doublings - whole);
    number_t f = nldexp(freq, (int)whole);
    number_t a = fbm_noise(color_px * f, color_py * f, (uint64_t)(int64_t)seed,
                           octaves, rough);
    if (between == 0)
        return a * much;
    f *= 2;
    number_t b = fbm_noise(color_px * f, color_py * f, (uint64_t)(int64_t)seed,
                           octaves, rough);
    double m = (1 - between) * (double)a + between * (double)b;
    m = 0.5 + (m - 0.5) / sqrt((1 - between) * (1 - between) +
                               between * between);
    m = m < 0 ? 0 : (m > 1 ? 1 : m);
    return (number_t)m * much;
}

/* 2021-02-09 MSUMMO calculate color functions */
const number_t V_MAX = DBL_MAX / SMUL;
/* Held within V_MAX, the largest value whose 256ths a double still holds;
 * says whether it had to be. */
static int color_f_normalize(number_t &v) {
    if (v < -V_MAX) {
        v = -V_MAX;
        return 1;
    } else if (v > V_MAX) {
        v = V_MAX;
        return 1;
    }
    return 0;
}

/* --- the speed change after ITER ------------------------------------------
 *
 * What the value handed to the speed function goes through first, so that the
 * function chosen -- a square root, a log, anything -- goes on as it was up to
 * an ITER and runs faster from there, or slower for an exponent below one.
 *
 * The count is cut at the ITERs of the groups that are on, and each stretch
 * has one exponent of its own: the initial one from the start to the lower
 * ITER, the lower group's from there to the higher ITER, and the higher
 * group's past that -- each the stretch's alone, not multiplied by the one
 * before. A group that is off is left out altogether: with one of them on, the
 * initial exponent runs up to its ITER and its own from there; with neither,
 * the initial exponent runs on for every count. The lower ITER goes first
 * whichever group holds it, and at the same ITER the second group replaces
 * the first.
 *
 * Each stretch has a speed besides, which its pace is multiplied by: the
 * exponent bends the stretch, the speed tilts it. The first stretch is
 * s0 v^p0, from the start, whether the groups are on or off. (A group at an
 * ITER of nought or below puts its own exponent and speed there instead, as a
 * single group's change always did.) Each stretch after it starts from the
 * count the one before it came to, f(I), and adds s I ((v/I)^e - 1), whose
 * pace is s e at I: an exponent of one goes on as the straight line
 * f(I) + s (v - I), whatever came before, and two squares again from its ITER.
 * So every stretch is C + K (v/J)^e -- the first with C = 0 and J = 1 -- and
 * that is all that is kept of it. Every speed at one, as they start, is the
 * change as it was before there were speeds.
 *
 * A group joins its stretch to the one before in one of two ways:
 *
 *   jump    at I, at once: the pace turns there, which shows as a ring where
 *           the bands change pace
 *   smooth  over a radius r on either side of I. From I - r to I + r the pace
 *           goes in a straight line from the one the count has at I - r to
 *           the one it has at I + r, half way at I, and a parabola about I is
 *           added to it that gives back what the line cut off the jump's
 *           corner -- so that from I + r on the count is the jump's again, to
 *           the digit, and the radius changes the colours inside it and
 *           nowhere else. Nothing moves before I - r either.
 *
 * The rounding used to blend the two stretches' paces as each would have gone
 * on: a stretch whose pace was still growing -- a square -- kept it growing
 * well into the rounding, so the pace hardly fell before I and was half way
 * down only half a radius after it; and what that gained moved every colour
 * past the rounding on -- by r^2/3 bands after a square.
 *
 * The parabola never takes the pace below the slower of its two ends: where
 * it would -- a slow pace turned very fast over a short radius, where the
 * straight line alone already adds up to more than the jump -- it is held
 * where the pace only just stays above it, and the count past the rounding
 * moved on by what that keeps back. Held anywhere lower, the pace dipped
 * towards nought just inside the rounding and rose again, a corner in all but
 * name. The first rounding stops short of nought, and two that would overlap
 * are cut down alike, so that each stays between its neighbours.
 *
 * 1.8.0 had a mode more, continuous, and a smooth of another kind; a position
 * from then is read into these two, near enough, by colorfun_from_legacy.
 * Its jump is the jump here, and with one group, an initial exponent of one
 * and every speed one color_power still works it, so that such a picture is
 * drawn to the bit as it was.
 *
 * Odd, as the functions are: a negative value is turned as its opposite is
 * and keeps its sign. Worked in number_t, as the formulas are: long double in
 * the ordinary build, __float128 in quad.
 */

/* 1.8.0's jump: untouched up to iter and iter (v/iter)^p past it, or with
 * iter at nought or below the plain power from the start. */
static number_t color_power(number_t v, number_t iter, number_t p)
{
    number_t a = v < 0 ? -v : v, r;
    if (!(a > 0))
        return 0;
    if (!(iter > 0))
        r = nexp(p * nlog(a));
    else if (a <= iter)
        r = a;
    else
        r = iter * nexp(p * nlog(a / iter));
    return v < 0 ? -r : r;
}

struct colorfun_group {
    number_t iter, after, speed, radius;
    int mode;
};

/* Stretch s at v, above nought, and its pace there. */
static inline number_t color_stretch(const colorfun_plan &plan, int s,
                                     number_t v)
{
    return plan.C[s] + plan.K[s] * nexp(plan.e[s] * nlog(v / plan.J[s]));
}
static inline number_t color_stretch_pace(const colorfun_plan &plan, int s,
                                          number_t v)
{
    return plan.K[s] * plan.e[s] / plan.J[s] *
           nexp((plan.e[s] - 1) * nlog(v / plan.J[s]));
}

/* The count inside the rounding of corner c: from where the rounding begins,
 * what the straight line of paces and the parabola add up to. The parabola's
 * part, r (t - t^3/3 + 2/3) for t from -1 to 1, is written in s = 1 + t as
 * r s^2 (3 - s) / 3, which has no difference of near equals in it: near the
 * start the other form lost what a large parabola left of the digits. */
static number_t color_rounded(const colorfun_plan &plan, int c,
                              number_t a)
{
    const number_t r = (plan.hi[c] - plan.lo[c]) / 2, u = a - plan.lo[c];
    const number_t s = u / r;
    return plan.from[c] + plan.pacelo[c] * u +
           (plan.pacehi[c] - plan.pacelo[c]) * u * u / (4 * r) +
           plan.bump[c] * r * s * s * (3 - s) / 3;
}

static void color_change_plan(colorfun_plan &plan, number_t initial,
                              number_t initialspeed, colorfun_group g1,
                              colorfun_group g2)
{
    plan = colorfun_plan();
    colorfun_group g[2];
    int n = 0;
    if (g1.mode)
        g[n++] = g1;
    if (g2.mode) {
        if (n && g2.iter == g[0].iter)
            g[0] = g2;
        else
            g[n++] = g2;
    }
    if (n == 2 && g[1].iter < g[0].iter) {
        colorfun_group t = g[0];
        g[0] = g[1];
        g[1] = t;
    }
    if (n == 1 && initial == 1 && initialspeed == 1 && g[0].mode == 1 &&
        g[0].speed == 1) {
        /* at an exponent of one there is nothing to do, as there was not */
        if (g[0].after == 1)
            return;
        plan.kind = 1;
        plan.iter = g[0].iter;
        plan.exponent = g[0].after;
        return;
    }

    /* the groups at nought or below are all at the start, and only the last
     * of them says anything: the exponent and the speed the first stretch
     * has */
    int k = 0;
    while (k < n && !(g[k].iter > 0))
        k++;
    plan.C[0] = 0;
    plan.J[0] = 1;
    plan.e[0] = k == 0 ? initial : g[k - 1].after;
    plan.K[0] = k == 0 ? initialspeed : g[k - 1].speed;
    int identity = plan.e[0] == 1 && plan.K[0] == 1;
    number_t radius[2] = {0, 0}, speed[2] = {1, 1};
    for (int j = k; j < n; j++) {
        const int c = plan.corners++;
        plan.at[c] = g[j].iter;
        plan.e[c + 1] = g[j].after;
        speed[c] = g[j].speed;
        if (g[j].mode == 2 && g[j].radius > 0)
            radius[c] = g[j].radius;
        if (g[j].after != 1 || g[j].speed != 1)
            identity = 0;
    }
    /* every exponent and speed one: the count as it was, rounded or not */
    if (identity)
        return;
    plan.kind = 2;

    /* the first rounding stops short of nought; two that would overlap are
     * cut down alike */
    if (plan.corners && radius[0] > plan.at[0] * 0.999L)
        radius[0] = plan.at[0] * 0.999L;
    if (plan.corners == 2 && radius[0] + radius[1] > plan.at[1] - plan.at[0]) {
        const number_t cut =
            (plan.at[1] - plan.at[0]) / (radius[0] + radius[1]);
        radius[0] *= cut;
        radius[1] *= cut;
    }

    /* each stretch from the count the one before it came to; and each
     * rounding from the count and pace at its two ends */
    for (int c = 0; c < plan.corners; c++) {
        const number_t I = plan.at[c], r = radius[c];
        plan.C[c + 1] = color_stretch(plan, c, I) - speed[c] * I;
        plan.K[c + 1] = speed[c] * I;
        plan.J[c + 1] = I;
        plan.lo[c] = I - r;
        plan.hi[c] = I + r;
        if (!(r > 0))
            continue;
        plan.from[c] = color_stretch(plan, c, plan.lo[c]);
        plan.pacelo[c] = color_stretch_pace(plan, c, plan.lo[c]);
        plan.pacehi[c] = color_stretch_pace(plan, c + 1, plan.hi[c]);
        /* what the straight line leaves of the jump's count, given back by a
         * parabola, whose area over the rounding is 4/3 r of its height */
        const number_t lost = color_stretch(plan, c + 1, plan.hi[c]) -
                                 plan.from[c] -
                                 r * (plan.pacelo[c] + plan.pacehi[c]);
        plan.bump[c] = 3 * lost / (4 * r);
        /* the lowest height that keeps a + b t + h (1 - t^2) between its two
         * ends from t = -1 to 1, b being half the line's rise: the one that
         * puts the parabola's lowest point on an end, h = -|b|/2 */
        const number_t held =
            -nfabs(plan.pacehi[c] - plan.pacelo[c]) / 4;
        if (plan.bump[c] < held) {
            plan.C[c + 1] += 4 * r * (held - plan.bump[c]) / 3;
            plan.bump[c] = held;
        }
    }
}

void color_change_prepare(fractal_context *c)
{
    color_change_plan(
        c->inplan, (number_t)c->incolorfuninitial,
        (number_t)c->incolorfuninitialspeed,
        {(number_t)c->incolorfuniter, (number_t)c->incolorfunexp,
         (number_t)c->incolorfunspeed, (number_t)c->incolorfunradius,
         c->incolorfunmode},
        {(number_t)c->incolorfun2iter, (number_t)c->incolorfun2exp,
         (number_t)c->incolorfun2speed, (number_t)c->incolorfun2radius,
         c->incolorfun2mode});
    color_change_plan(
        c->outplan, (number_t)c->outcolorfuninitial,
        (number_t)c->outcolorfuninitialspeed,
        {(number_t)c->outcolorfuniter, (number_t)c->outcolorfunexp,
         (number_t)c->outcolorfunspeed, (number_t)c->outcolorfunradius,
         c->outcolorfunmode},
        {(number_t)c->outcolorfun2iter, (number_t)c->outcolorfun2exp,
         (number_t)c->outcolorfun2speed, (number_t)c->outcolorfun2radius,
         c->outcolorfun2mode});
}

/* A group as a position from before the radius holds it, in the four modes
 * of then -- off, jump, continuous, smooth -- read into the three there are,
 * near enough; the radius is left as it is unless the group is to be smooth.
 *
 * Off and jump are what they were. 1.8.0's smooth rounded the corner in the
 * logarithm of the count, a little on either side of ITER however far off;
 * the smooth here at a radius of three quarters of ITER keeps within four
 * hundredths of ITER of it for the exponents up to one and a half, and
 * within a tenth up to two. Continuous went on from ITER at the pace the
 * count came with, which no stretch here starts at: it is drawn as smooth at
 * three quarters of ITER, with the exponent that keeps nearest it,
 * 0.66 + p/4 + 0.09 p^2 -- one at one, where continuous changed nothing
 * either -- within twelve hundredths of ITER of it up to a half. Found by
 * setting the curves side by side over four ITERs. */
void colorfun_from_legacy(number_t iter, number_t *exponent, int *mode,
                          number_t *radius)
{
    switch (*mode) {
    case 1:
        break;
    case 2:
        *mode = 2;
        *exponent = (number_t)0.66 + *exponent / 4 +
                    (number_t)0.09 * *exponent * *exponent;
        if (iter > 0)
            *radius = iter * 3 / 4;
        break;
    case 3:
        *mode = 2;
        if (iter > 0)
            *radius = iter * 3 / 4;
        break;
    default:
        *mode = 0;
        break;
    }
}

/* The count v as the plan changes it. */
static number_t color_change(number_t v, const colorfun_plan &plan)
{
    if (plan.kind == 1)
        return color_power(v, plan.iter, plan.exponent);
    number_t a = v < 0 ? -v : v, r;
    if (!(a > 0))
        return 0;
    int s = 0;
    for (int c = 0; c < plan.corners; c++) {
        if (a > plan.lo[c] && a < plan.hi[c]) {
            r = color_rounded(plan, c, a);
            return v < 0 ? -r : r;
        }
        if (a >= plan.at[c])
            s = c + 1;
    }
    r = color_stretch(plan, s, a);
    return v < 0 ? -r : r;
}

/* v bent by speed function func, as colorfun[] names them; held is set when
 * a value had to be held within V_MAX on the way. Ident hands v back as it
 * came. */
static number_t color_f_bend(number_t v, const int func, int &held)
{
    switch (func) {
    case 1:
        // LOG
        if (v > 0)
            v = nlog(v + 1);
        else
            v = -nlog(-v + 1);
        break;
    case 2:
        // ATAN
        v = natan(v);
        break;
    case 3:
        // SIN
        v = nsin(v);
        break;
    case 4:
        // COS
        v = ncos(v);
        break;
    case 5:
        // SQR
        if (v > 0)
            v = v * v;
        else
            v = -(v * v);
        held |= color_f_normalize(v);
        break;
    case 6:
        // CUBE
        v = v * v * v;
        held |= color_f_normalize(v);
        break;
    case 7:
        // SQRT
        if (v > 0)
            v = nsqrt(v);
        else
            v = -nsqrt(-v);
        break;
    case 8:
        // CBRT
        if (v > 0)
            v = ncbrt(v);
        else
            v = -ncbrt(-v);
        break;
    case 9:
        /* LOG LOG: log once more, for the deep zooms where the count runs to
         * millions and log alone still packs the bands too close -- at a
         * million it gives 13.8 bands and this 2.7 */
        if (v > 0)
            v = nlog(1 + nlog(v + 1));
        else
            v = -nlog(1 + nlog(-v + 1));
        break;
    case 10:
        /* EXP: the other way from log, e^v - 1, which opens the small values
         * out. For the quantities that stay small -- the components that true
         * colour and the decompositions read; at a count of ninety it is
         * already past 10^38. */
        if (v > 0)
            v = nexpm1(v);
        else
            v = -nexpm1(-v);
        held |= color_f_normalize(v);
        break;
    case 11: {
        /* TRIANGLE: the sine with its curves straightened -- the same period
         * and height, the palette walked forwards and back at a steady pace
         * instead of lingering at the turns */
        number_t s = v * 2 / N_PI;
        number_t t = s + 1 - 4 * nfloor((s + 1) / 4);
        v = 1 - nfabs(t - 2);
        break;
    }
    }
    return v;
}

static number_t color_f_func(number_t v, const float speed, const int func,
                             const int shift, const int inset)
{
    int held = color_f_normalize(v);
    /* The speed change after ITER, on the value before the function has it.
     * Off -- the default -- or with every exponent at one, nothing at all is
     * done to it, so that the value goes on to the function to the bit as it
     * did before there was such a thing. */
    const colorfun_plan &plan = inset ? cfractalc.inplan : cfractalc.outplan;
    if (plan.kind) {
        v = color_change(v, plan);
        held |= color_f_normalize(v);
    }
    v = color_f_bend(v, func, held);
    v = v * speed;
    v += shift;
    held |= color_f_normalize(v);
    /* Handed back whole. It was a float, which past 65536 entries kept less
     * than a 256th of one -- the gradient between two entries in steps --
     * and past sixteen million skipped whole entries. Held within V_MAX by
     * color_f_normalize, it goes on to color_wrap, which takes any value
     * round the palette.
     *
     * Except a value that had to be held there -- an infinity, or a number
     * past what a double's 256ths hold, coming in or on the way -- and is
     * still past what a float's 256ths hold: such a value, made 256ths in a
     * float, came out infinite and landed on the first entry, where a value
     * that is no number lands. It is made one, and lands there still. What
     * colour an infinity takes is anybody's guess, and this way it takes the
     * one it always took: the middle of a potential gone to infinity, the
     * ground round a formula that escapes at once. */
    if (held && (v > FLT_MAX / SMUL || v < -FLT_MAX / SMUL))
        v = (number_t)NAN;
    return v;
}

/* Where a colour value lands in a palette cycle of cycle entries -- in 256ths
 * of an entry, as the colouring counts.
 *
 * It was converted to an int and wrapped, which is right as far as an int
 * reaches and undefined past it: 2^31 256ths is 8.4 million bands, which cube
 * reaches at 204 iterations, and every pixel past that took one and the same
 * entry, a flat patch where the bands should have gone on. Inside an int this
 * is that conversion and wrap exactly; outside it, the same wrap worked in
 * number_t, which holds the whole value. */
static unsigned int color_wrap(number_t v, unsigned int cycle)
{
    if (!(v == v) || cycle == 0)
        return 0;
    /* Rounded to a 256th of a 256th before it is cut to a 256th: the
     * colouring tells nothing finer, and a value that lies on the edge of a
     * 256th -- the real part of a Newton orbit come to rest on a root, at one
     * or at minus a half -- otherwise falls either side of it with the last
     * digits of the convergence. One way or the other it is a 256th of an
     * entry, which nobody sees, until an edge filter draws every such pixel
     * as an edge and the contours come out toothed like gears. A float used
     * to round most of them away, and left notches. */
    v = nfloor(v * 256 + (number_t)0.5) / 256;
    if (v > -(number_t)2147483647 && v < (number_t)2147483647) {
        int i = (int)v;
        if (i < 0) {
            i = cycle - ((-i) % cycle) - 1;
            if (i < 0)
                i = 0;
        }
        return ((unsigned int)i) % cycle;
    }
    number_t r = nfmod(ntrunc(v < 0 ? -v : v), (number_t)cycle);
    if (v >= 0)
        return (unsigned int)r % cycle;
    return (cycle - (unsigned int)r - 1) % cycle;
}

static void color_precalc(number_t &iter, int inset) {
    const float speed = inset ? cfractalc.incolorspeed : cfractalc.outcolorspeed;
    const int func = inset ? cfractalc.incolorfun : cfractalc.outcolorfun;
    const int shift = inset ? cfractalc.incolorshift : cfractalc.outcolorshift;
    iter = (color_f_func(iter / SMUL, speed, func, shift, inset) * SMUL);
}
static void color_precalc(number_t &zre, number_t &zim, number_t &pre,
                          number_t &pim, int inset) {
    const float speed = inset ? cfractalc.incolorspeed : cfractalc.outcolorspeed;
    const int func = inset ? cfractalc.incolorfun : cfractalc.outcolorfun;
    const int shift = inset ? cfractalc.incolorshift : cfractalc.outcolorshift;
    zre = color_f_func(zre, speed, func, shift, inset);
    zim = color_f_func(zim, speed, func, shift, inset);
    pre = color_f_func(pre, speed, func, shift, inset);
    pim = color_f_func(pim, speed, func, shift, inset);
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

/* --- the decompositions -----------------------------------------------------
 *
 * Where in a turn the angle of (re, im) lies, measured as the decompositions
 * always have -- atan2(re, im), from the imaginary axis, a quarter turn on --
 * and turned by rotation degrees: in turns, with any whole number of them
 * over. */
static number_t decomp_turn(number_t re, number_t im, number_t rotation)
{
    return natan2(re, im) / (N_PI + N_PI) + (number_t)0.75 + rotation / 360;
}

/* The walk round a turn: where in its repeat a place v in the turn, from
 * nought to one, lies, also from nought to one -- that many repeats of the
 * turn, each walked right, left, or right and back again. */
static number_t decomp_walk(number_t v, number_t repeats, int way)
{
    number_t x = v * repeats;
    number_t w = x - nfloor(x);
    if (way == 1)
        w = 1 - w;
    else if (way == 2)
        w = 1 - nfabs(2 * w - 1);
    return w;
}

/* Which of n sectors a place in a repeat, from nought to one, falls in. One
 * itself, where walking left starts, is in the last. */
static int decomp_index(number_t w, int n)
{
    int s = (int)(w * n);
    if (s < 0)
        s = 0;
    return s < n ? s : n - 1;
}

/* The parity of a band number, which may be past what an int holds. */
static int decomp_odd(number_t band)
{
    band = nfloor(band);
    if (!(band > -(number_t)1e18 && band < (number_t)1e18))
        return 0;
    return (int)((long long)band & 1);
}

/* What an odd band is moved along the palette by: half a cycle of it, in
 * 256ths, which sets it as far from its neighbours as the palette goes. */
static number_t decomp_half(void)
{
    return (number_t)(((unsigned int)(cpalette.size - 1)) << 8) / 2;
}

/* A colour decomposition's value, for a turn.
 *
 * Where in the turn it lies, from the seam the mode has always had, a quarter
 * turn on from where decomp_turn counts; walked round its repeats; the
 * sector of the repeat it falls in, if there are sectors; and spread over
 * speed times the palette a turn has always taken, 20000, from the quarter of
 * it on that the mode has always started at. At one, right, one, it is the
 * turn as it was. */
static number_t decomp_colour(number_t turn, int sectors, number_t repeats,
                              int way, number_t speed)
{
    number_t u = turn - (number_t)0.25;
    u -= nfloor(u);
    number_t w = decomp_walk(u, repeats, way);
    if (sectors > 0)
        w = (number_t)decomp_index(w, sectors) / sectors;
    return ((number_t)0.25 + speed * w) * 20000;
}

/* The binary decompositions colour alternate sectors, so the number is
 * even: two, the half-planes, where nothing is said. */
static int decomp_binary_sectors(int sectors)
{
    return sectors > 0 ? sectors + (sectors & 1) : 2;
}

/* Whether a turn falls in an odd sector of a binary decomposition: walked
 * round its repeats from where decomp_turn counts, which puts the two
 * half-planes on either side of the real axis as they always were. */
static int decomp_binary_odd(number_t turn, int sectors, number_t repeats,
                             int way)
{
    return decomp_index(decomp_walk(turn - nfloor(turn), repeats, way),
                        decomp_binary_sectors(sectors)) &
           1;
}

/* --- the tessellations ------------------------------------------------------
 *
 * What a point of a tessellation is coloured, inside or out, once its tile is
 * known: turn, the tile's own amount from nought to one, the hash that names
 * it; out, how far out of its middle the point stands, from nought to one;
 * and (rx, ry), where it stands seen from that middle. colour and fill are
 * what the colour comes from and how the tile stands apart, as
 * intessellcolor says; (are, aim) is what the angle of TESSELL_ANGLE is of,
 * zmag the value of TESSELL_ZMAG, and the rest the side's decomposition
 * settings, with ring_odd() saying whether the point is in an odd ring of
 * what the angle is of.
 *
 * much, the intensity, multiplies every amount the tile has of its own --
 * its turn, the bevel and the edge -- and nothing else: at one, as they
 * always were, to the bit, since a product by one is exact; at nought what
 * the colour comes from alone. */
template <class Ring>
static number_t tessell_colour(int colour, int fill, number_t much,
                               number_t turn, double out, double rx,
                               double ry, number_t are, number_t aim,
                               number_t zmag, number_t rotation, int sectors,
                               number_t repeats, int way, number_t speed,
                               number_t radial, Ring ring_odd)
{
    number_t i_f;
    /* A tile's own amount, along the palette: for an angle as much of it as
     * the angle's turn has always taken, so that moving a decomposition
     * setting keeps the tiles' colours where they were; for the rest, and
     * for a flat tile, the whole cycle. */
    const number_t cycle = 2 * decomp_half();
    if (fill == TESSELL_FLAT)
        i_f = much * (turn * cycle);
    else if (colour == TESSELL_ANGLE || colour == TESSELL_MIDDLE) {
        /* the angle, of z or about the tile's middle, cut and walked as the
         * decompositions are; the rings, of the modulus for z and of how far
         * out of the middle the point stands for the tile, which makes them
         * the tile's shape shrunk */
        const int middle = colour == TESSELL_MIDDLE;
        i_f = decomp_colour(decomp_turn(middle ? (number_t)rx : are,
                                        middle ? (number_t)ry : aim, rotation),
                            sectors, repeats, way, speed) +
              much * (turn * 20000);
        if (middle ? radial > 0 && decomp_odd(radial * (number_t)out)
                   : ring_odd())
            i_f += decomp_half();
    } else if (colour == TESSELL_ZMAG)
        i_f = zmag + much * (turn * cycle);
    else
        i_f = much * (turn * cycle);
    /* shaded: an eighth of what a turn of the angle takes, from the middle to
     * the edge -- a bevel and not a rainbow, as the modes' own spans are,
     * whatever the size of the palette; edged: the outer part of the tile in
     * the other half of the palette */
    if (fill == TESSELL_SHADED)
        i_f += much * ((number_t)out * 2500);
    else if (fill == TESSELL_EDGED && out > TESSELL_EDGE)
        i_f += much * decomp_half();
    return i_f;
}

/* The tile of the outside's tiling under (x, y), and where in it the point
 * stands when something reads that; turn is nought where no tile is to be
 * had, past the reach of a tiling that never repeats. */
static number_t tessell_turn(int tiling, int colour, int fill, number_t x,
                             number_t y, double &out, double &rx, double &ry)
{
    number_t turn = 0;
#ifdef USE_SFFE
    uint64_t name;
    if (colour == TESSELL_MIDDLE || fill == TESSELL_SHADED ||
                fill == TESSELL_EDGED
            ? randsctile_place(tiling, x, y, &name, &out, &rx, &ry)
            : randsctile_name(tiling, x, y, &name))
        turn = (number_t)(name >> 11) / (number_t)((uint64_t)1 << 53);
#else
    (void)tiling, (void)colour, (void)fill, (void)x, (void)y, (void)out,
        (void)rx, (void)ry;
#endif
    return turn;
}

/* A value of the outside, through the colouring speed and function, onto
 * the palette: what color_output ends with, for the tessellation modes,
 * which may look up two values for one pixel. */
static unsigned int color_lookup(number_t i_f)
{
    color_precalc(i_f, 0);
    unsigned int iter = color_wrap(i_f, ((unsigned int)(cpalette.size - 1))
                                            << 8);
    if ((cpalette.type & (C256 | SMALLITER)) || !(iter & 255))
        return (cpalette.pixels[1 + (iter >> 8)]);
    unsigned int i1 = cpalette.pixels[1 + (iter >> 8)], i2;
    if ((int)(iter >> 8) == cpalette.size - 2)
        i2 = cpalette.pixels[1];
    else
        i2 = cpalette.pixels[2 + (iter >> 8)];
    return (interpoltype(cpalette, i2, i1, iter & 255));
}

/* The outside's tessellation modes: the tiles of one of randsctile's tilings
 * laid over the pixel's place in the plane, frequency to a unit -- not over
 * the escape point, which jumps from band to band and would give grain, as
 * it did for the fbm motion -- worn by the smooth count, by the whole one,
 * or alone. iter is the whole count in 256ths, smooth the smooth one.
 *
 * Sized by the count, the frequency doubles factor times every ten
 * iterations, as the fbm marks do, and for the same reason it is not laid
 * over the plane as it stands, which would stretch the tiles: the tiles are
 * read at the two powers of two either side of where the count stands, each
 * a tiling of its own at a fixed size, and the two pixels are blended by
 * where between them the count stands. Tiles have no middle value to blend
 * towards, as the motion has, so it is the colours that are blended, a
 * dissolve from one size to the next; a palette of 256 entries, which has
 * no colours between its own, takes the nearer. */
static unsigned int tessell_output(number_t zre, number_t zim,
                                   unsigned int iter, number_t smooth)
{
    const int colour = cfractalc.outtessellcolor;
    const int fill = cfractalc.outtessellfill;
    const number_t base =
        cfractalc.coloringmode == OutColormodeType::ColOut_tessell_smooth
            ? smooth
            : cfractalc.coloringmode == OutColormodeType::ColOut_tessell_iter
                  ? (number_t)iter
                  : 0;
    /* zmag has no use outside, where the modulus is past the bailout: the
     * log of it, as iter+log(mag) adds, is what tells the points apart */
    const number_t zmag =
        colour == TESSELL_ZMAG ? nlog(zre * zre + zim * zim + 1) * SMUL * 4
                               : 0;
    /* the bands of the decomposition settings' iterations per band, as the
     * colour decomposition has them */
    auto ring_odd = [&]() {
        return cfractalc.outdecompradial > 0 &&
               ((iter >> SHIFT) / cfractalc.outdecompradial) & 1;
    };
    auto at = [&](number_t f) {
        double out = 0, rx = 0, ry = 0;
        const number_t turn =
            tessell_turn(cfractalc.outtessellation, colour, fill,
                         color_px * f, color_py * f, out, rx, ry);
        return color_lookup(
            base + tessell_colour(colour, fill, cfractalc.outtessellintensity,
                                  turn, out, rx, ry, zre, zim, zmag,
                                  cfractalc.outdecomprotation,
                                  cfractalc.outdecompsectors,
                                  cfractalc.outdecomprepeats,
                                  cfractalc.outdecompway,
                                  cfractalc.outdecompspeed,
                                  (number_t)cfractalc.outdecompradial,
                                  ring_odd));
    };
    number_t f = cfractalc.outtessellfrequency;
    if (cfractalc.outtessellscaling == FBM_SCALE_FIXED ||
        !(cfractalc.outtessellfactor != 0))
        return at(f);
    const number_t count = cfractalc.outtessellscaling == FBM_SCALE_SMOOTH
                               ? smooth / 256
                               : (number_t)(iter >> SHIFT);
    number_t doublings = cfractalc.outtessellfactor * count / 10;
    /* past a thousand doublings either way the tiles are long gone below a
     * pixel or above the picture, and the power of two has to fit an int */
    if (!(doublings > -1000))
        doublings = -1000;
    if (doublings > 1000)
        doublings = 1000;
    const number_t whole = nfloor(doublings);
    const double between = (double)(doublings - whole);
    f = nldexp(f, (int)whole);
    const unsigned int a = at(f);
    if (between == 0)
        return a;
    const unsigned int b = at(2 * f);
    if (cpalette.type & (C256 | SMALLITER))
        return between < 0.5 ? a : b;
    int w = (int)(between * 256);
    w = w < 0 ? 0 : (w > 255 ? 255 : w);
    return (interpoltype(cpalette, b, a, w));
}

static unsigned int color_output(number_t zre, number_t zim, unsigned int iter);
static unsigned int color_output(number_t zre, number_t zim, unsigned int iter)
{
    number_t i_f;
    /* The smooth count SMOOTHOUTPUT handed on, or, where none was worked
     * out, the whole one in its place. */
    const number_t smooth =
        color_smooth_held ? color_smooth : (number_t)iter * 256;
    color_smooth_held = 0;
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
            if (decomp_plain(&cfractalc, 0)) {
                if (zim > 0)
                    i_f = ((cfractalc.maxiter << SHIFT) - iter);
                break;
            }
            {
                /* sectors of the repeats, and bands of so many iterations
                 * across them; the odd ones moved from the count towards the
                 * count backwards by the speed -- all the way at one, as they
                 * always were, not at all at nought, twice as far at two */
                int odd = decomp_binary_odd(
                    decomp_turn(zre, zim, cfractalc.outdecomprotation),
                    cfractalc.outdecompsectors, cfractalc.outdecomprepeats,
                    cfractalc.outdecompway);
                if (cfractalc.outdecompradial > 0)
                    odd ^= ((iter >> SHIFT) / cfractalc.outdecompradial) & 1;
                if (odd) {
                    number_t back = (number_t)(cfractalc.maxiter << SHIFT) -
                                    (number_t)iter;
                    i_f = iter + cfractalc.outdecompspeed * (back - iter);
                }
            }
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
            if (decomp_plain(&cfractalc, 0)) {
                i_f = ((atan2l(zre, zim) / (M_PI + M_PI) + 0.75) *
                       20000);
                break;
            }
            i_f = decomp_colour(decomp_turn(zre, zim,
                                            cfractalc.outdecomprotation),
                                cfractalc.outdecompsectors,
                                cfractalc.outdecomprepeats,
                                cfractalc.outdecompway,
                                cfractalc.outdecompspeed);
            if (cfractalc.outdecompradial > 0 &&
                ((iter >> SHIFT) / cfractalc.outdecompradial) & 1)
                i_f += decomp_half();
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
        case OutColormodeType::ColOut_fbm_iter:
            /* the plain count worn, for a formula or a taste that does not
             * want the smoothing */
            i_f = fbm_fold(iter + fbm_at_pixel(0, fbm_count(iter >> SHIFT,
                                                            smooth)) *
                                      SMUL);
            break;
        case OutColormodeType::ColOut_fbm:
            /* the motion alone: no count at all, so what is drawn outside the
             * set is the noise itself -- though its marks may be sized by one */
            i_f = fbm_fold(fbm_at_pixel(0, fbm_count(iter >> SHIFT, smooth)) *
                           SMUL);
            break;
        case OutColormodeType::ColOut_tessell_smooth:
        case OutColormodeType::ColOut_tessell_iter:
        case OutColormodeType::ColOut_tessell:
            /* a pixel of their own, which may be two blended */
            return tessell_output(zre, zim, iter, smooth);
        case OutColormodeType::ColOut_abs_real_minus_abs_imag:
            /* how lopsided the escape point is, with no count: the level sets
             * of the difference rather than of the modulus */
            i_f = ((myabs(zre) - myabs(zim)) * SMUL * 200);
            break;
        default:
            break;
    }

    color_precalc(i_f, 0);
    iter = color_wrap(i_f, ((unsigned int)(cpalette.size - 1)) << 8);
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
    number_t i_f = iter;
    /* The argument a of every periodic function a mode uses -- a sine, a
     * cosine, frac, the parity of a stripe -- goes in as k*f(a): bent by the
     * speed function f the inside periodic function names, so that it comes
     * round faster in some places than in others, and multiplied by the
     * frequency k, how often it comes round. A periodic function inside the
     * argument of another has its own argument made k*f(a) as well. The
     * constant of a period written into a mode -- the eighth of frac(mag),
     * the fortieth of squares -- belongs to the function, and stays outside.
     * At ident and one, k*f(a) is a, to the bit. */
    const number_t k = cfractalc.incolorfrequency;
    const int g = cfractalc.incolorperiodicfun;
    auto periodic = [k, g](number_t a) {
        int held = 0;
        return k * (g ? color_f_bend(a, g, held) : a);
    };
    if (cfractalc.coloringmode.ColorMode) {
        color_f_normalize(zre);
        color_f_normalize(zim);
        color_f_normalize(pre);
        color_f_normalize(pim);
    }
    /* The decompositions with their settings moved: the angle of z, z-c or
     * z*c, turned and cut into sectors, and rings of its modulus across
     * them, through k*f. Only where the settings are not where they start,
     * which keeps the modes' own code below to the bit. */
    const int plain = decomp_plain(&cfractalc, 1);
    number_t qre = zre, qim = zim;
    if (cfractalc.indecompof == 1)
        qre = zre - pre, qim = zim - pim;
    else if (cfractalc.indecompof == 2)
        qre = zre * pre - zim * pim, qim = zre * pim + zim * pre;
    auto ring_odd = [&]() {
        return cfractalc.indecompradial > 0
                   ? decomp_odd(cfractalc.indecompradial *
                                periodic(nsqrt(qre * qre + qim * qim)))
                   : 0;
    };
    auto decomposition = [&]() {
        number_t v = decomp_colour(
            decomp_turn(qre, qim, cfractalc.indecomprotation),
            cfractalc.indecompsectors, cfractalc.indecomprepeats,
            cfractalc.indecompway, cfractalc.indecompspeed);
        return ring_odd() ? v + decomp_half() : v;
    };

    switch (cfractalc.incoloringmode) {
    case INCOLORING_FBM:
        /* the motion alone, which gives the inside a surface where it had one
         * flat tone */
        i_f = fbm_fold(fbm_at_pixel(1) * SMUL);
        break;
    case INCOLORING_FBM_ZMAG:
        /* zmag worn: the inside keeps what it says about the orbit and stops
         * being clean about it */
        i_f = fbm_fold(((zre * zre + zim * zim) *
                            (number_t)(cfractalc.maxiter >> 1) * SMUL +
                        SMUL) +
                       fbm_at_pixel(1) * SMUL);
        break;
    case INCOLORING_FBM_DECOMP:
        /* and the same over the decomposition */
        if (plain)
            i_f = fbm_fold(((atan2l(zre, zim) / (M_PI + M_PI) + 0.75) * 20000) +
                           fbm_at_pixel(1) * SMUL);
        else
            i_f = fbm_fold(decomposition() + fbm_at_pixel(1) * SMUL);
        break;
    case 1: /* zmag */
        i_f = (((zre * zre + zim * zim) *
                    (number_t)(cfractalc.maxiter >> 1) * SMUL +
                SMUL));
        break;
    case 2: /* real */
        if (!plain) {
            i_f = decomposition();
            break;
        }
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
    case 5: /* the bands that choose between the two as well as the two */
        if (((int)(periodic(zre * zre + zim * zim) * 10)) % 2)
            i_f = (cosl(periodic(zre * zim * pre * pim)) * 256 * 256);
        else
            i_f = (sinl(periodic(zre * zim * pre * pim)) * 256 * 256);
        break;
    case 6:
        i_f = ((zre * zre + zim * zim) * cosl(periodic(zre * zre)) * 256 *
               256);
        break;
    case 7:
        i_f = (sinl(periodic(zre * zre - zim * zim)) * 256 * 256);
        break;
    case 8:
        i_f = (atanl((zre * zim * pre * pim)) * 256 * 64);
        break;
    case 9: /* squares a fortieth of a unit on a side, at a frequency of one */
        if ((abs((int)(periodic(zre) * 40)) % 2) ^
            (abs((int)(periodic(zim) * 40)) % 2))
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
        i_f = ((nsin(periodic(zre)) * nsin(periodic(zim)) + 1) * 20000);
        break;
    case 19: /* sign(imag): two flat tones. A decomposition of the inside,
              * where binary decomposition does the same to the outside. */
        if (plain) {
            i_f = (zim > 0 ? 20000 : 45000);
            break;
        }
        /* the two tones as far apart as the speed says: at one the 25000 of
         * 20000 and 45000, at nought one tone */
        i_f = (decomp_binary_odd(
                   decomp_turn(qre, qim, cfractalc.indecomprotation),
                   cfractalc.indecompsectors, cfractalc.indecomprepeats,
                   cfractalc.indecompway) ^
               ring_odd())
                  ? 20000
                  : 20000 + cfractalc.indecompspeed * 25000;
        break;
    case 20: /* frac(mag): contours at every eighth of the modulus, at a
              * frequency of one */
        {
            number_t bands = periodic(nsqrt(zre * zre + zim * zim)) * 8;
            i_f = ((bands - nfloor(bands)) * 60000);
        }
        break;
    case INCOLORING_COS_ABS:
        /* cos(|z|): rings about where the orbit settled, eight to a unit of
         * the modulus at a frequency of one */
        i_f = ((ncos(periodic(nsqrt(zre * zre + zim * zim)) * 16 * N_PI) + 1) *
               20000);
        break;
    case INCOLORING_HEXAGONS:
        /* hexagons: three waves a third of a turn apart, six to a unit, which
         * add up to a lattice of hexagons */
        {
            static const number_t h = (number_t)0.86602540378443864676372317L;
            number_t s = ncos(periodic(zim) * 12 * N_PI) +
                         ncos(periodic(-h * zre - zim / 2) * 12 * N_PI) +
                         ncos(periodic(h * zre - zim / 2) * 12 * N_PI);
            i_f = ((s + (number_t)1.5) / (number_t)4.5 * 40000);
        }
        break;
    case INCOLORING_SPIRAL:
        /* spiral: five arms, wound once to every factor of e in the modulus.
         * The five stay five -- a number of turns that is not whole would
         * leave a seam where the angle wraps -- and the frequency and the
         * function act on the winding. */
        i_f = ((nsin(5 * natan2(zim, zre) +
                     periodic(nlog(zre * zre + zim * zim +
                                   (number_t)1e-30) / 2) *
                         2 * N_PI) +
                1) *
               20000);
        break;
    case INCOLORING_COS_ABS_ZC:
        /* cos(|z-c|): rings about the parameter, which tell the bulbs
         * apart as |z-c| does */
        i_f = ((ncos(periodic(nsqrt((zre - pre) * (zre - pre) +
                                    (zim - pim) * (zim - pim))) *
                     16 * N_PI) +
                1) *
               20000);
        break;
    case INCOLORING_TESSELLATION:
        /* tessellation: the tiles of one of randsctile's tilings, laid over z
         * forty to a unit as squares' squares are, through k*f as the
         * argument of every periodic function of the inside is. Each tile
         * has an amount of its own, the hash that names it, by which it is
         * told from its neighbours whatever their number and shape; a point
         * no tile is to be had for -- past the reach of a tiling that never
         * repeats -- has nought. See intessellcolor for what is coloured
         * and how. */
        {
            const int colour = cfractalc.intessellcolor;
            const int fill = cfractalc.intessellfill;
            /* where in its tile the point stands, measured only when
             * something reads it, the edge costing more than the name */
            double out = 0, rx = 0, ry = 0;
            const number_t turn = tessell_turn(
                cfractalc.intessellation, colour, fill, periodic(zre) * 40,
                periodic(zim) * 40, out, rx, ry);
            if (plain && colour == TESSELL_ANGLE && fill == TESSELL_TURNED &&
                cfractalc.intessellintensity == 1) {
                /* what the mode drew before it had a choice, to the bit */
                i_f = ((atan2l(zre, zim) / (M_PI + M_PI) + 0.75 + turn) *
                       20000);
                break;
            }
            i_f = tessell_colour(
                colour, fill, cfractalc.intessellintensity, turn, out, rx, ry,
                qre, qim,
                (zre * zre + zim * zim) * (number_t)(cfractalc.maxiter >> 1) *
                        SMUL +
                    SMUL,
                cfractalc.indecomprotation, cfractalc.indecompsectors,
                cfractalc.indecomprepeats, cfractalc.indecompway,
                cfractalc.indecompspeed, cfractalc.indecompradial, ring_odd);
        }
        break;
    case 21: /* log(mag): zmag with the contrast moved onto the small values,
              * which is where an attracting orbit spends its time */
        i_f = (nlog(zre * zre + zim * zim + (number_t)1 / 1000000) * 8000 +
               30000);
        break;
    };
    color_precalc(i_f, 1);
    iter = color_wrap(i_f, ((unsigned int)(cpalette.size - 1)) << 8);

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
/* the orbit multiplied by c: counted at the bailout, see SMOOTH_FAR */
#define SMOOTHNEAR
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
#define SMOOTHNEAR
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
#define SMOOTHNEAR
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
/* a squared modulus, carried past the bailout as the circle's is */
#define SMOOTHFAR (cfractalc.bailoutmode == BAILOUT_CIRCLE)
#define SMOOTHFARMEASURE (zpr * zpr + zpm * zpm)
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
/* escaping on the circle, carried past the bailout as the built-in ones are */
#define SMOOTHFAR (!newtok && cfractalc.bailoutmode == BAILOUT_CIRCLE)
#define SMOOTHFARMEASURE (zre * zre + zim * zim)
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
