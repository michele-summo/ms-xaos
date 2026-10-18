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

/* How many inside colouring modes there are, true colour included, and where
 * true colour sits: last, which three places used to spell as the number 10
 * and would have gone on spelling as 10 after modes were added in front of
 * it. */
#define INCOLORING 31
/* True colour used to be the last inside mode and was found by counting from
 * the end. The three fractional Brownian motion modes are written after it so
 * that every number already saved in a position still means what it meant, so
 * where it stands is now written down rather than counted. */
#define INCOLORING_TRUECOLOR 22
/* The fbm modes, which the flat list in the menu leaves out as it leaves out
 * true colour: they have a dialog of their own with the settings they take.
 * The three, in the order they are written in incolorname. */
#define INCOLORING_FBM 23
#define INCOLORING_FBM_ZMAG 24
#define INCOLORING_FBM_DECOMP 25
/* And the periodic ones written after them, for the same reason. */
#define INCOLORING_COS_ABS 26
#define INCOLORING_HEXAGONS 27
#define INCOLORING_SPIRAL 28
#define INCOLORING_COS_ABS_ZC 29
/* The tiles of one of randsctile's tilings, which took the place of hexagon
 * tiles before that was ever released: its number is free to take. */
#define INCOLORING_TESSELLATION 30
/* The ways the outside's fbm marks can be sized; see outfbmscaling. */
#define FBM_SCALE_FIXED 0
#define FBM_SCALE_ITER 1
#define FBM_SCALE_SMOOTH 2
#define OUTCOLORING (OutColormodeClass::ColOut_MAXMode + 1)
/* How many of them the flat list in the menu offers: the ones up to true
 * colour, which has a submenu of its own. The fbm modes are not among them --
 * they live in a submenu of their own with the settings they take. */
#define OUTCOLORING_LISTED 17
/* How many modes the per-formula symmetry tables below describe.
 *
 * They are written out mode by mode, one row apiece, and they were written
 * when there were ten inside modes and eleven outside ones with true colour
 * last of each. True colour has moved along as modes were added in front of
 * it; its row has not, so it is looked up where it has always been.
 *
 * A mode with no row of its own must not simply read the row that happens to
 * sit at its index: a row the compiler zero-filled does not mean "nothing
 * known" but "symmetric about both axes", and one belonging to another mode
 * means whatever that mode is symmetric about. Either would fold the picture.
 * combine_methods gives them no symmetry, which is always safe and forgoes an
 * optimisation rather than drawing the wrong thing. */
#define INCOLORING_DESCRIBED 10
#define OUTCOLORING_DESCRIBED 11

#define TCOLOR 15
/* One more than the speed functions colorfun[] names. New ones go on the end:
 * the number is what a saved position carries. */
#define COLORFUN 13

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
                    /* added after the others, and so before true colour,
                     * which must stay last: see INCOLORING_TRUECOLOR for the
                     * same rule on the inside modes */
                    ColOut_iter_plus_angle,
                    ColOut_iter_plus_log_mag,
                    ColOut_iter_plus_real_times_imag,
                    ColOut_max_real_imag,
                    ColOut_iter_banded,
                    ColOut_abs_real_minus_abs_imag,
                    ColOut_True_color,
                    /* Written after true colour, not before it. Every mode
                     * added before had pushed true colour along by one, and
                     * the mode is saved in a position as a number, so a
                     * picture saved with true colour came back as whatever
                     * had taken its place. These three leave every number
                     * where it is. */
                    ColOut_fbm_smooth,
                    ColOut_fbm_iter,
                    ColOut_fbm,
                    ColOut_MAXMode = ColOut_fbm
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
/* Appended rather than slotted in beside the other triangles: the number is
 * what a saved position carries, so inserting one would move every shape
 * after it and quietly change what old files draw. */
#define BAILOUT_TRIANGLE180 12
#define BAILOUTMODES 13
#define BAILOUT_MAXSIDES 8

/* The speed change after ITER of one side, as color_change works it: the
 * stretches of the count and how each is joined to the one before. In long
 * double whatever the build, as the speed functions are, so that the two
 * builds colour alike.
 *
 * kind      0 nothing to do; 1 one group in jump, an initial exponent of one
 *           and every speed one, worked the way 1.8.0 worked it -- iter,
 *           exponent -- so that a picture from then draws to the bit as it
 *           did; 2 the stretches
 * corners   how many ITERs cut the count -- nought, one or two
 * at        each corner's ITER
 * lo, hi    where its rounding begins and ends: the ITER less and plus the
 *           radius, or the ITER itself for a jump
 * from, pacelo, pacehi, bump
 *           for a rounding, the count where it begins, the paces at its two
 *           ends, and the height at ITER of the parabola added to the
 *           straight line of paces between them
 * C, K, J, e
 *           the stretches, before the first corner, between, past the last:
 *           C + K (v/J)^e */
/* What the two groups start as -- the initial exponent and the speeds at one,
 * the modes off. The first group's ITER and exponent are the ones 1.8.0 had;
 * the radius is how many iterations a smooth group rounds over on either side
 * of its ITER. */
#define COLORFUN_ITER 100
#define COLORFUN_EXP 2
#define COLORFUN_RADIUS 10
#define COLORFUN2_ITER 50
#define COLORFUN2_EXP 2
#define COLORFUN2_RADIUS 10

struct colorfun_plan {
    int kind;
    number_t iter, exponent;
    int corners;
    number_t at[2], lo[2], hi[2], from[2], pacelo[2], pacehi[2], bump[2];
    number_t C[3], K[3], J[3], e[3];
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
    /* And ln ln of the bailout, which the smooth count carried past it reads
     * at every pixel (see smooth_count_far), with the bailout it is of. */
    number_t bailoutloglog, bailoutloglogof;
    OutColormodeType coloringmode;
    int incoloringmode;
    int intcolor, outtcolor;
    //MSUMMO BEGIN HACK
    number_t incolorspeed, outcolorspeed;
    int incolorfun, outcolorfun;
    int incolorshift, outcolorshift;
    /* How often the inside modes that repeat do so: what the argument of
     * their periodic part is multiplied by -- the sines and cosines, the
     * squares, the contours of frac(mag). One, the default, draws them as
     * they always were, to the bit; the other modes never read it. */
    number_t incolorfrequency;
    /* And the speed function, by its number in colorfun[], that the argument
     * of every periodic function they use -- the sines and cosines, frac, the
     * parity of the squares and of the bands of cos(mag) -- is bent by before
     * the frequency multiplies it: k*f(a). Nought, ident, leaves it as it
     * was; the modes that do not repeat never read it. */
    int incolorperiodicfun;
    /* The decompositions' settings, a side at a time: binary decomposition
     * and color decomposition outside, sign(imag), Decomposition-like and
     * fbm + decomposition inside.
     *
     * sectors   how many each repeat of the turn is cut into. Nought, where
     *           it starts, is the mode's own: the two half-planes of the
     *           binary ones, the unbroken turn of the others. The binary ones colour alternate
     *           sectors, so an odd number is taken up to the even one above
     *           it: two sectors of a colour would meet.
     * rotation  the degrees the sectors are turned by, anticlockwise.
     * of        inside only, whose angle: z, z-c or z*c.
     * radial    bands across the sectors, alternate ones told apart. Outside
     *           a band is that many iterations; inside that many rings to a
     *           unit of the modulus of what the angle is of, which goes
     *           through k*f as the argument of every periodic function of
     *           the inside does. Nought, none.
     *
     * And how the turn is walked:
     *
     * repeats   how many times the turn repeats, the sectors in each and the
     *           palette across each; any number, a part of one leaving a seam
     *           where the angle wraps.
     * way       0 right, as it always went; 1 left, the other way round; 2 both,
     *           forwards and back again in each repeat, a mirror with no seam
     *           where repeats is whole.
     * speed     for the colour ones, how much of the palette a repeat spans:
     *           more colours in each, not more of them. For the binary ones,
     *           how far apart the colours of alternate sectors are: one as
     *           they always were, nought none, two twice.
     *
     * The cut at nought and the walk at one, right, one, the modes draw as they
     * always did, to the bit. */
#define DECOMP_MAXSECTORS 3600
    int indecompsectors, outdecompsectors;
    number_t indecomprotation, outdecomprotation;
    int indecompof;
    int indecompradial, outdecompradial;
    number_t indecomprepeats, outdecomprepeats;
    int indecompway, outdecompway;
    number_t indecompspeed, outdecompspeed;
    /* The speed change after ITER, a side at a time as the speed is: the
     * count handed to the speed function is cut at up to two ITERs and each
     * stretch driven by an exponent of its own, so that the colours speed up
     * or slow down from an ITER and not before.
     *
     * initial    the exponent from the start to the first ITER of a group
     *            that is on, or for every count with both off: the plain
     *            power v^initial, the same whether the groups are on or off.
     *            One, the default, leaves it as it was
     * initialspeed
     *            what the first stretch's pace is multiplied by:
     *            initialspeed v^initial. One, the default, leaves it as it was
     *
     * and two groups, fun and fun2, of
     *
     * iter       where the group's stretch begins, in the units the function
     *            is handed -- iterations, for the modes that count them.
     *            Nought or below puts it at the start: a plain power.
     * exponent   the stretch's exponent, from iter to the next ITER or on.
     *            The stretch starts from the count the one before it came to
     *            and adds iter ((v/iter)^exponent - 1), whose pace is the
     *            exponent at iter: above one the colours speed up from
     *            there, below one they slow down, at one they go on at a
     *            pace of one
     * speed      what the stretch's pace is multiplied by, so that it adds
     *            speed iter ((v/iter)^exponent - 1) -- at an exponent of one a
     *            straight line of slope speed. One leaves it as it was
     * mode       0 off, which is the default and leaves everything as it
     *            was; 1 jump, the pace turning at iter at once; 2 smooth, the
     *            pace going over from the one at iter - radius to the one at
     *            iter + radius, and the count the jump's again past it. The
     *            names are colorfun_modes[]; 1.8.0's four are read into these
     * radius     for smooth, how many iterations it rounds over on either
     *            side of iter
     *
     * A group that is off is left out, its ITER and exponent with it. The
     * stretches go in the order of their ITERs, whichever group holds the
     * lower; at the same ITER the second group replaces the first. Each
     * exponent is its stretch's own, never multiplied by the one before. See
     * color_change in formulas.cpp. */
    number_t incolorfuninitial, outcolorfuninitial;
    number_t incolorfuninitialspeed, outcolorfuninitialspeed;
    number_t incolorfuniter, outcolorfuniter;
    number_t incolorfunexp, outcolorfunexp;
    number_t incolorfunspeed, outcolorfunspeed;
    int incolorfunmode, outcolorfunmode;
    number_t incolorfunradius, outcolorfunradius;
    number_t incolorfun2iter, outcolorfun2iter;
    number_t incolorfun2exp, outcolorfun2exp;
    number_t incolorfun2speed, outcolorfun2speed;
    int incolorfun2mode, outcolorfun2mode;
    number_t incolorfun2radius, outcolorfun2radius;
    /* What the two groups make of the count, worked out once a picture by
     * color_change_prepare and not saved. */
    struct colorfun_plan inplan, outplan;
    /* The fractional Brownian motion the fbm colouring modes are drawn from,
     * kept apart for the two sides as the speed and the shift are.
     *
     * intensity  how many bands of colour the noise moves the value by
     * frequency  cells of the lattice to a unit of the plane: the scale of
     *            the marks, and what has to be raised as one zooms in
     * octaves    how many are summed, each at twice the frequency
     * roughness  what each octave keeps of the one before it; a half is the
     *            plain fractional Brownian motion, more is grittier
     * seed       the same picture every time, and a different one on demand
     */
    number_t infbmintensity, outfbmintensity;
    number_t infbmfrequency, outfbmfrequency;
    number_t infbmroughness, outfbmroughness;
    int infbmoctaves, outfbmoctaves;
    int infbmseed, outfbmseed;
    /* How the marks of the outside's motion are sized: by the frequency
     * alone (FBM_SCALE_FIXED, where it starts and what every position saved
     * before there was a choice draws), or by the frequency raised with the
     * count, whole (FBM_SCALE_ITER) or smooth (FBM_SCALE_SMOOTH), so that
     * the marks grow finer where the orbit took longer to leave. factor is
     * how many times the frequency doubles over ten iterations; at nought
     * the count changes nothing. The inside has no count to go by -- every
     * point of it ran to maxiter -- and so has neither. */
    int outfbmscaling;
    number_t outfbmfactor;
    /* Which of randsctile's tilings, 1 to 45, the inside's tessellation mode
     * lays over z; hexagons where it starts. */
    int intessellation;
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

/* Whether a side's decomposition settings are where they start, and the
 * decompositions draw as they always did: the old code, to the bit, and the
 * symmetry their rows of the formula tables describe. */
static inline int decomp_plain(const fractal_context *f, int inset)
{
    if (inset)
        return !f->indecompsectors && f->indecomprotation == 0 &&
               !f->indecompof && !f->indecompradial &&
               f->indecomprepeats == 1 && !f->indecompway &&
               f->indecompspeed == 1;
    return !f->outdecompsectors && f->outdecomprotation == 0 &&
           !f->outdecompradial && f->outdecomprepeats == 1 &&
           !f->outdecompway && f->outdecompspeed == 1;
}
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
/* Where the pixel being coloured stands on the plane.
 *
 * The colouring is handed the point the orbit left by and not the point it
 * started from, and the fbm modes want the second: the escape point jumps
 * from pixel to pixel, so noise read at it comes out as grain rather than as
 * marks. The parser keeps the same thing in sffe_position, but only when the
 * parser is built in, and a colouring mode may not depend on that. */
extern thread_local number_t color_px, color_py;

void set_fractalc(fractal_context *, struct image *img);

/* Works out inplan and outplan from the two groups of each side. set_fractalc
 * calls it once a picture; anything else that colours with a context it has
 * filled in itself has to as well. */
void color_change_prepare(fractal_context *c);

/* A group in the four modes positions named before the radius -- off, jump,
 * continuous, smooth -- turned into the three there are: mode, and exponent
 * and radius where they have to change. See formulas.cpp. */
void colorfun_from_legacy(number_t iter, number_t *exponent, int *mode,
                          number_t *radius);

/* The cell of the palette one point of the plane is drawn from, or nought for
 * the inside colour, which is no cell. Returns 0 when it cannot be had. See
 * fractal.cpp. */
int fractal_palette_cell(fractal_context *context, struct image *img,
                         number_t x, number_t y, int *cell);
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
