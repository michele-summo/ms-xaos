#ifndef FBM_NOISE_H
#define FBM_NOISE_H

/* The fractional Brownian motion, one for the whole program: the colouring
 * modes of that name (engine/formulas.cpp) and fbm() in the parser
 * (sffe/sffe_cmplx_gsl.cpp) both draw it from here.
 *
 * Gradient noise: every corner of a lattice is given a direction, picked by a
 * hash of the corner, and a point takes from each corner the slope that
 * direction gives at its distance, the four blended by a quintic that is flat
 * at both ends. Octaves of it are summed, each at twice the frequency and
 * keeping a share of the height of the one before -- the roughness -- which is
 * what makes it read as wear rather than as a pattern.
 *
 * It was value noise, a number at each corner blended by a smoothstep, and
 * that draws a grid. The smoothstep has no slope at either end, so the field
 * goes flat along every line of the lattice; and every octave's lattice held
 * the lines of the first, so all of them went flat there together. Measured
 * over FBM_ERROR.xpf, the field was a third as steep on the lines as between
 * them: squares of eight pixels, a picture that looked worked out at a lower
 * resolution than it was shown at. Gradient noise has the same slope on the
 * lines as between them -- measured, 0.99 of it at the default settings and
 * 1.08 at those of FBM_ERROR.xpf, where across a line it had been nought --
 * and each octave is shifted against the one before, so no two share a line
 * or a corner either.
 *
 * Five things are kept from before, so that the numbers saved with a picture
 * go on meaning what they meant:
 *
 * - the value runs from nought to one and never outside, which is what lets
 *   the callers add it to something without taking that under nought;
 * - the contrast: one octave of gradient noise swings only seven tenths as far
 *   from its middle as one of value noise, and fbm_noise_contrast takes it
 *   back onto the distribution value noise had, octave by octave;
 * - the marks are the size they were: gradient noise at the same lattice
 *   draws them smaller, and the lattice runs at FBM_NOISE_SCALE of the
 *   frequency asked for so that they come out as large as the value noise
 *   drew them -- measured by where the field stops resembling itself, 0.56
 *   of a cell at the default settings, before and after;
 * - the seed alone decides the marks, with no clock and no pass in it;
 * - and the cost, near enough: the colouring modes take 159 ns a pixel at four
 *   octaves and 452 at eighteen, against 160 and 479; fbm() 182 and 498
 *   against 161 and 438; randsc 167 against 147 -- the best of eight runs
 *   each, old and new taken by turns.
 */

#include <cstdint>

#include "config.h"

/* Past this there is no cell to stand in: a quarter of what an int64 holds,
 * as randsc keeps, which leaves room for the doubling between octaves. Every
 * octave after the first that reaches it is finer still, so what has been
 * summed by then is all there is to have. */
#define FBM_NOISE_LIMIT ((number_t)2.0e18)
#define FBM_NOISE_INDEX_LIMIT ((int64_t)2000000000000000000LL)

/* The lattice, in cells to a unit of the frequency asked for. See above. */
static const number_t FBM_NOISE_SCALE = (number_t)3 / 5;

/* The step from one octave to the next: twice the point, shifted by these.
 * They are 2 - phi and phi - 1, and being irrational no number of doublings
 * brings a line or a corner of one octave onto one of another.
 *
 * The octaves were turned against one another as well, by the angle of a
 * 3-4-5 triangle, as is often done to keep the lattice from showing. It did
 * nothing here that could be measured: over 400000 points the slope differed
 * by a tenth of a per cent from one direction to another, turned or not --
 * sixteen directions are round enough already. And it cost twelve
 * nanoseconds an octave, since a turn has to divide the cell by five where a
 * doubling only adds it to itself. */
static const double FBM_NOISE_SHIFT_X = 0.381966011250105151795;
static const double FBM_NOISE_SHIFT_Y = 0.618033988749894848205;

/* Sixteen directions, half a step off the axes.
 *
 * Everything inside a cell is worked in double, whatever the build: the
 * place in the cell is a fraction, which double holds to far below anything a
 * colour can show, and long double arithmetic is the x87's, which made the
 * motion a quarter to two fifths slower than the value noise it replaced
 * (207 ns against 164 at four octaves, 688 against 483 at eighteen). What the
 * zoom needs the precision for is which cell, and that is an integer. */
static const double FBM_NOISE_DIR[16][2] = {
    {+0.980785280403230430579, +0.195090322016128275839},
    {+0.831469612302545235671, +0.555570233019602177649},
    {+0.555570233019602177649, +0.831469612302545235671},
    {+0.195090322016128275839, +0.980785280403230430579},
    {-0.195090322016128275839, +0.980785280403230430579},
    {-0.555570233019602177649, +0.831469612302545235671},
    {-0.831469612302545235671, +0.555570233019602177649},
    {-0.980785280403230430579, +0.195090322016128275839},
    {-0.980785280403230430579, -0.195090322016128275839},
    {-0.831469612302545235671, -0.555570233019602177649},
    {-0.555570233019602177649, -0.831469612302545235671},
    {-0.195090322016128275839, -0.980785280403230430579},
    {+0.195090322016128275839, -0.980785280403230430579},
    {+0.555570233019602177649, -0.831469612302545235671},
    {+0.831469612302545235671, -0.555570233019602177649},
    {+0.980785280403230430579, -0.195090322016128275839}};

/* Two-dimensional gradient noise with directions of unit length reaches the
 * square root of a half at most, so this takes it onto nought to one. */
static const double FBM_NOISE_GAIN = 0.70710678118654752440;

/* The contrast curve, over the value, nought to one in sixty-four steps: in
 * each, the cubic in the place x within the step, lowest power first.
 *
 * It is the map that takes the values of one octave of gradient noise onto
 * the distribution one octave of the old value noise had -- measured over
 * eight million points by tools/fbm-contrast-table.py, which writes it -- as
 * a monotone cubic Hermite table, so that it never turns back on itself and
 * never leaves nought to one. A polynomial fitted to the same map did both,
 * by a thousandth, near the ends. Being monotone and smooth it moves no
 * outline and draws no crease: the field is the same shapes, only further
 * apart in value.
 *
 * Written out over the value rather than over the distance from the middle,
 * where it was fitted, so that looking it up costs a multiplication, a
 * truncation and three steps of Horner. Over the distance it wanted the
 * distance, the side and the Hermite basis as well, and at seventy cycles
 * one after another it came to seventeen nanoseconds an octave.
 *
 * The table, and what it measured, are in fbm_contrast_table.h. */
#include "fbm_contrast_table.h"

static inline double fbm_noise_contrast(double v)
{
    double s = v * 64;
    int k = (int)s; /* v is not below nought, so this is its floor */
    k = k < 63 ? k : 63;
    double x = s - k;
    const double *c = FBM_NOISE_CURVE[k];
    return c[0] + x * (c[1] + x * (c[2] + x * c[3]));
}

static inline uint64_t fbm_noise_mix(uint64_t h)
{
    h ^= h >> 33;
    h *= 0xFF51AFD7ED558CCDULL;
    h ^= h >> 33;
    h *= 0xC4CEB9FE1A85EC53ULL;
    h ^= h >> 33;
    return h;
}

/* Each octave has a key of its own, so that none agrees with another. */
static inline uint64_t fbm_noise_key(int octave, uint64_t seed)
{
    return fbm_noise_mix((uint64_t)octave * 0x9E3779B97F4A7C15ULL ^ seed);
}

/* floor(x) as an index, with the fraction left over, or 0 past the limit. */
static inline int fbm_noise_cell(number_t x, int64_t *cell, number_t *frac)
{
    if (!(x > -FBM_NOISE_LIMIT && x < FBM_NOISE_LIMIT))
        return 0;
    int64_t c = (int64_t)x;
    number_t t = (number_t)c;
    if (t > x) {
        c--;
        t -= 1;
    }
    *cell = c;
    *frac = x - t;
    return 1;
}

/* The slope a corner's direction gives at (dx, dy) from it. */
static inline double fbm_noise_slope(int64_t cx, int64_t cy, uint64_t key,
                                     double dx, double dy)
{
    uint64_t h = fbm_noise_mix((uint64_t)cx * 0x9E3779B97F4A7C15ULL ^
                               (uint64_t)cy * 0xC2B2AE3D27D4EB4FULL ^ key);
    const double *d = FBM_NOISE_DIR[h >> 60];
    return d[0] * dx + d[1] * dy;
}

/* One octave, nought to one, in the cell (cx, cy) at (u, v) inside it, as
 * gradient noise gives it: before the contrast curve. */
static inline double fbm_noise_raw(int64_t cx, int64_t cy, double u, double v,
                                   uint64_t key)
{
    double fu = u * u * u * (u * (u * 6 - 15) + 10);
    double fv = v * v * v * (v * (v * 6 - 15) + 10);
    double a = fbm_noise_slope(cx, cy, key, u, v);
    double b = fbm_noise_slope(cx + 1, cy, key, u - 1, v);
    double c = fbm_noise_slope(cx, cy + 1, key, u, v - 1);
    double d = fbm_noise_slope(cx + 1, cy + 1, key, u - 1, v - 1);
    double lo = a + (b - a) * fu;
    double hi = c + (d - c) * fu;
    return 0.5 + (lo + (hi - lo) * fv) * FBM_NOISE_GAIN;
}

/* One octave with the contrast value noise had, which is what randsc is. */
static inline double fbm_noise_octave(int64_t cx, int64_t cy, double u,
                                      double v, uint64_t key)
{
    return fbm_noise_contrast(fbm_noise_raw(cx, cy, u, v, key));
}

/* A cell and a place in it, with the place brought back into [0, 1).
 * Without a branch: the place is as likely a little under nought as not, and
 * a branch on it would be mispredicted every other time. */
static inline void fbm_noise_settle(int64_t *cell, double *frac)
{
    int64_t whole = (int64_t)*frac;
    whole -= (double)whole > *frac; /* truncation went up: one below */
    *cell += whole;
    *frac -= (double)whole;
}

/* The next octave's point: twice this one, shifted. Worked on the cell and
 * the place apart, so the cell is integer arithmetic and exact however far
 * out, and only the place is a double. Doubling the whole coordinate in
 * number_t, as value noise did, is exact too, but it is the x87's arithmetic
 * and a shift after it would round the place away far enough out -- long
 * double holds a coordinate of 10^15 to a ten-thousandth of a cell.
 *
 * Nothing overflows below FBM_NOISE_INDEX_LIMIT: twice it is inside what an
 * int64 holds. */
static inline void fbm_noise_step(int64_t *cx, int64_t *cy, double *u,
                                  double *v)
{
    *cx *= 2;
    *cy *= 2;
    *u = 2 * *u + FBM_NOISE_SHIFT_X;
    *v = 2 * *v + FBM_NOISE_SHIFT_Y;
    fbm_noise_settle(cx, u);
    fbm_noise_settle(cy, v);
}

/* The motion at (x, y), with x and y already multiplied by the frequency:
 * nought to one. octaves is held between one and twenty-four, and a
 * roughness that is not above nought is taken as a half. */
static inline number_t fbm_noise(number_t x, number_t y, uint64_t seed,
                                 int octaves, number_t rough)
{
    if (octaves < 1)
        octaves = 1;
    if (octaves > 24)
        octaves = 24;
    if (!(rough > 0))
        rough = (number_t)1 / 2;

    int64_t cx, cy;
    number_t pu, pv;
    if (!fbm_noise_cell(x * FBM_NOISE_SCALE, &cx, &pu) ||
        !fbm_noise_cell(y * FBM_NOISE_SCALE, &cy, &pv))
        return 0;
    double u = (double)pu, v = (double)pv;

    /* The octaves first and the curve after, in a loop of its own. The
     * curve is a chain of some forty cycles, each step waiting on the last;
     * put inside the loop that walks the octaves it waited on the octave
     * too, and the processor could not start one octave's curve until it
     * was nearly done with the one before. Apart, the curves of all the
     * octaves are independent and run side by side. Measured over fbm() at
     * eighteen octaves: 614 ns the one way, 498 the other, and 465 with no
     * curve at all. */
    double raw[24];
    int n = 0;
    for (;;) {
        raw[n] = fbm_noise_raw(cx, cy, u, v, fbm_noise_key(n, seed));
        if (++n == octaves)
            break;
        fbm_noise_step(&cx, &cy, &u, &v);
        if (!(cx > -FBM_NOISE_INDEX_LIMIT && cx < FBM_NOISE_INDEX_LIMIT &&
              cy > -FBM_NOISE_INDEX_LIMIT && cy < FBM_NOISE_INDEX_LIMIT))
            break;
    }
    double sum = 0, amp = 1, norm = 0, share = (double)rough;
    for (int i = 0; i < n; i++) {
        sum += amp * fbm_noise_contrast(raw[i]);
        norm += amp;
        amp *= share;
    }
    return (number_t)(sum / norm);
}

#endif /* FBM_NOISE_H */
