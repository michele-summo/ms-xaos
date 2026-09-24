/* What the thumbnails in the Tilings tab of the formula reference show, shared
 * by the program that draws them (randsctile-thumbnails.cpp) and by the test
 * that checks they are still what randsctile draws (formula_help_test.cpp).
 *
 * The pictures are kept in the tree rather than drawn each time the tab is
 * opened, and a picture kept is a copy that can go stale: a tiling redrawn or
 * a hash changed, and the gallery goes on showing the old one. So each is
 * stored with a fingerprint of what randsctile drew over its square, and the
 * test draws the same points again and compares.
 */
#ifndef RANDSCTILE_THUMBNAILS_H
#define RANDSCTILE_THUMBNAILS_H

#include <cstdint>
#include <cstdio>

/* The square each thumbnail shows: eight units a side centred on the origin,
 * which with one tile to the unit of area is some sixty tiles -- enough to see
 * how the shapes repeat, and the largest stars whole. y goes up, as the plane
 * does in XaoS. */
#define RANDSCTILE_THUMB_SPAN 8
/* Pixels a side as stored: twice the 128 the tab shows, so that a screen
 * scaled up to 200% is still drawn from something as fine as itself. */
#define RANDSCTILE_THUMB_SIDE 256
/* Samples a pixel, each way, which is what smooths the outlines. */
#define RANDSCTILE_THUMB_SAMPLES 4
/* The seed only colours the tiles, and any would do; this one is the tests'. */
#define RANDSCTILE_THUMB_CALL "randsctile(%d,7)"

/* How many samples a side the whole square is drawn from. */
#define RANDSCTILE_THUMB_FINE (RANDSCTILE_THUMB_SIDE * RANDSCTILE_THUMB_SAMPLES)

/* The point sample (row, column) of the fine grid stands for. */
static inline void randsctile_thumb_point(int row, int column, double *x,
                                          double *y)
{
    *x = ((column + 0.5) / RANDSCTILE_THUMB_FINE - 0.5) * RANDSCTILE_THUMB_SPAN;
    *y = (0.5 - (row + 0.5) / RANDSCTILE_THUMB_FINE) * RANDSCTILE_THUMB_SPAN;
}

/* One sample in sixteen each way: 64 by 64 points an eighth of a unit apart,
 * a tenth of a second for all forty-five where the pictures take half a
 * minute, and fine enough that a tiling redrawn or a hash changed cannot slip
 * between them all. */
#define RANDSCTILE_THUMB_PRINT_STEP 16

/* The fingerprint of a thumbnail: FNV-1a over the values randsctile gives at
 * the sparse grid above, each cut to 24 bits. The values are 53-bit fractions
 * and the same at every precision; the cut is there so that the rounding of
 * anything they are later summed with can never tell the two builds apart.
 * value_at is the call through the parser, passed in so that this header
 * needs none of the parser's own. */
template <typename Eval>
static inline uint64_t randsctile_thumb_fingerprint(Eval value_at)
{
    uint64_t h = 1469598103934665603ULL;
    for (int r = RANDSCTILE_THUMB_PRINT_STEP / 2; r < RANDSCTILE_THUMB_FINE;
         r += RANDSCTILE_THUMB_PRINT_STEP)
        for (int c = RANDSCTILE_THUMB_PRINT_STEP / 2; c < RANDSCTILE_THUMB_FINE;
             c += RANDSCTILE_THUMB_PRINT_STEP) {
            double x, y;
            randsctile_thumb_point(r, c, &x, &y);
            double v = value_at(x, y);
            uint32_t q = v <= 0 ? 0 : v >= 1 ? 0xFFFFFF : (uint32_t)(v * 16777216.0);
            for (int b = 0; b < 3; b++) {
                h ^= (q >> (8 * b)) & 0xFF;
                h *= 1099511628211ULL;
            }
        }
    return h;
}

#endif /* RANDSCTILE_THUMBNAILS_H */
