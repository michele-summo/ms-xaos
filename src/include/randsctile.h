#ifndef RANDSCTILE_H
#define RANDSCTILE_H

#include <cstdint>

#include "config.h"

/* How many tilings randsctile knows, numbered from one. */
#define RANDSCTILE_TILINGS 45

/* The tile of tiling kind -- randsctile's first argument, 1 to 45 -- that
 * holds the point (x, y), in the units randsctile lays its tiles of unit
 * area on: a hash that names it, the same for every point of the tile and
 * for no other tile near it. Says whether there was one: not for a kind
 * outside 1 to 45, nor past the reach of the tiling, where randsctile goes
 * flat. For the inside's tessellation mode, which colours by the tiles
 * without the parser's arguments around them. */
int randsctile_name(int kind, number_t x, number_t y, uint64_t *name);

/* The same tile and, besides its name, where the point stands in it: out,
 * nought at the tile's deepest and one at its edge, and (rx, ry), the point
 * seen from the tile's middle -- what randsctile's skew and rosette read. It
 * costs more than the name alone, the edge being measured, so it is for when
 * they are wanted. */
int randsctile_place(int kind, number_t x, number_t y, uint64_t *name,
                     double *out, double *rx, double *ry);

#endif
