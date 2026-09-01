# MS XaoS

A fork of [XaoS](https://github.com/xaos-project/XaoS) 4.3.3. Everything the
original does, it still does; this describes what has been added or changed,
and why. Version 1.1.

## Two binaries

`XaoS.exe` computes at 64 bits of mantissa, `XaoS-quad.exe` at 113. Both are
built from the same sources — the precision is a compile-time type, not a
setting — and the quad one calls itself **MS XaoS Quad** so the two are
distinguishable once running.

Quad zooms about fifteen orders of magnitude deeper. Rendering the same
position at 1e-18: the ordinary binary keeps 21% of its detail, the quad one
92%. It is not the default because `__float128` is emulated in software:
roughly eight times slower on the built-in formulas and fifteen on
user-defined ones.

A position saved by the quad binary carries `(precision 113)`. Opening it with
a narrower build says so and draws it anyway; opening it with a version that
does not know the command is refused, which is the intent — a picture computed
at 113 bits cannot be reproduced at 64, and failing to open is more honest than
drawing something else.

## Saved positions

The view is written at the full precision of the build — 21 significant digits
ordinarily, 36 for quad. It used to be written with a digit count derived from
the zoom depth, which under-reported and capped at 17, so a deep position
reopened a fraction of a pixel from where it was saved. At a high iteration
count that is a visibly different picture.

Everything a position can hold now survives a save and load unchanged: all 169
shipped examples and 141 generated positions covering every formula, plane,
colouring mode and filter round-trip pixel for pixel.

## User formulas

**Reference.** Help → User formula reference lists every function, variable and
notation the formula language accepts, in a window that can be left open beside
the formula being written. It is grouped by what the functions are for, with the
special functions — `erf`, `gamma`, `lambertw` — in a section of their own,
since none of them is elementary and one does not reach for them by accident.
The list is checked against the parser's own table by a test, so it cannot
drift.

**`ifiterf(a; b)`** — evaluates `a` on every pass but the final one and `b`
on that. The final pass is the last the iteration limit allows: a formula has
no way of knowing which pass will be the one that escapes, that depending on
the value it has not produced yet. Only the chosen one is evaluated.

**`ifiterr(a; b; n)`** — evaluates `a` while the pass number is below `n`
and `b` from `n` onwards. Only the chosen one is evaluated, as with the other
two, though the threshold had to be taught to the parser first: an argument
may now be marked as read by the selector rather than chosen by it, and is
then evaluated before the choice is made.

**Removed.** `powi`, `powdc` and `logcn` were second names for `pow` and
`logn`; `rad`, `deg` and `sign` were listed with no implementation behind them,
and a second `trunc` was shadowed by the working one. A formula using a removed
alias should use the name that remains — the function is identical.

**`erf(z)`** — the error function over the complex plane. Accurate to about
three ulp inside a bailout of two, where a fractal actually iterates.

**`randsc(seed; size; degradation; kaleidoscope; mode)`** — coherent
noise over the point, giving blobs rather than per-pixel snow. `size`
(default `1+i`) is the average width of a blob along the real axis and its
height along the imaginary one.
`degradation` (default `1+i`) shrinks them as the iteration proceeds: the
size is multiplied by it at every pass, component by component, so `0.5+0.2i`
over `1+i` gives `1+i` on the first pass, then `0.5+0.2i`, then `0.25+0.04i`.
A zero in either component of either argument returns zero rather than dividing
by zero. The last two are the kaleidoscope, below. Only the seed is required.

**`randscq(...)`** — the same field without the interpolation: a mosaic of flat
square cells instead of blobs. Same arguments, same meaning.

**`randscp(...)`** — the same field again with the curves taken out but not the
irregularity: one seed is scattered inside each cell and every position takes
the value of the nearest seed, which draws a Voronoi diagram — flat convex
polygons with straight edges, no two the same shape.

**`randsch(...)`** — the same again on hexagons: a honeycomb of flat cells. Of
the regular polygons that tile the plane this is the one without a grain —
every cell has the same six neighbours at the same six angles — so it reads as
a material rather than as a grid.

**`randsct(...)`** — equilateral triangles, alternating in orientation. A
triangular mosaic does have a grain, and that is what one asks for by choosing
it.

| | |
| --- | --- |
| `randsc` | soft, curved |
| `randscq` | hard, regular, square |
| `randscp` | hard, irregular |
| `randsch` | hard, regular, hexagonal |
| `randsct` | hard, regular, triangular |

All five lay one cell over each unit square of the size in force, so `size`
means the same thing throughout: changing one letter changes the shape of the
cells and not the scale of the picture. A hexagon of circumradius one covers
2.6 unit squares and a triangle of side one covers 0.43, so those two grids are
scaled to match the rest rather than left as they come.

Each pass gets a field of its own. The iteration is hashed along with the
position, so a formula calling one of these once per iteration gets a new
value every time — including at a degradation of one, which leaves the size
alone. Degradation sets how big the blobs are on a given pass; it is not
what makes the pass different.

There is a floor under the size. A cell is found by dividing the position by
the size and the answer has to land in an integer, so the size cannot usefully
go below about the position over what an integer holds — a billionth of a
billionth of it. A degradation of a half reaches that in some sixty passes.
Past it there is no cell structure left to resolve and the field is one flat
value over the whole plane; it still changes on every pass and still differs
between the five functions, so a formula subtracting one from another does not
settle on zero and iterate to the limit for nothing. For a fade that stays a
picture the whole way, use a degradation near one.

All of them hash the position of the point and the iteration, never `z`, and
never any global state. So the same picture comes back on every redraw, at
any thread count, and in Mandelbrot or Julia mode alike. `randsc` is
continuous, so the two precisions agree to about 1e-19; the mosaics are step
functions and they can disagree on a hairline along the cell edges, which is
inherent in asking for hard edges.

A formula calling any of them turns boundary tracing off. That optimisation
fills a region it found one colour around without computing it — true of a
fractal, false of a noise field.

The old `rand` is unchanged and still depends on how many times it has been
called, so the same position does not redraw the same way. Prefer `randsc`.

### Kaleidoscopes

The two arguments after the degradation fold the plane before the field is
sampled from it. The first is how many wedges to fold it into, `1` — the
default — leaving it alone; the second is which mirror does the folding:

| mode | |
| --- | --- |
| `0` (and anything else) | the far half of each wedge mirrors the near half, so every wedge is symmetric about its own bisector |
| `1` | the same the other way about, the near half mirroring the far one |

    randsc(13;{0.6,0.6};{1,1};6;0)     six wedges, each a mirror of itself
    randsc(13;{0.6,0.6};{1,1};3;1)     three wedges, folded the other way

Both folds are continuous where the wedges meet, so the noise stays
coherent and the two precisions go on agreeing; a fold that met itself unevenly
would show as a seam.

This is the only part of the family that costs trigonometry, and only a level
of two or more reaches it. A call that says nothing about it pays one
comparison.

### Getting a picture out of the noise

`randsc` returns a value in `[0, 1)`. Used alone the iteration never leaves a
bailout of 4, so nothing escapes and the image is flat. Multiply it, or let it
perturb an iteration that does escape:

    z+randsc(13;{0.35,0.35};{1,1})*0.25       blobs
    randsc(13;{0.15,0.15};{0.5,0.5})*2.5      brownian motion
    randscq(13;{0.15,0.15};{0.5,0.5})*2.5     scattered squares, shrinking
    z+randscp(13;{0.3,0.3};{1,1})*0.25        irregular polygons
    z+randsch(13;{0.3,0.3};{1,1})*0.25        a honeycomb
    z+randsct(13;{0.3,0.3};{1,1})*0.25        triangles
    z^2+c+randsc(13;{0.25,0.25};{1,1})*1.2    the set itself deformed

A degradation of 0.5 halves the blobs every pass, so after twenty iterations
they are a millionth of their size and below a pixel. For a slow fade over a
long run, use something near 1 — `0.97^60` is still 0.16.

An integer seed is exact in both builds; a fractional one is quantised to its
leading bits, which agree except for about one seed in fifty million.

## Watching the orbit

`trap` and `stripe` keep one number about a whole orbit rather than about the
point, which is the kind of quantity the engine's colouring modes cannot reach:
a calculation loop there is compiled to hand back a colour, and only the last
`z` and the one before it survive to be coloured by. Done in the formula it
costs the engine nothing.

Both hand back the value they were given until the last pass the iteration
limit allows, and what they gathered on that one. So a whole formula is

    trap(z^2+c; 3)

and nothing else: the fractal iterates as it would, and on the last pass the
value becomes the trap, which the inside colouring modes then draw. A point
that escapes never reaches that pass and keeps its ordinary outside colour, so
what these draw is the inside.

**`trap(a; shape; centre; size)`** — the nearest the orbit ever came to a
shape. `shape` defaults to 0, `centre` to the origin, `size` to 1.

| shape | |
| --- | --- |
| `0` | the centre itself, a point |
| `1` | a horizontal line through it |
| `2` | a vertical one |
| `3` | both, a cross |
| `4` | a ring of radius `size` |
| `5` | a square of half-side `size` |
| `6` | a diamond of half-diagonal `size` |

**`stripe(a; density)`** — the average of `(sin(density * arg a) + 1) / 2`
along the orbit. `density` defaults to 4 and is how many stripes go round a
turn; a whole number, or they do not meet where the turn closes. An average
over the orbit changes smoothly with the point even where the iteration count
jumps, which is what draws the fibres the method is known for.

    trap(z^2+c;0)                    how near the orbit passed the origin
    trap(z^2+c;3)                    a cross, which draws rays
    trap(z^2+c;4;{0,0};{0.5,0})      rings inside the cardioid
    stripe(z^2+c;6)                  six stripes to a turn

The running quantity lives on the call site, so two traps in one formula keep
their own and a thread cannot disturb another. Both turn boundary tracing off,
for the same reason the noise does: two neighbours that take the same number of
passes can still have seen quite different orbits, and a region filled without
being computed would be wrong.

These are not colouring modes, and cannot be: the modes live inside each
formula's compiled loop, which returns a colour directly, and there are
thirty-one of those with a positional symmetry table apiece. What is here works
for a user formula only.

## Bailout shape

Calculation → Bailout → Bailout mode. The iteration stops when `z` leaves a
region, which has always been a circle. The shape is visible in the bands
outside the set, so it is as much a drawing tool as a numerical one.

| shape | escapes when |
| --- | --- |
| Circle (classic) | \|re\|² + \|im\|² passes the bailout |
| Square | either component alone does |
| Diamond | \|re\| + \|im\| does |
| Real axis | the real component alone |
| Imaginary axis | the imaginary component alone |
| Both axes | both components at once |
| Triangle 0° / 90° / 180° / −90° | the point leaves a triangle of that orientation |
| Hexagon 0° / 90° | the point leaves a hexagon of that orientation |
| Octagon | the point leaves an octagon |

The polygons are measured by their **apothem** — the distance from the centre
to the middle of a side — so a polygon and the circle of the same bailout
touch along the sides rather than at the corners. Orientation is part of the
shape: a triangle turned by 120° is the same triangle, so 0°, 90°, 180° and
−90° are four different ones, while a hexagon repeats every 60° and so has two,
and an octagon every 45°, which leaves it nothing to choose. The triangle and
the hexagon therefore have a submenu of their own for the orientation. Their normals are worked out once per frame, so the inner loop costs
three to eight multiply-adds rather than a sine and a cosine per side.

Against the circle, the diamond changes about 39% of pixels and the real axis
44%. The silhouette of the set does not move — that is decided by what never
escapes — so look at the bands.

Written to a saved position only when it is not the circle, so a position that
does not ask for a shape stays loadable by any earlier version.

## Zooming

**Selection zoom** (View menu) replaces the continuous zoom with a rectangle:
drag with the left button and release to go there, as Ultra Fractal does. Zoom
in only. The rectangle is forced to the aspect ratio of the window by growing
its short side, never by cropping the long one, and it works correctly with the
fractal rotated. Fast julia takes the drag while it is running, and gives it
back on the way out.

**Fixed steps**, also in the View menu:

| | |
| --- | --- |
| Zoom in 2× | `Ctrl` `+` |
| Zoom out 2× | `Ctrl` `-` |
| Zoom in 10× | `Ctrl` `Alt` `+` |
| Zoom out 10× | `Ctrl` `Alt` `-` |

They scale about the centre of the view, so they undo one another exactly.

## Antialiasing in linear light

Filters → Antialiasing in linear light, or `-linearaa`. Ordinary antialiasing
averages the sample colours as stored, which darkens edges because the stored
values are not proportional to light. This converts to linear light, averages
there, and converts back, so an edge keeps the brightness of what it is
between.

## Building

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Both binaries land in `bin/`. Options:

| | |
| --- | --- |
| `-DQUAD=OFF` | skip the second binary |
| `-DDEEPZOOM=ON` | make the main binary quad, and build no second one |
| `-DDIAG=ON` | also build `XaoS-diag`, which reports overlay and calculation ordering to `xaos-diag.txt` |
| `-DOPENGL=ON` | the OpenGL driver |

`cmake --build build --target deploy` fills `bin/` with the Qt and compiler
runtime so it runs on a machine without Qt.

## Tests

    ctest --test-dir build

293 tests: the formula parser at both precisions, the accuracy of `gamma`,
`lambertw`, `erf` and `randsc` against exact references where any exist, the
iteration loops against recorded checksums, saved positions round-tripping, the
overlay save/restore, and the help reference against the parser's own table.

One of them measures speed rather than answers. Timing in nanoseconds would
assert which machine the suite is running on, so it asserts a ratio instead:
what a degradation costs the noise functions, counted in multiplications. It
should be two, one per component, because the size is multiplied by the
degradation once per pass; raising it to the power of the pass instead reads
six or seven and climbs with the iteration limit. The test fails on the two
implementations this one replaced.
