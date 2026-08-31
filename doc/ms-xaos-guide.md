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
the formula being written. The list is checked against the parser's own table
by a test, so it cannot drift.

**Removed.** `powi`, `powdc` and `logcn` were second names for `pow` and
`logn`; `rad`, `deg` and `sign` were listed with no implementation behind them,
and a second `trunc` was shadowed by the working one. A formula using a removed
alias should use the name that remains — the function is identical.

**`erf(z)`** — the error function over the complex plane. Accurate to about
three ulp inside a bailout of two, where a fractal actually iterates.

**`randsc(seed; size; degradation)`** — coherent noise over the point, giving
blobs rather than per-pixel snow. `size` (default `1+i`) is the average width
of a blob along the real axis and its height along the imaginary one.
`degradation` (default `1+i`) shrinks them as the iteration proceeds: the size
in force is `size * degradation^n`, taken component by component, so `0.5+0.2i`
over `1+i` gives `1+i` on the first pass, then `0.5+0.2i`, then `0.25+0.04i`.
A zero in either component of either argument returns zero rather than dividing
by zero. Only the seed is required.

**`randscq(...)`** — the same field without the interpolation: a mosaic of flat
square cells instead of blobs. Same arguments, same meaning.

Both hash the position of the point and the iteration, never `z`, and never any
global state. So the same picture comes back on every redraw, at any thread
count, and in Mandelbrot or Julia mode alike. `randsc` is continuous, so the
two precisions agree to about 1e-19; `randscq` is a step function and they can
disagree on a hairline along the cell edges, which is inherent in asking for
hard edges.

A formula calling either turns boundary tracing off. That optimisation fills a
region it found one colour around without computing it — true of a fractal,
false of a noise field.

The old `rand` is unchanged and still depends on how many times it has been
called, so the same position does not redraw the same way. Prefer `randsc`.

### Getting a picture out of the noise

`randsc` returns a value in `[0, 1)`. Used alone the iteration never leaves a
bailout of 4, so nothing escapes and the image is flat. Multiply it, or let it
perturb an iteration that does escape:

    z+randsc(13;{0.35,0.35};{1,1})*0.25       blobs
    randsc(13;{0.15,0.15};{0.5,0.5})*2.5      brownian motion
    randscq(13;{0.15,0.15};{0.5,0.5})*2.5     scattered squares, shrinking
    z^2+c+randsc(13;{0.25,0.25};{1,1})*1.2    the set itself deformed

A degradation of 0.5 halves the blobs every pass, so after twenty iterations
they are a millionth of their size and below a pixel. For a slow fade over a
long run, use something near 1 — `0.97^60` is still 0.16.

An integer seed is exact in both builds; a fractional one is quantised to its
leading bits, which agree except for about one seed in fifty million.

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
| Triangle 0° / 90° / −90° | the point leaves a triangle of that orientation |
| Hexagon 0° / 90° | the point leaves a hexagon of that orientation |
| Octagon | the point leaves an octagon |

The polygons are measured by their **apothem** — the distance from the centre
to the middle of a side — so a polygon and the circle of the same bailout
touch along the sides rather than at the corners. Orientation is part of the
shape: a triangle turned by 120° is the same triangle, so 0°, 90° and −90° are
three different ones, while a hexagon repeats every 60° and an octagon every
45°. Their normals are worked out once per frame, so the inner loop costs
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

257 tests: the formula parser at both precisions, the accuracy of `gamma`,
`lambertw`, `erf` and `randsc` against exact references where any exist, the
iteration loops against recorded checksums, saved positions round-tripping, the
overlay save/restore, and the help reference against the parser's own table.
