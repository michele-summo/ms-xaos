# MS XaoS

A fork of [XaoS](https://github.com/xaos-project/XaoS) 4.3.3. Everything the
original does, it still does; this describes what has been added or changed,
and why. Version 1.7.2.

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

## Numbers in the dialogs

A field that asks for a number takes one written any way it can be read, 1e-18
as readily as 0.5. The single-field dialogs used a spin box, which keeps two
decimal places unless told otherwise, so 0.01 survived and 0.001 became
nothing — a Newton convergence of a millionth could not be entered at all.

What is shown back is the shortest text that reads as the very number held, so
a millionth of a millionth looks like 1e-18 rather than like
9.99999999999999999978e-19, which is what printing every digit the build
carries makes of a decimal a binary float cannot hold exactly. A coordinate
that needs all twenty-one digits still gets them: needing them is the same
thing as not reading back unchanged without them.

The exponent form is used from a thousand billion upwards and below a hundred
billionths, and plain decimal between — so twenty is twenty, a hundred
thousand million is written out, and 1e-8 is not. Left to itself the printing
would choose by comparing the exponent against the number of digits asked for,
which is no use when the digits are chosen for shortness: one digit is all it
takes to read back as twenty, and twenty came out as 2e+01.

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
A fourth tab lists the numbers that appear as arguments and mean something
particular — which shape a trap measures against, which mirror the
kaleidoscope folds with — rather than burying them in the description of the
function that takes them. A fifth shows a picture of every tiling `randsctile`
takes, since forty-five numbers with a line of words each say little about the
shapes.

Each row of the function list says what the call takes and in what order —
"a, b", or "seed, size, degradation" where the position means something. The
list and those counts are both checked against the parser’s own table by a
test, so a function that gains or loses an argument cannot go on being
described with the one it used to have. An argument that may be left out is
shown in brackets with the value it takes when it is: `a, [b=1], [c=1]`.

**Arguments are separated by a comma, and by nothing else.** A semicolon used
to serve as well: the parser turned it into a comma before reading anything, so
`LOGN(5;Z^2)+C` and `LOGN(5,Z^2)+C` were one formula, and the same went for the
two parts of a complex number. That is gone. A semicolon is now refused where it
stands — *Invalid operator: ;* — and the comma is the separator throughout, as
it is everywhere else in the syntax.

### The palette probe

**View → Palette probe (click a point)** turns the pointer into a probe: click
anywhere on the fractal and the message line says **where in the palette that
point was drawn from**, as a number rather than as a colour —
*Palette place 12.375*, and nothing else.

The palette editor lays **31 colours** down, and the palette is built by walking
from each to the next **eight cells at a time**. So the number is the cell over
eight: `12.000` is exactly the editor's twelfth colour, `12.375` three eighths
of the way from the twelfth to the thirteenth. A cell is an eighth of a colour,
finer than anyone reads off the editor.

* With the **custom palette** in use the number is brought round to **`0` to
  `31`**, which is the editor's own numbering.
* With any other palette it is **left where it falls and runs far past 31** — a
  palette of 65534 cells reaches 8191. The number out of range is what says the
  palette is not the one the editor draws.
* A point that took the **inside colour** reads *Palette place none*: that
  colour is one flat tone no cell of the palette stands for.

What it is for is less "what colour is this" than **which part of the palette
the picture actually uses**, and the answer is often sobering. With the plain
`iter` colouring at 60 iterations the cells run 1 to 59 — `0.125` to `7.375` —
so the picture never touches the editor's colours 8 to 30 and editing them
changes nothing. A colouring that spreads, `iter + real` say, walks the whole
ring instead, several times over.

It reads the point rather than the screen, so it is exact: the same colour can
sit in two cells, and a truecolour picture blends two neighbouring cells
together, so no colour names a cell of its own. Like the selection zoom, the
probe hands the pointer back to fast julia mode while that is running.

**Under the bar there is a line saying what the call you are in takes**, with
the argument you are writing in bold. Put the cursor inside `randsc(` and it
reads *randsc(**seed**, [size=1+i], [degradation=0.5+0.5i], [kaleidoscope=1],
[mode=0], [skew=0])*; move past a comma and the bold moves with you. It needs
no more of the formula than an open bracket with a name in front of it, so it
is there while the formula is half written, which is when it is wanted. Put the
cursor outside every call and the line is empty.

Six of the shipped positions were written the old way and have been rewritten
with commas: `circle`, `heart`, `helloween`, `pentafrac`, `warriormask` and
`burnship`, all under **File → Load position → examples/Malczak**. **A position
of your own that has a semicolon in its formula will not load** until the
formula is written with commas — the message names the character, and the fix
is to open the formula and replace them.

**Leaving a place empty.** An argument shown in brackets may be left out in the
middle of a call as well as at the end, by writing nothing between the two
separators: `julian(z, ,3)` gives the first and the third and lets the function
say what the second is. Spaces make no difference — `f(z, ,5)` and `f(z,,5)`
are the same call — and an argument that is not in brackets must still be
written, so `poly( ,1,2)` and `sin(,)` are refused rather than guessed at.

What an empty place means is the function's own business. Most take the default
they declare; a coefficient of `poly` left empty is a term that is not there,
so `poly(z,1, ,1)` is `z^2 + 1`; and a branch of `ifiter` or `ifiterl` left
empty repeats the one before it, which is how a branch is given more than one
pass in the cycle. `ifiter(f(z), , ,g(z), , , , )` runs `f` for three passes
and `g` for five, and costs no more to evaluate than writing them out would:
where each choice leads is settled once, when the formula is parsed.

A user formula that has not been written yet says `z^2+c`. It used to say the
burning ship, which is a fractal of its own and a puzzle to meet as a starting
point.

**`ifiterf(a, b)`** — evaluates `a` on every pass but the final one and `b`
on that. The final pass is the last the iteration limit allows: a formula has
no way of knowing which pass will be the one that escapes, that depending on
the value it has not produced yet. Only the chosen one is evaluated.

**`ifiterr(a, b, n)`** — evaluates `a` while the pass number is below `n`
and `b` from `n` onwards. Only the chosen one is evaluated, as with the other
two, though the threshold had to be taught to the parser first: an argument
may now be marked as read by the selector rather than chosen by it, and is
then evaluated before the choice is made.

**The variables.** `z` is the running value; `c` is the point, which is the
pixel in mandelbrot mode and the constant in julia mode; `x` is the plain
coordinate, which is the pixel in either mode and does not change from pass to
pass; `n` is the iteration number. In the *User initialization* `z` is where `z`
would have started had there been no initialization, which is the same value
`x` is, so an initialization of `z` alone says exactly what saying nothing
says.

That last part is new, and so is being able to say it: `z` was registered on
neither of the two parsers the initialization is run through, and `x` on the
wrong one, so an initialization naming either was refused by the thread that
had to compute it -- while the dialog, whose parser has both, accepted it and
gave no sign that it was then thrown away. `n` was left holding whatever the
pixel before it had left behind, which is to say an initialization reading it
read the order the pixels happened to be computed in.

**`p1`, `p2` ... `p9999`** — the value `z` had on an earlier pass: `p1` on the
pass before this one, `p2` the one before that, `p9999` nine thousand nine
hundred and ninety-nine passes back. `p` is another name for `p1`. Before there
have been that many passes they stand at whatever the history starts from,
which *Fractal → Set p values on first iteration* decides: the point being
iterated when it is on, zero when it is off.

They used to stop at `p6` and were kept by shifting a six-place array along by
one every pass — a copy per place per pass, paid whether or not the formula
named any of them. That is why they stopped at six: at nine thousand it would
have been 320 KB of copying per pass, measured at twenty-four times the whole
cost of iterating `z^2+c`. The history is a ring now, written once per pass,
and only the places a formula actually names are read out of it, so `p9999`
costs what `p1` costs and a formula that names none costs nothing at all.
The reference window lists them among the variables.

**Suffixes on a variable.** `z`, `c`, `x`, `p` and `p1` to `p9999` take
suffixes that stand for the calls made on them most often:

| suffix | stands for |
| --- | --- |
| `_b` | `bship(…)` |
| `_bi` | `bshipi(…)` |
| `_br` | `bshipr(…)` |
| `_pM` | `parchment(…, M)` |
| `_paM` | `parchmenta(…, M)` |

They are read from left to right, each wrapping what the ones before it made:
`z_p3` is `parchment(z,3)`, `c_b_p2` is `parchment(bship(c),2)` and `p12_p2_p3`
is `parchment(parchment(p12,2),3)`. *M* is part of the name, so it is written in
figures — a whole number, one or more; `z_p0` is refused with a message that
says so. Not on `n`, which is a count and has no second component for any of
them to work on.

They are spelled out as the calls before the formula is read, so a suffixed
variable draws exactly what the calls written out draw and costs what they
cost; the formula is saved as it was written. Anything that is not one of the
five is left alone and refused as the unknown name it is.

**Removed.** `powi`, `powdc` and `logcn` were second names for `pow` and
`logn`; `rad`, `deg` and `sign` were listed with no implementation behind them,
and a second `trunc` was shadowed by the working one. A formula using a removed
alias should use the name that remains — the function is identical.

**`erf(z)`** — the error function over the complex plane. Accurate to about
three ulp inside a bailout of two, where a fractal actually iterates.

**`poly(z, k1, k2, ..., km)`** — a polynomial in `z`:

    k1*z^(m-1) + k2*z^(m-2) + ... + k(m-1)*z + km

The first coefficient written multiplies the highest power and the last stands
alone, so the call reads in the order one says the polynomial. Worked out by
Horner's rule, which is m−1 multiplications rather than the m(m−1)/2 that
raising each power separately would take, and the more accurate of the two into
the bargain.

    poly(z,1,0,0)+c            the Mandelbrot, written out
    poly(z,1,0,0,{0.7,0.2})    z^3 + 0.7+0.2i

**`randsc(seed, size, degradation, kaleidoscope, mode, skew, skew_mode, selfsim)`** — coherent
noise over the point, giving blobs rather than per-pixel snow. `size`
(default `1+i`) is the average width of a blob along the real axis and its
height along the imaginary one.
`degradation` (default `0.5+0.5i`, halving them each pass) shrinks them as
the iteration proceeds: the
size is multiplied by it at every pass, component by component, so `0.5+0.2i`
over `1+i` gives `1+i` on the first pass, then `0.5+0.2i`, then `0.25+0.04i`.
A zero in either component of either argument returns zero rather than dividing
by zero. Then come the two kaleidoscope arguments, `skew` (default `0`),
`skew_mode` and `selfsim`, each described below. Only the seed is required.

### The skew, and why these fields are hard to colour

The engine colours with the two components of the orbit: `real` and `zmag` read
one, `imag`, `angle` and `real / imag` the other. These fields hand back one
number for a whole cell, on the real axis, so every mode draws **one tone a
cell** — `zmag` and `iter + real` look like the value truncated, and the other
three have nothing at all to read. Measured over nine hundred pixels: 14 to 126
values for the first two, and exactly **one** for the rest.

Varying the value inside a cell cannot be done with one real number. The same
number decides whether the point leaves, so a value that moves across a cell
takes half the cell out and leaves the other half in, and the cell comes out
**cut in two by a straight line**. That was tried twice and cut the mosaics both
times.

`skew` does it with the second component instead, and it **turns** the value
rather than scaling it. How far it turns follows **how far out of the middle of
its cell** the point stands: with `out` running from nought in the middle to one
at the edge, the value is turned by twice the arc tangent of
`skew_re*out + skew_im`.

`out` is measured **in the cell's own geometry**, so what the colour draws
inside a cell is the shape the field is cut into, not a set of lines laid across
it — nought in the middle of a cell, one along the whole of its edge:

| field | what the colour draws inside a cell |
| --- | --- |
| `randsc` | its own blobs — the turn follows the level, so the contours are the field's |
| `randscq` | squares, about the middle of the square |
| `randsch` | hexagons, about the middle of the hexagon |
| `randsct` | triangles, about the middle of the triangle |
| `randscp` | its own polygon, shrunk in step by step — no two cells the same shape |

The **real part** of the skew is that gradient across a cell. The **imaginary
part** is a flat turn the whole cell shares: it shifts a cell's colour without
drawing anything inside it, so `0.4+0.2i` gives the shape's contours and a
per-cell shift together.

### What the skew draws: `skew_mode`

The seventh argument says **what varies inside the cell**. The numbers **add
together** rather than choosing one another out, so a call may ask for two of
them and get both.

| | |
| --- | --- |
| **1** (and `0`) | **the shape** — the cell's own outline, shrunk step by step. What a call that names no mode gets, and what is used when neither 1 nor 2 is named |
| **2** | **the rosette** — the angle round the middle of the cell, folded into as many turns as the kaleidoscope has wedges, so each cell carries the picture's own symmetry. Plain spokes when nothing is folded |
| **4** | **across the wedge** — a turn that follows where the point stands in its wedge of the kaleidoscope, the same in every wedge |
| **8** | **radial** — the modulus moves as well as the angle |

`3` is the shape and the rosette at once, which spirals. `11` is those two with
the modulus moving as well.

**`4` follows the kaleidoscope's own geometry.** It turns the value by
nothing on the edges of each wedge and by the whole of the turn the imaginary
part gives on its own — twice its arc tangent, 127 degrees for `2` — down the
line through its middle, and by the turn of the imaginary part times how far
across in between. Every wedge gets the same, so a kaleidoscope of six is still
a kaleidoscope of six, mirrors and all, and the turn is continuous across the
edges and the middles, which are where the fold mirrors. It does nothing when
nothing is folded, or when the imaginary part of the skew is nought. No
trigonometry, and it costs nothing at 64 bits of mantissa and a tenth of a call
at 113.

It was first **a different turn in each wedge** — the imaginary part times the
wedge's number plus one — meant to make the copies stop being identical. That
is the one thing a kaleidoscope's copies cannot be: a six-fold picture came out
as six unrelated slices with a straight cut at every join. The multiple also
ran up against half a turn, so four slices of six sat nearly on the real axis,
and there `real / imag` drew them as snow however `selfsim` averaged the noise.

What is left of that is a property of `real / imag`, which divides by the
imaginary part of *z*: wherever a turn carries it across nought the colouring
jumps from one end of the palette to the other, and the last scraps of the fine
passes decide the side pixel by pixel. A turn that sweeps across the wedge
crosses there for most settings. On a six-fold `randsct` coloured by
`real / imag`, with `selfsim` at one, an imaginary part of `0.25` left no grain
at all and `2` left bands of it beside the edges; a colouring that does not
divide by one component has none.

**`8` is the one that answers `zmag`.** A turn leaves the modulus where it is,
and `zmag` and the bailout read the modulus and nothing else — so no other mode
can reach them. The price is the one the skew was designed to avoid: the escape
moves with the modulus, so the figure is cut wherever the bailout falls inside a
cell. Pair it with a circular bailout and a colouring that reads the components
apart, or accept the cutting on purpose.

A skew of nought leaves every mode doing nothing, so `skew_mode` alone changes
no picture.

A straight ramp was written first — `skew_re*du + skew_im*dv` across the cell —
and it is what this replaces. It coloured, but it drew the same diagonal,
vertical or horizontal bands across every cell alike whatever the field was cut
into, and the bands ruined the picture.

At **`0`**, which is what a call that does not name it gets, the turn is
nothing, the imaginary part stays at nought, and every value is the number it
always was **to the bit** — the golden checksums from before the argument
existed pass unchanged, and a saved position renders identically. The `talc`
position, three noise calls and a `c` over 345600 values, comes out to the same
signature as the build before the argument.

**How much colour it gives.** The engine's index is
`(iter + quantity) * speed + shift`, so a spread of one in the quantity is one
band's worth at a speed of one. The spread of `imag` inside a cell:

| skew | 0.02 | 0.05 | 0.1 | 0.2 | 0.4 | 0.6 | 1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `randsc`, `randsct` | 0.08 | 0.21 | 0.43 | 0.83 | 1.5 | 1.9 | 2.1 |
| the three mosaics | 0.06 | 0.14 | 0.28 | 0.53 | 0.96 | 1.2 | 1.4 |

So **`0.3` to `0.6` is where a band appears at a speed of one**, which is the
thing these fields were missing. Counting distinct values would say `0.01` is
enough — nine hundred of them — but a thousand values inside a hundredth of a
band are all one colour.

**What it costs.** Scaling moved the modulus of the value, and a round bailout
looks at exactly that: a cell whose level sat near where the bailout falls came
out cut in two, and the mosaics lost their shapes. A turn leaves the modulus
where it is, so **under a circular bailout the skew is free at any strength** —
over 90000 pixels, at 0.05, 0.3 and 1, and over all five fields: not one pixel
leaves differently.

A bailout polygon is another matter, because it does not read the modulus but
the components, and a turn walks the value round a circle that can cross a side.
How far it can cross by is how far the shape's corners stand out past its sides,
and the numbers follow that — the fraction of the picture whose escape changes,
over the five fields:

| bailout shape | corners over apothem | skew 0.01 | skew 0.05 | skew 0.3 | skew 1 |
| --- | --- | --- | --- | --- | --- |
| circle | 1.00 | **0%** | **0%** | **0%** | **0%** |
| hexagon | 1.15 | **0%** | 0.0–0.8% | 0.7–2.2% | 0.4–2.4% |
| square | 1.41 | **0%** | 0.0–1.0% | 0.8–2.8% | 1.2–12.8% |
| triangle | 2.00 | 0.05–1.3% | 0.7–5.0% | 2.6–6.4% | 7.4–21.6% |

The **number** the bailout is set to moves it as well, and for the same reason:
what can change side is a value whose modulus falls between the polygon's
apothem and its corners, and where that ring sits among the cells is what the
number decides.

So the two wants pull against each other under a polygon, and there is no
arranging otherwise: the escape reads the two components of the value and so
does the colouring, so anything that gives the colour something to read is
something the escape can see. **Pair the skew with a circular bailout** and it
costs nothing at any strength. Under a polygon, either take the few per cent or
keep the skew at `0.05` and raise the colour speed instead.

What a turn cannot touch is **`zmag`**, which reads the modulus. Colour with
something that reads the two components apart: `iter + real`, `iter + imag`,
`angle`, `real / imag` outside, and `real` or `real / imag` in the incolouring.

Two things follow. The value is complex while the skew is not nought, so
`randsc(7,,,,,0.4)*z` turns `z` as well as scaling it. And the kaleidoscope
folds either way: measured over five fields, two, three, five and six wedges and
both mirrors, a turn of one wedge leaves every value where it was.

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

**`randsctile(tiling, seed, ...)`** — the same field over any of **forty-five
tilings**, which the first argument chooses; everything after it is what the
family takes, one place further along, so `selfsim` is the ninth. Every tiling
is scaled to a tile of unit area on average, so the first argument changes the
shape of the cells and not their scale. A number outside 1 to 45 draws nought.

| | |
| --- | --- |
| 1–3 | the regular tilings: squares, triangles, hexagons |
| 4–11 | the eight Archimedean ones, regular polygons with every corner alike: 4.8.8 (octagons and squares), 3.6.3.6, 3.4.6.4, 3.12.12, 4.6.12, 3.3.3.4.4, 3.3.4.3.4 (snub square), 3.3.3.3.6 (snub hexagonal) |
| 12–19 | their duals, which are not regular: tetrakis square, rhombille (the cubes), deltoidal trihexagonal (kites), triakis triangular, kisrhombille, and the prismatic, Cairo and floret pentagons |
| 20–32 | bricks, Flemish bond, herringbone, basketweave, Pythagorean (squares of two sizes), chevrons, squares and rhombi, houses (a pentagon, up and hanging by turns), rows of squares and triangles in two rhythms, hexagons among triangles, Greek crosses, T tetrominoes |
| 33–37 | Islamic stars: eight-pointed with crosses, six-pointed with hexagons, eight-pointed from 4.8.8, twelve-pointed from 3.12.12 and from 4.6.12 |
| 38 | Voronoi cells, irregular — `randscp`'s, with a value of its own |
| 39–43 | tilings that never repeat, from de Bruijn's multigrids: Penrose's rhombs, Penrose's kites and darts, Ammann–Beenker's squares and rhombi, and rhombi in twelve and in seven directions |
| 44–45 | by substitution: the pinwheel, whose triangles face every way there is, and the chair, Ls cut into Ls |

The **Tilings tab** of Help → User formula reference has a picture of each,
captioned as the Values tab lists it, and a double click on one copies the
start of its call — `randsctile(18, ` — to the clipboard, ready for the seed.
The pictures are drawn by `randsctile` itself, eight units square round the
origin, by `tools/randsctile-thumbnails.cpp`, and kept in the binary rather
than drawn each time the tab is opened. Kept pictures can go stale, so each is
stored with a fingerprint of what the function drew, and a test draws the same
points again and fails, saying which tiling, until
`cmake --build <build directory> --target randsctile-thumbnails` has drawn it
anew.

The periodic ones are **tables**, written by `tools/randsctile-tables.py`, which
builds each tiling from its geometry and checks it — twenty thousand points at
random, each of which must fall in exactly one tile — before writing it. The
Islamic stars are made from the Archimedean tilings: every polygon of eight
sides or more becomes a star with its points on the polygon's corners, and the
thin triangle between each notch of the star and the side goes to the tile
across that side, so two stars sharing a side make a rhombus between them and a
small polygon takes a notch on each side it shares with a star. Tiles may be
concave — stars, crosses, the T — so a point is placed in its tile by the
crossing number rather than side by side. The point is split into whole periods
in the working precision, as `randscq` splits it into cells, so a periodic
tiling is exact however far out.

The ones that never repeat have **no table** and are worked out afresh at every
point. Penrose's rhombs, Ammann–Beenker and the twelve- and sevenfold rhombi
come from **de Bruijn's multigrid**: families of parallel lines in a second
plane, a tile wherever two cross, named by the two lines — four whole numbers,
exact whoever asks. The point is placed by trying the crossings next to where
it lands in that plane, the families whose lines pass nearest first, which
finds it in four to eight tries on average. The kites and darts come from
Robinson's triangles cut level by level out of a wheel of ten round the origin,
the pinwheel from a right triangle cut into five, the chair from an L cut into
four; at each level a line or two decides which child the point is in.

**What it costs**, measured per call against `randscq`'s 220 ns on this
machine: the periodic tilings 255 to 380 ns, the ones that never repeat 410 to
680. They are worked in `double` in both builds, the place within a period
included — in long double the ones that never repeat took three times as
long — and at that the two builds drew all forty-five tilings alike to the
pixel, forty thousand points of each. Past 2³² cells from the origin, where
`double` would start to lose their tiles, the tilings that never repeat go flat
as the family does past its grid; the periodic ones hold out as far as the rest
of the family does.

The skew measures a tile as `randscp` measures its polygon, against the radius
of the largest circle the tile holds, so the contours are the tile shrunk,
concave ones included; the rosette turns about the tile's middle.

    z+randsctile(18,13,{0.3,0.3},{1,1})*0.25        Cairo pentagons
    z+randsctile(36,13,{0.3,0.3},{1,1})*0.25        twelve-pointed stars
    z+randsctile(39,13,{0.3,0.3},{1,1},5)*0.25      Penrose, folded five ways

**`fbm(value, seed, [intensity=4], [frequency=8], [octaves=4],
[roughness=0.5])`** — a fractional Brownian motion: the same noise the family
above is built from, summed in octaves, each at twice the frequency of the one
before and keeping `roughness` of its height. What it draws is **wear** rather
than a pattern — stains and dents — because no octave is large enough to see on
its own and none is small enough to disappear.

Where `randsc` reads the position and can only read the position, **this reads
whatever you write in front of it**. That is the whole point of having it as a
function:

    fbm(z,7)          moves with the orbit
    fbm(x,7)          stands still on the plane
    fbm(z*3+c,7)      whatever that is
    z^2+c+fbm(z,7,0.1)*i    the set with a rough edge

The value and the seed are required; the rest have defaults. The motion runs
from **nought to `intensity` and never below**, so adding it to something cannot
pull that under nought. The seed is read the way `randsc` reads one, so the same
number means the same field in both. `octaves` is held between one and
twenty-four, and it is worth raising only when `roughness` is: at the plain half
the eighth octave carries a two hundred and fiftieth of the whole.

**It has no kaleidoscope of its own**, where the `randsc` family does, because
one is already there to compose with: `parchment(a,b)` folds the angle of `a`
into `|b|` sectors and `parchmenta` mirrors the alternate halves, so
`fbm(parchmenta(z,6),7)` says it, and says it where the reader of the formula
can see it.

The colouring modes under **Fractal → Outcoloring mode → Other coloring mode →
Fractional Brownian Motion** apply the same motion to a colour. This one puts it
where a formula can use it.

| | |
| --- | --- |
| `randsc` | soft, curved |
| `randscq` | hard, regular, square |
| `randscp` | hard, irregular |
| `randsch` | hard, regular, hexagonal |
| `randsct` | hard, regular, triangular |
| `randsctile` | hard, any of forty-five tilings |

All six lay one cell over each unit square of the size in force, so `size`
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
between the six functions, so a formula subtracting one from another does not
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

    randsc(13,{0.6,0.6},{1,1},6,0)     six wedges, each a mirror of itself
    randsc(13,{0.6,0.6},{1,1},3,1)     three wedges, folded the other way

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

    z+randsc(13,{0.35,0.35},{1,1})*0.25       blobs
    randsc(13,{0.15,0.15},{0.5,0.5})*2.5      brownian motion
    randscq(13,{0.15,0.15},{0.5,0.5})*2.5     scattered squares, shrinking
    z+randscp(13,{0.3,0.3},{1,1})*0.25        irregular polygons
    z+randsch(13,{0.3,0.3},{1,1})*0.25        a honeycomb
    z+randsct(13,{0.3,0.3},{1,1})*0.25        triangles
    z^2+c+randsc(13,{0.25,0.25},{1,1})*1.2    the set itself deformed

A degradation of 0.5 halves the blobs every pass, so after twenty iterations
they are a millionth of their size and below a pixel. For a slow fade over a
long run, use something near 1 — `0.97^60` is still 0.16.

An integer seed is exact in both builds; a fractional one is quantised to its
leading bits, which agree except for about one seed in fifty million.

### Every pass so far: `selfsim`

Each pass of these functions is a field of its own — the iteration goes into
the hash — with cells the degradation times those of the pass before. A formula
that calls one on every pass therefore reads a new and finer field each time,
and a colouring mode that reads the last pass reads the last field alone. At a
degradation of 0.75 the cells of a view five wide are below a pixel by the
sixteenth pass, and past that the inside is snow: a hundred passes of
`randsct(672,,{0.75,0.75},,,{2,2},6)+c` in julia mode put 70% of neighbouring
inside pixels a sixteenth of the palette or more apart.

The eighth argument hands back instead **the passes so far, averaged**, standing
to one another as the octaves of a fractional Brownian motion do: pass *n*
weighs d^(nH), where *d* is the degradation — the geometric mean of its two
components taken without their signs, so `0.75+0.75i` is 0.75 — and *H* is the
argument. The weights are divided by their total, so the answer stays among the
values averaged: in [0, 1] without a skew.

| selfsim | |
| --- | --- |
| left out | off — the one pass, to the bit. What a call that names nothing gets, and a place left empty |
| `0` | the plain average: every pass weighed alike, and what every *H* near nought gives — not off |
| `1` | the plain motion: each pass weighs *d* times the one before |
| `0.5` | rougher: the fine passes keep more of their weight, so more detail shows |
| larger | smoother: the first few passes take nearly all of it |
| negative | the fine passes outweigh the coarse ones, which is the snow back again |
| complex, `{1,2}` | the real part weighs as above, the imaginary part turns the passes — see below |

    randsct(672,,{0.75,0.75},,,{2,2},6,1)+c     the same, self-similar

On that formula the snow goes: neighbouring inside pixels come out 1.8 palette
cells apart on average rather than 35, and 0.6% of them more than sixteen apart
rather than 70%. What is left is the triangles at every scale at once, the large
ones carrying the small.

* **`selfsim` may be complex.** Then d^(nH) is complex too: d^(n·Hr) in size
  and a turn of n·Hi·ln d, so every pass is turned that much further than the
  one before and then averaged as a real *H* averages it — a spiral of
  octaves. At 0.75 an imaginary part of one is some sixteen degrees a pass.
  The weights are divided by **the sum of their sizes**, not by their complex
  sum: that keeps the answer within the largest of the values averaged. The
  complex sum would come near nought at some pass and blow every value up at
  that pass for every pixel at once, the weights having nothing to do with the
  pixel. The turns give the value an imaginary part even without a skew, so
  `imag` and `angle` have something to read. A degradation of one has a
  logarithm of nought and turns nothing; `{1,0}` is `1` to the bit; and
  `{0,2}` is every pass weighed alike, and turned.
* **Nought is not off.** d^(n·0) is one for every pass, so `0` is the plain
  average, and so is every *H* near it — from either side, and complex too.
  Nought was the switch at first, which made it the one value the curve does
  not pass through: `0.000001` was the plain average and `0` the last pass
  alone, a quarter of the range apart. The last pass alone is the far end of
  the curve, *H* running to minus infinity. To have no average, leave the
  argument out, or its place empty.
* **It is the skewed value that is averaged**, so the skew and every bit of
  `skew_mode` go on doing what they did, and the answer is complex when the
  skew makes it so.
* **An average draws the values together.** At 0.75 and *H* of one the spread
  is about two fifths of what one pass has, so raise the colour speed to get the
  contrast back.
* **A degradation of one** weighs every pass alike and gives their plain
  average, which settles toward a flat middle as the passes add up. That is
  what the definition says, and it is done as said.
* **Passes a call did not see are worked out when it is asked.** A call on a
  branch the formula takes only from the seventh pass still gets the first
  seven, and the answer at a pass is the same whichever route reached it. A new
  pixel is told by its position as well as by its pass: the pass alone, which is
  what `trap` goes by, would let such a call carry on the average the pixel
  before it left.
* **Nothing overflows.** The average is kept as one, each pass taking a share
  of it that follows from the share of the pass before, rather than as a sum and
  a total — so a negative *H*, or a degradation above one, whose weights would
  run past what the type holds, has nowhere to do it.
* **It costs** about a tenth more per call at 64 bits of mantissa and a seventh
  at 113. Off, nothing that can be measured.

### Figures instead of noise

`sierpinskyt`, `sierpinskyc` and `snowflake` stand where the noise stands, and
nothing about them is random: the figure is where it looks like it is, and
mandelbrot mode and julia mode draw the same picture, a shape in the plane being
in the same place either way.

**They are fractals of their own, not fields to multiply into one.** Written
alone — `sierpinskyt()` and nothing else — each draws its figure, the way
Fractal → More Formulae draws its **Sierpinski**, **Sierpinski Carpet** and
**Koch Snowflake**.

All three do it by **carrying the point to its parent**. Every part of one of
these figures has one: a hole in a gasket sits inside a bigger hole one level up,
a hole in a carpet inside the cell that was cut the same way, a triangle of a
snowflake on a side of the hexagon at its middle or on the free edge of another
triangle. The step onto the parent is a motion of the plane — doubling away from
the nearest corner for the gasket, blowing a cell up by the number of cells for
the carpet. The topmost part has no parent inside the figure, so its step carries
it out of the bailout, and a part *n* levels down takes *n* steps to get there:
**the pass a point leaves on is the level it stands at**, and the iteration count
is the picture.

What comes back on the pass a point leaves on matters as much as the pass
itself: **every outside colouring mode but the iteration count reads it**. So the
topmost part of each figure is thrown out rather than simply declared gone — the
gasket doubles its middle hole away from a corner and it lands outside the
triangle, the carpet scales its middle cell about the cell beside it, the
snowflake throws its middle hexagon out the way that sector faces. Each carries
where it came from with it, so `real`, `imag`, `angle`, smooth colouring and the
rest all have something to work on. Only a point that was never in the figure is
handed a number standing in for "gone".

A snowflake is read **from a hexagon out**, which is the thing to know about it.
It comes apart exactly into a regular hexagon at the middle, six triangles of
that hexagon's own side standing on its six sides, twelve of a third that on
their free edges, forty-eight of a ninth on theirs, and so on for ever. Three of
the six are the corners of the triangle the figure grew from and three are the
first bumps put on its edges, and **nothing tells them apart**: from the hexagon
out a snowflake has six-fold symmetry.

Read from the first triangle instead, five eighths of the figure is level one —
one flat triangle filling the picture, with the snowflake only in the fringe
around it, and nothing for the outside colour to say. Read from the hexagon, the
first two levels are five twelfths each and the rest come down evenly. A
triangle standing on the hexagon steps onto the largest triangle the hexagon
holds — corner on corner, three times the area — so it is blown up by the square
root of three and turned a twelfth of a turn; one standing on the free edge of
another triangle is blown up by three and turned to face the way that edge
faces.

**Every step expands**, and that is a rule rather than an accident. The gasket
doubles, the carpet blows a cell up by the number of cells, the snowflake by
three or by the root of three. A triangle standing on the hexagon could instead
be folded flat across the side it stands on, which lands it on one of the six the
hexagon is made of and is tidier; but a fold is an isometry, and an isometry has
no sensitivity to where the point started. Written inside a larger formula the
figure would then hand that formula a rigid picture of itself over five twelfths
of its area, and the formula would draw a flat region there. A test asserts the
three expansions.

The point really travels, and that is what makes these usable rather than merely
correct: `z` along the way is a point of the plane like any other, so the
colouring modes that read it, smooth colouring, and writing a figure inside a
larger formula all mean something.

A gasket and a carpet fill their shape, so every point in one has a level and
leaves on it. A snowflake does not fill its hexagon: it leaves six corners of
**ground**, and that ground is no part of the figure and has no level. It is
**turned half about**, which lands ground on ground, so it never leaves and comes
out in the **inside colour** — the empty space stays empty instead of taking a
band of the outside colour off the figure.

Turned, and not handed back where it stood, because **a fixed point is a stopped
orbit and not a bounded one**. Standing still, the loop ran every pass to the
limit with nothing changing; the incolouring modes that compare a pass with the
one before had nothing to compare; and `snowflake()` written inside a larger
formula handed that formula back the point it was already holding, where
`sierpinskyt()` and `sierpinskyc()` move every point they are given. Half about
is the choice over a sixth of a turn, which the figure's symmetry would also
allow: both land ground on ground, but half about keeps the point inside a square
and a circle as well as inside the hexagon. It leaves the distance from the
centre alone, so the incolouring modes that read that are untouched, and the
picture `snowflake()` draws is unchanged to the pixel.

That is the one place these figures cost something. A ground pixel is asked the
membership question once a pass for as long as it is looked at, where the other
two figures and the body of the snowflake leave on the pass their level names.
Drawn straight through the iteration loop, a whole snowflake costs about 2.4
times what it cost when the ground was given a band of its own. Boundary tracing
is left on for these formulas — their bands really are solid, which is the
condition for it — and the ground is one solid region, so what reaches the
screen is less than that number suggests.

Written inside a larger formula — `snowflake()^-2+snowflake()` and the like — the
three do **not** behave alike, and the reason is worth stating plainly because no
amount of work on the figure will change it.

A gasket is a **dust**: it has no area, its points never leave, and every one of
them is doubled every pass. A formula built on it has something chaotic under it
everywhere, at every depth. A snowflake is **solid**: ten twelfths of it is the
hexagon and the six triangles, which leave on the first and second pass, so a
formula built on it runs out of figure almost at once and draws whatever it draws
on its own. The detail is along the Koch boundary, where the levels run deep, and
the flat regions are the interiors, where they do not.

**What a pass knows is what that pass paid for.** The snowflake decides whether a
point is ground by walking down the Koch curve, and that walk is allowed **one
level per pass**. It used to run all twenty-four of its levels in a single
evaluation, so a formula written around the figure was handed the whole boundary
on its first pass and drew it in full however few passes it was given — at
`Iterations: 2` you still got every bump of the Koch curve, which no other
formula in XaoS does and which is not what an iteration count is for. Running out
of levels is not the same as being outside: a point out of the frame or above the
hull is ground and known to be ground now, while an undecided one is carried on
as part of the figure and asked again next pass with a level more to spend. The
picture `snowflake()` itself draws is unchanged — a point in the figure never
needed the deep walk, only the assurance that it was not ground.

What is left is the honest difference between the two figures, and no work on
either will close it. A gasket is a **dust**: no area, its points never leave,
every one of them doubled every pass, so a formula built on it has something
chaotic under it at every depth. A snowflake is **solid**, and ten twelfths of it
is the hexagon and the six triangles, which leave on the first and second pass.
Measured at one zoom — the fraction of pixels standing on a band edge —
`sierpinskyt()^-2+sierpinskyt()` reaches 17 per cent by twelve passes and
`snowflake()^-2+snowflake()` 2.5; both grow with the pass count, which is the
point, and where they stop is what the figures are.

`radius` means what `bailout` means and is read the same way: **as the square of
the distance**. What it draws is the shape a bailout of that number draws,
**inscribed in it exactly** — set the bailout to the matching shape and the same
number and the figure and the bailout lie over one another:

| | |
| --- | --- |
| `sierpinskyt` | the triangle of bailout **triangle −90°**, corner for corner |
| `sierpinskyc` | the square of bailout **square**, side for side |
| `snowflake` | the hexagon of bailout **hexagon 0°**, its six points on the hexagon's six corners |

The thing to know is that a bailout polygon stands its **sides** the square root
of the bailout from the centre — that is its apothem, not its circumradius — so
its corners are further out than the number says: twice as far for a triangle,
and by a seventh for a hexagon. `sierpinskyt(4)` therefore reaches 4 to its
corners while `sierpinskyc(4)` reaches 2 to its sides, and both are right.

Set the bailout larger than the radius and the figure no longer fills the shape:
there is a **margin** between the two, and set it smaller and the figure is cut
off. Both are what was asked for and neither is corrected. But the margin is no
part of the figure and has no level, so it is not given one: it is **turned by
the figure's own symmetry** — a third of a turn for a gasket, a quarter for a
carpet, half about for a snowflake — which keeps it outside the figure and inside
the bailout, so it never leaves and is drawn in the **inside colour** with the
incolouring modes working on it. It used to leave on the first pass and take the
lowest band of the outside colour, which drew the whole margin in one flat tone
and took a band off the figure.

A turn and not a mirror, though a mirror would serve as well for keeping it
where it is: what a figure hands back has to turn with the point, so that turning
the picture turns what is drawn on it, and a mirror turns it the other way. For
the same reason the middle of a carpet — the block of cells the border ring
encloses, which has no parent and so has to be thrown out — is thrown the way its
own **quarter** of the square faces, the quarters cut by the diagonals. Thrown
always the same way, as it was, it drew that middle as a ramp from one side to
the other with nothing across it; the square has four-fold symmetry and so should
what is drawn on it. A test asserts that each figure answers a turn of its own
order with the same turn.

**`sierpinskyt([radius=4], [kaleidoscope=1], [mode=0])`** — the Sierpinski
gasket, in that triangle.

**`sierpinskyc([radius=4], [squares=3], [kaleidoscope=1], [mode=0])`** — the
Sierpinski carpet, in that square. The square is cut into `squares` by
`squares`, the ring of cells along the border is kept, everything that ring
encloses is thrown away, and the same is done to each cell that was kept. So
the picture is **one square in the middle and 4×`squares`−4 around it**: eight
at three, which is the carpet as it is usually drawn, twelve at four, sixteen
at five. Two is the one number with no ring to speak of, and there the far
corner goes instead, which is a gasket again — a square cut in four with one
corner taken away is what a gasket is.

**`snowflake([radius=4], [kaleidoscope=1], [mode=0])`** — the Koch snowflake,
its points on that hexagon's corners, banded by generation from its middle
out: the hexagon on the first
pass, the six triangles on its sides on the second, the twelve on their free
edges on the third, the forty-eight on the fourth. The six corners of ground the
figure does not cover are not banded at all — they never leave, and are drawn in
the inside colour.

Every argument has a default, so `snowflake()` is a call, and so is
`sierpinskyc( ,5)` — a lacier carpet at the default size.

**The kaleidoscope is the noise family's**, written last on each of the three
and meaning the same: how many wedges the plane is cut into, and which mirror
folds them. `1` — what a call that says nothing gets — folds nothing and leaves
the figure the number it always was, to the bit.

What it folds is different, though, and has to be. The noise folds the
**position**, which stands still for the whole orbit; a figure has no use for
the position and reads **z**, the point it is carrying down towards the part it
stands in, so z is what is folded. The picture comes out the same way — the
figure drawn in one wedge and repeated round the origin — and it holds all the
way down, every step being taken on the folded point and folded again on the
pass after. Inside the half wedge the fold leaves alone, the figure is the
figure.

    sierpinskyt(4,6,0)      the gasket six times round the origin
    sierpinskyc(4,3,5,1)    the carpet in five wedges, the near mirror
    snowflake(4,5)          five wedges of snowflake

A snowflake already has the six-fold symmetry with the mirror down each
bisector, so folding one into six wedges that way is the identity on the
picture: `snowflake(4,6,0)` comes out **pixel for pixel** the same as
`snowflake(4)` over 90000 pixels. Ask for five wedges, or four, to see it fold.

One thing to know about the carpet: what it cuts away it throws one of four
ways, the quarters cut by the diagonals of the square, and on a diagonal itself
the choice is a tie. A fold can land a whole ray of the plane on that diagonal,
and there the two sides of the tie are thrown opposite ways — a hairline along
the diagonal, which is the step the carpet has always had rather than anything
the fold does. Measured over ten thousand points at six wedges: every
disagreement sat on a diagonal, and the folded points themselves agreed to a
part in 10^16.

    sierpinskyt()        with Fractal -> Bailout shape -> triangle -90
    sierpinskyc()        with bailout shape square
    sierpinskyc( ,5)     a lacier carpet
    snowflake()          with bailout shape hexagon 0

None of them costs more than the noise beside them, pass for pass: measured
against `randsc`, the gasket 0.44, the carpet 0.54 and the snowflake 0.50 of it
at 64 bits of mantissa, and 0.50, 0.87 and 0.65 at 113. Reading a snowflake from
its middle made it cheaper as well as better to look at: six sevenths of the
figure is the hexagon or a triangle standing on it, and those are answered by
four multiplications and a comparison without walking the curve at all. What a snowflake costs
over a whole picture is the paragraph above, and is about how many passes the
ground takes rather than about what a pass costs. That was not free, and at 113 bits it
was very nearly lost — a square root and a division are both software there, and
either costs about what a whole figure costs. So the reciprocal of the square
root of the radius is kept on the call site rather than taken again every pass,
and the root itself comes back by multiplying the radius into that reciprocal
rather than by dividing; the gasket compares its three barycentric weights
scaled by a positive number, which changes neither which is largest nor which is
negative; the carpet counts in cells rather than in the plane, where the
multiplication by `squares` and the division by it cancel; and the snowflake
takes the child to walk onto from the edge whose bump it was found under, which
it had to find anyway. It walks the Koch curve in `double` besides, since the
figure stands at a fixed size and 52 bits of a fraction of one edge is more than
the 24 levels it draws can use. A test asserts the ratios, so a rewrite that
loses them says so.

## Palettes

Fractal → Palette asks for an algorithm number, and says what each one is
called after the number. 1 to 3 are XaoS's own; 4 to 7 were rebuilt after
measuring what the first three are made of, and 8 after measuring paintings.

| | |
| --- | --- |
| 1 | **dark and colour** — a near-black with a tint of its own and a colour, by turns |
| 2 | **black, colour, white** |
| 3 | **warm over dark green** — bright violets, reds and oranges over dark greens and teals, black or white every third |
| 4 | **smog** — dark and melancholy: greys, ochre, slate and dusty violet, with oxblood and scarlet among them |
| 5 | **warm over night** — rose, red, burgundy, brown, orange, ochre and yellow over teal, petrol, blue and indigo, in 3's skeleton |
| 6 | **favoured pairs** — the pairs of colours the first three put side by side most, set down whole between black and white |
| 7 | **random colours** — every stop a colour truly at random, and nothing else |
| 8 | **Kandinsky** — the gradients of his paintings, from mat to vivid and now and then to all but grey |

### What the first three are made of

Twenty thousand palettes of each, made as the program makes them, every colour
named by the nearest of three dozen reference colours in CIELAB, the colours
and the pairs of colours counted:

* **Anchors.** Black and white are a third to a half of every stop, and the
  colours sit between them with a period of two or three stops, so the pattern
  shows in the first few bands — which, in a palette some three thousand stops
  long, are all most pictures ever use.
* **1 and 2 hold 256 colours and no more.** They take the bottom byte of the
  generator, and the bottom byte of that generator is a generator of period 256
  on its own: red decides green and blue, and every colour is always followed by
  the same one. Out of twenty thousand palettes of 2, no colour was once
  followed by a different one. That is where their recurring pairs come from —
  pink beside yellow twenty-four times as often as chance, burgundy beside olive
  six, red beside petrol eight.
* **3 draws from two windows**: bright colours from violet round through red to
  orange, and dark ones from green round to teal. Its commonest pairs are a
  bright warm colour on a dark green.

### The four after them

Each was measured against the first three before it stayed: the Jensen-Shannon
divergence between the colours two algorithms show, by name, and between the
pairs of neighbours they show, over five thousand palettes — 0 for the same
palette, 1 for nothing in common. 1 and 2 stand 0.11 apart by that measure, 1
and 3 0.16.

* **4, smog**, is dark and spent where 1 is vivid on black. Deep stops are a
  colour barely above black, light ones an overcast grey with a tint and never
  white, and every fourth stop a plain grey. Its hues are walked from a ring of
  256 as 1 and 2 walk theirs, so neighbours recur as theirs do, but only over
  the hues of smog — ochre and olive, and steel, slate and a dusty violet — and,
  for a quarter of the ring, crimson to scarlet: oxblood when deep and a dimmed
  scarlet when light, kept saturated, since a dark red that is not reads as
  brown. One stop in five reads as red, against one in thirty in 1. It stands
  0.43 to 0.56 from the first three. It was first 1's mechanism on another ring
  of 256 colours, chosen as far from the classic ring as a ring can be, and
  measured 0.05 from 1 — nearer than 2 is: the likeness was the mechanism,
  black half the time and a colour from anywhere the other half, not the colours.
* **5, warm over night**, is 3's skeleton with the temperature turned round:
  warm colours of every brightness, from rose through red, burgundy, brown and
  orange to ochre and yellow, over night colours kept dark, from teal through
  petrol and blue to indigo. Both windows are 125 degrees wide, near 3's 135;
  they were a third of that at first and read as two colours. Navy is its
  commonest colour, and navy with olive, burgundy or brown its commonest pairs.
* **6, favoured pairs**, takes every pair of colours one of the first three puts
  one after the other at least three times as often as chance — sixty-six of
  them, with the two colours as they actually came out — and sets each pair down
  whole, the two colours touching, between a black stop and a white one.
* **7, random colours**, is every stop a colour drawn at random from the top of
  the generator, where 1 and 2 draw from the bottom: seventy-nine thousand
  different colours in twenty thousand palettes where 2 has 256, and no pair of
  them more than a fifth more common than chance. It had black and white every
  third stop at first, and measured 0.00 from 2: the rhythm of very dark and
  very light was all anyone could see. Without it, it has less range from dark
  to light than the anchored ones, which is what random means.

The four they replace were a spectrum, a duotone, a triad and a complementary
pair, reasoned out from the colour circle, and on the screen each was one
colour. They spread their hues along the whole palette, and the palette is some
three thousand stops long in the program — four to nine in the test that was
meant to catch it, which is why it did not — so a picture saw the first hue
alone. And they kept every stop a colour, deep or pale, with no black and no
white, so nothing in them had an edge. The test now makes its palettes as the
program does, and asks of all seven that the part a picture sees go from dark
to light and hold more than one colour; 8, which came after, is asked a little
less, and why is below.

A position saved with 4 to 7 before this change comes back in different
colours: the number and the seed are all a position records of its palette, and
what the number means has changed. One saved with 1 to 3 is unaffected, to the
bit, and one saved with 4 to 7 names an algorithm the original XaoS does not
have and will refuse.

### 8, Kandinsky

Made from two sets of pictures in his manner: four vivid vignettes — black,
vermilion, yellow, cerulean, violet and pink on cream paper — and one mat
picture, burnt orange, ochre, brown, sage and petrol on the same paper. The
vignettes give eighteen colours, by k-means in CIELAB, and each has a **mat
form**: the colour of the mat picture nearest it in lightness and hue, chroma
left out, since chroma is what the two differ in. Vermilion becomes burnt
orange, yellow ochre, burgundy brown, black a dark brown, cerulean sage, blue
petrol. The mat picture has no grey, so the grey stays grey. The aqua it seems
to have is a sage, 128,146,121: it reads as aqua only beside the orange.

Every palette draws a **mood**. Four in five fall evenly between the mat forms
and the vivid ones; the fifth go past the mat forms towards grey, down to
fifteen hundredths of their chroma, and one palette in twenty or so is all but
grey.

**Anchors come by accident.** One palette in ten may turn a stop here and there
to plain black or plain white: each stop after the first on its own throw, at a
rate the palette draws between one stop in twenty and one in ten, black or
white as the throw falls — no rhythm and no taking turns, unlike the anchors of
1 to 3. The other nine palettes in ten are exactly what they were before the
anchors came in: the throw was already made for every stop, and left unused.

What follows a colour is where **his gradients** take it. Wherever the paint
of the vignettes changes softly — over eight pixels, with no step larger than a
third of the change — the colours at the two ends were counted, and every stop
draws the next from those counts: paper goes on into a wash of yellow and then
yellow, black into vermilion (three times in ten) and vermilion into burgundy,
blue into mauve and mauve into pink, cerulean into grey aqua. The hard edges
were counted too and left out: most of them are paper beside paper of another
shade, and palettes made from them came out paper and grey. A palette starts on
the paper or on the black.

Over twenty thousand palettes, those of the top tenth of moods stand 0.06 from
the vignettes, by the same divergence as above, and those of the tenth above
nought 0.07 from the mat picture; from the other seven palettes they stand 0.39
to 0.67.

Because a gradient stays a while in one part of the colour circle, and 8 has no
black and white every few stops to fall back on, the tests ask less of it than
of the others. Where the palette is laid out long, forty entries a stop and
more, a picture sees ten stops or so, and one palette in a hundred is then a
single field — ochre with sage in it, petrol with vermilion — short of a third
of the way from black to white or of a second colour. And three stops, which is
what a short palette of 256 entries often is, are a single gradient.

## Watching the orbit

`trap` and `stripe` keep one number about a whole orbit rather than about the
point, which is the kind of quantity the engine's colouring modes cannot reach:
a calculation loop there is compiled to hand back a colour, and only the last
`z` and the one before it survive to be coloured by. Done in the formula it
costs the engine nothing.

Both hand back the value they were given until the last pass the iteration
limit allows, and what they gathered on that one. So a whole formula is

    trap(z^2+c, 3)

and nothing else: the fractal iterates as it would, and on the last pass the
value becomes the trap, which the inside colouring modes then draw. A point
that escapes never reaches that pass and keeps its ordinary outside colour, so
what these draw is the inside.

**`trap(a, shape, centre, size)`** — the nearest the orbit ever came to a
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

**`stripe(a, density)`** — the average of `(sin(density * arg a) + 1) / 2`
along the orbit. `density` defaults to 4 and is how many stripes go round a
turn; a whole number, or they do not meet where the turn closes. An average
over the orbit changes smoothly with the point even where the iteration count
jumps, which is what draws the fibres the method is known for.

    trap(z^2+c,0)                    how near the orbit passed the origin
    trap(z^2+c,3)                    a cross, which draws rays
    trap(z^2+c,4,{0,0},{0.5,0})      rings inside the cardioid
    stripe(z^2+c,6)                  six stripes to a turn

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

## Fractional Brownian Motion

A menu of its own on each side — **Fractal → Incoloring mode → Fractional
Brownian Motion** and the same under Outcoloring mode — holding three modes and
the five numbers they are drawn from.

What it does is make a clean colouring look **used**: the bands keep their shape
and their order and stop being perfect, as though the thing had been left out in
the weather. It is not randomness pixel by pixel, which reads as grain and was
tried and thrown away; it is a noise field with a scale, so what it draws are
marks.

**Outside** the three modes are `fbm + smooth`, `fbm + iter` and `fbm` alone —
the smooth count worn, the banded count worn, and the noise on its own, which
draws clouds around the set with no count in them at all.

**Inside** they are `fbm`, `fbm + zmag` and `fbm + decomposition`. The first is
the interesting one: the inside is one flat tone by default, and this gives it a
surface.

### The five numbers

One dialog a side, under **Settings** in the same menu. They are saved with the
position.

| | |
| --- | --- |
| **intensity** | how many bands of colour the noise moves the value by. `4` is the default; below one it is a whisper, above ten the bands stop being bands |
| **frequency** | cells of the lattice to a unit of the plane — the size of the marks. **This is the one to raise as you zoom in**: at the first view `8` reads well, at a span of `0.02` it takes about `200` |
| **octaves** | how many are summed, each at twice the frequency of the one before |
| **roughness** | what each octave keeps of the height of the one before it. `0.5` is the plain motion; higher is grittier |
| **seed** | the same picture every time. Change it for a different one of the same character |

**Octaves and roughness hold each other up.** At a roughness of `0.5` the eighth
octave carries a two hundred and fiftieth of the whole, so octaves past four or
five change nothing you can see — measured, four and eighteen came out the same
picture. Raise the roughness and the fine octaves have something to spend, and
then the count starts to matter.

### What it costs, and the one thing it will not do

The noise is worked out **once a pixel**, in the colouring, not once an
iteration. Over a picture of 90000 points at a span of `0.02` and 400
iterations, `fbm + smooth` took a fifth longer than plain `smooth` at four
octaves and a quarter longer at eighteen; at 4000 iterations, where the orbit is
doing more of the work, the same four octaves cost a tenth.

The marks are anchored **to the plane**, not to the view: it is dirt on the
thing, not on the lens. So zooming in magnifies them, and if you go far enough
you end up inside one mark and the picture is clean again. Raising the frequency
is the answer, and it is why the frequency is a setting rather than a constant.
It cannot be tied to the zoom instead: the engine reuses pixels computed at the
previous scale, and a colouring that moved with the view would leave them
carrying the wrong marks.

The modes are numbered **after** true colour, so every mode number a saved
position already holds still means what it meant.

## More colouring modes

Calculation is not what decides the colour of a pixel: the colour comes from a
formula in the final `z`, the parameter `c` and the iteration count, and there
is one such formula per mode. Eighteen more of them, in a submenu of their own
under Fractal → Incoloring mode and Fractal → Outcoloring mode, above the
true-colour submenus.

**Inside**, where the orbit never escaped and only the final `z` and `c` say
anything:

| mode | |
| --- | --- |
| `real/mag` | the cosine of the angle: shading that closes on itself where `atan2` has a seam |
| `max(\|real\|,\|imag\|)` | the square norm, so the bands are squares where `zmag`'s are circles |
| `\|real\|+\|imag\|` | the same bands turned by an eighth of a turn |
| `min(\|real\|,\|imag\|)` | near zero along either axis, so rays |
| `\|z-c\|` | how far the orbit settled from the parameter, which tells the bulbs apart |
| `\|z*c\|` | the two moduli against each other |
| `angle(z)-angle(c)` | the settled direction measured against the parameter |
| `real*imag` | a saddle, four lobes about the axes |
| `sin(real)*sin(imag)` | a grid laid over where the orbit settled |
| `sign(imag)` | two flat tones: a decomposition of the inside, as binary decomposition does to the outside |
| `frac(mag)` | contours at every eighth of the modulus |
| `log(mag)` | `zmag` with the contrast moved onto the small values, where an attracting orbit spends its time |

**Outside**, where the count and the escape point are what there is:

| mode | |
| --- | --- |
| `iter+angle` | the direction it left by, blended into the count: the bands acquire a twist |
| `iter+log(mag)` | how far past the bailout it went, a cheaper smoothing than the log of a log |
| `iter+real*imag` | biomorphs made continuous rather than thresholded |
| `max(\|real\|,\|imag\|)` | which side it left by, with no count at all |
| `iter banded` | the count folded into eight bands, showing the shape of the level sets rather than their number |
| `\|real\|-\|imag\|` | how lopsided the escape point is |

The numbers the older modes have do not move, so a saved position goes on
loading as it did.

Colouring speed and colouring shift now reach smooth and smooth log, which
they never did: those two returned a pixel of their own rather than going
through the step the other modes end with, so both controls did nothing at all
while either was chosen. At the settings they come with — speed one, shift
nothing — the step is the identity, so no picture that was saved before moves.

Smooth and smooth log follow the bailout shape. They interpolate between one
pass and the next by asking where the orbit crossed the bailout, and they used
to ask it of the modulus whatever shape had been chosen — so escaping on a
square and interpolating on a circle put better than a fifth of the escaping
points, and half of them with a diamond, outside the band they belong to. Each
shape now hands over the quantity it actually tests and the threshold it tests
it against. For the circle the two are the same thing, so nothing that was
drawn before moves.

They also work in Newton mode, where an orbit converges rather than escapes: the
quantity is then how far the last step moved and the threshold the convergence
limit, and the interpolation is the same question asked of a number going down
instead of up. The two Newton fractals, which had no smooth variant at all, have
one.

**Newton convergence**, in the Calculation menu, is read at last. The limit was
written into the iteration as a millionth, so the setting was offered, saved and
reloaded while nothing ever consulted it.

**Smooth and smooth log** have moved into that submenu on the outside menu.
They interpolate between one iteration and the next using how far past the
bailout the orbit went, which needs a formula that escapes on the bailout —
twenty of the thirty-one do, and the other eleven either converge (the Newtons)
or stop on a test of their own (the Sierpinskis, Koch, the Clock). For those
the mode does nothing at all, which is not something the first list should
offer. The user formula, which could always have had it and did not, has it
now.

What these modes cannot do is anything that needs the orbit rather than its
end: an orbit trap, a stripe average, a Lyapunov exponent, distance
estimation. Those need a value carried from one iteration to the next, which
in this engine means a second compiled copy of every calculation loop. `trap`
and `stripe` in a user formula do it instead.

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
