# MS XaoS

A fork of [XaoS](https://github.com/xaos-project/XaoS) 4.3.3. Everything the
original does, it still does; this describes what has been added or changed,
and why. Version 2.5.

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
function that takes them.

Each row of the function list says what the call takes and in what order —
"a; b", or "seed; size; degradation" where the position means something. The
list and those counts are both checked against the parser’s own table by a
test, so a function that gains or loses an argument cannot go on being
described with the one it used to have. An argument that may be left out is
shown in brackets with the value it takes when it is: `a; [b=1]; [c=1]`.

**Leaving a place empty.** An argument shown in brackets may be left out in the
middle of a call as well as at the end, by writing nothing between the two
separators: `julian(z; ;3)` gives the first and the third and lets the function
say what the second is. Spaces make no difference — `f(z, ,5)` and `f(z,,5)`
are the same call — and an argument that is not in brackets must still be
written, so `poly( ;1;2)` and `sin(,)` are refused rather than guessed at.

What an empty place means is the function's own business. Most take the default
they declare; a coefficient of `poly` left empty is a term that is not there,
so `poly(z;1; ;1)` is `z^2 + 1`; and a branch of `ifiter` or `ifiterl` left
empty repeats the one before it, which is how a branch is given more than one
pass in the cycle. `ifiter(f(z); ; ;g(z); ; ; ; )` runs `f` for three passes
and `g` for five, and costs no more to evaluate than writing them out would:
where each choice leads is settled once, when the formula is parsed.

A user formula that has not been written yet says `z^2+c`. It used to say the
burning ship, which is a fractal of its own and a puzzle to meet as a starting
point.

**`ifiterf(a; b)`** — evaluates `a` on every pass but the final one and `b`
on that. The final pass is the last the iteration limit allows: a formula has
no way of knowing which pass will be the one that escapes, that depending on
the value it has not produced yet. Only the chosen one is evaluated.

**`ifiterr(a; b; n)`** — evaluates `a` while the pass number is below `n`
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

**Removed.** `powi`, `powdc` and `logcn` were second names for `pow` and
`logn`; `rad`, `deg` and `sign` were listed with no implementation behind them,
and a second `trunc` was shadowed by the working one. A formula using a removed
alias should use the name that remains — the function is identical.

**`erf(z)`** — the error function over the complex plane. Accurate to about
three ulp inside a bailout of two, where a fractal actually iterates.

**`poly(z; k1; k2; ...; km)`** — a polynomial in `z`:

    k1*z^(m-1) + k2*z^(m-2) + ... + k(m-1)*z + km

The first coefficient written multiplies the highest power and the last stands
alone, so the call reads in the order one says the polynomial. Worked out by
Horner's rule, which is m−1 multiplications rather than the m(m−1)/2 that
raising each power separately would take, and the more accurate of the two into
the bargain.

    poly(z;1;0;0)+c            the Mandelbrot, written out
    poly(z;1;0;0;{0.7,0.2})    z^3 + 0.7+0.2i

**`randsc(seed; size; degradation; kaleidoscope; mode)`** — coherent
noise over the point, giving blobs rather than per-pixel snow. `size`
(default `1+i`) is the average width of a blob along the real axis and its
height along the imaginary one.
`degradation` (default `0.5+0.5i`, halving them each pass) shrinks them as
the iteration proceeds: the
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
first two levels are five twelfths each and the rest come down evenly. So the
step onto the parent is a **fold** for a triangle standing on the hexagon —
across the side it stands on, which lands it exactly on one of the six the
hexagon is made of — and, for one standing on the free edge of another triangle,
a blowing-up by three **and a turn**, each free edge facing a different way.

The point really travels, and that is what makes these usable rather than merely
correct: `z` along the way is a point of the plane like any other, so the
colouring modes that read it, smooth colouring, and writing a figure inside a
larger formula all mean something.

A gasket and a carpet fill their shape, so every point in one has a level and
leaves on it. A snowflake does not fill its hexagon: it leaves six corners of
**ground**, and that ground is no part of the figure and has no level. It is
handed back exactly where it stands, so it never leaves and comes out in the
**inside colour** — the incolouring modes work there, on the point itself, and
the empty space stays empty instead of taking a band of the outside colour off
the figure.

That is the one place these figures cost something. A ground pixel is asked the
membership question once a pass for as long as it is looked at, where the other
two figures and the body of the snowflake leave on the pass their level names.
Drawn straight through the iteration loop, a whole snowflake costs about 2.4
times what it cost when the ground was given a band of its own. Boundary tracing
is left on for these formulas — their bands really are solid, which is the
condition for it — and the ground is one solid region, so what reaches the
screen is less than that number suggests.

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

**`sierpinskyt([radius=4])`** — the Sierpinski gasket, in that triangle.

**`sierpinskyc([radius=4]; [squares=3])`** — the Sierpinski carpet, in that
square. The square is cut into `squares` by `squares`, the ring of cells along
the border is kept, everything that ring encloses is thrown away, and the same
is done to each cell that was kept. So the picture is **one square in the middle
and 4×`squares`−4 around it**: eight at three, which is the carpet as it is
usually drawn, twelve at four, sixteen at five. Two is the one number with no
ring to speak of, and there the far corner goes instead, which is a gasket again
— a square cut in four with one corner taken away is what a gasket is.

**`snowflake([radius=4])`** — the Koch snowflake, its points on that hexagon's
corners, banded by generation from its middle out: the hexagon on the first
pass, the six triangles on its sides on the second, the twelve on their free
edges on the third, the forty-eight on the fourth. The six corners of ground the
figure does not cover are not banded at all — they never leave, and are drawn in
the inside colour.

Every argument has a default, so `snowflake()` is a call, and so is
`sierpinskyc( ;5)` — a lacier carpet at the default size.

    sierpinskyt()        with Fractal -> Bailout shape -> triangle -90
    sierpinskyc()        with bailout shape square
    sierpinskyc( ;5)     a lacier carpet
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

Fractal → Palette asks for an algorithm number. There were three, and all three
work the same way: some four to nine colours are chosen, every second or third
of them pinned to black or white, and the palette is interpolated between them.
That is what gives XaoS its banded, high-contrast look, and it was the only look
it had.

Four more choose their colours in relation to each other:

| | |
| --- | --- |
| 1–3 | colours scattered between black and white anchors — as before |
| 4 | **spectrum** — right round the hue circle, one turn, darkening at both ends |
| 5 | **duotone** — two hues, one owning the shadows and the other the highlights, as a press does it with two inks, banded |
| 6 | **triad** — three hues spread round the circle, taken in turn, deep and bright alternating |
| 7 | **complementary** — two hues from opposite sides of the circle, alternating |

They cost what the others cost: a palette is made once, when it is asked for,
and none of these does more per colour than one conversion out of hue,
saturation and value. Each is driven from the seed the dialog holds, so an
algorithm and a seed give the same palette every time — which is all a saved
position records of its colours.

Three things had to be settled by hand, and a test now keeps them settled.

A hue cycle at one brightness has no dark anywhere in it, and a fractal shown
in it has hue where it should have shape — so the spectrum swells from deep to
bright and back, with saturation running the other way. Holding saturation
steady instead, the channel sum could not span more than twice the value, and
the palette was a fifth of what it can be.

Every segment keeps its colour. The three older ways get their character by
pinning segments to black and white, and a segment pinned to either has no hue
at all; done in a palette whose whole point is which hues it holds, it leaves
the hues invisible. Two of these four went that way — alternating near-black
with a washed-out near-white — and came out grey with a tint, which is what
"monochrome" means when it is a complaint. Contrast here is between a deep
colour and a bright one, never between nothing and nothing.

And which segment is deep and which is bright is settled by where the segment
sits, not by the dice: left to the dice, a palette four segments long comes out
all one weight more often than not.

Deep and bright alternate **band by band**, as they do in the three older ways,
rather than swelling once from one end of the palette to the other. A palette
that swells once changes colour over a hundred entries where a banded one
changes over thirty, and a fractal drawn in the first has no edges to its
rings — the gradient reads as slow. The test measures how far a palette travels
in brightness against how far it reaches, and asks that the new four travel as
far as the old three do, seed for seed.

A position saved with one of the new four names an algorithm the original XaoS
does not have and will refuse to load; one saved with 1 to 3 is unaffected.

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
