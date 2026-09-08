# MS XaoS 1.6

A fork of [XaoS](https://github.com/xaos-project/XaoS) 4.3.3, by Michele Summo.
This is its first release, so what follows is everything it adds over 4.3.3
rather than the changes of one version.

The guide is [doc/ms-xaos-guide.md](ms-xaos-guide.md); it says what each of
these does and, where a choice was not obvious, why it was made that way.

## Two binaries

`XaoS.exe` is built at `long double`, `XaoS-quad.exe` at 128 bits. The second
zooms about four decades further before the picture goes to blocks and is some
twenty times slower, so it is there for when the first has run out rather than
for every day.

## Writing formulas

**Help → User formula reference** lists every function a user formula may call,
what it takes, and what it does, with a tab of its own for the numbers that
appear as arguments and mean something particular. It is checked against the
parser's own table by a test, so a function that gains or loses an argument
cannot go on being described with the one it used to have.

**Under the formula bar there is a line** saying what the call the cursor is in
takes, with the argument being written in bold. It needs no more of the formula
than an open bracket with a name in front of it, so it is there while the
formula is half written, which is when it is wanted.

**Arguments are separated by a comma, and by nothing else.** A semicolon used
to serve as well; it is now refused where it stands. Six of the shipped
positions were written the old way and have been rewritten — a position of your
own that has a semicolon in its formula will not load until the formula is
written with commas.

**An argument may be left out**, in the middle of a call as well as at the end,
by writing nothing between the two separators: `julian(z, ,3)`. What an empty
place means is the function's own business — most take the default they
declare, a coefficient of `poly` left empty is a term that is not there, and an
empty branch of `ifiter` repeats the one before it.

**`p1` to `p9999`** reach back into the orbit as far as you like, and an
initialization can say where `z` starts.

## New formula functions

| | |
| --- | --- |
| `randsc`, `randscq`, `randscp`, `randsch`, `randsct` | coherent noise over the point: soft blobs, square cells, irregular polygons, hexagons, triangles. Each takes a size, a degradation that shrinks the cells as the iteration proceeds, a kaleidoscope, and a skew that gives the colouring modes something to read inside a cell |
| `fbm` | a fractional Brownian motion over a point the formula names: octaves of the same noise, each at twice the frequency of the one before |
| `sierpinskyt`, `sierpinskyc`, `snowflake` | the gasket, the carpet and the Koch snowflake as fractals of their own, banded by generation, each inscribed in the bailout shape of the same number |
| `trap`, `stripe` | how near the orbit ever came to a shape, and the average of a sine along it |
| `poly` | a polynomial in `z`, written out by its coefficients |
| `ifiterf`, `ifiterr` | one branch on the last iteration and another on the rest; one below a threshold and another above it. Only the chosen branch runs |

## Colouring

**Eighteen more colouring modes**, inside and out, in a submenu of their own
above the true-colour submenus.

**Fractional Brownian Motion**, at the foot of Other coloring mode on each
side: three modes a side and the five numbers they are drawn from — intensity,
frequency, octaves, roughness and seed, kept apart for the inside and the
outside and saved with the position. What it does is make a clean colouring
look used: the bands keep their shape and their order and stop being perfect.

**The colouring speed and the shift now reach smooth**, which they never did:
that mode returns a pixel of its own and skipped them.

## Bailout shapes

Fractal → Bailout shape offers the circle and five polygons — square, triangle,
hexagon and the rest — measured by the apothem, so a polygon of a given bailout
stands its sides where the circle of that bailout would be. User formulas
honour the shape too, and smooth colouring follows it.

## Palettes

Four more palettes, and a **palette probe** under View: click a point and the
message line says where in the palette that point was drawn from, as a number
rather than as a colour. What it is most useful for is seeing **which part of
the palette a picture actually uses** — with the plain count colouring at sixty
iterations the answer is the first eight of the editor's thirty-one colours,
and editing the others changes nothing.

## Getting about

**Selection zoom**: drag a rectangle. And four fixed steps, in and out by two
and by ten.

The **A** key no longer starts the autopilot.

## Elsewhere

Antialiasing can average its samples in **linear light**. Numbers in the
dialogs are shown and read at the precision the build has, and the notation is
chosen by the size of the number rather than by the count of its digits.

## Tests

343 of them, run at both precisions: the iteration loops against golden
checksums, the noise fields and the figures, the parser against a corpus that
includes every formula XaoS ships, the reference against the parser's own
table, and the line under the formula bar.

## Compatibility

* Positions saved by XaoS 4.3.3 load unchanged, **except** one whose user
  formula separates arguments with a semicolon.
* The colouring mode is saved as a number. The modes added here are written
  after true colour rather than before it, so every number an existing position
  holds still means what it meant.
