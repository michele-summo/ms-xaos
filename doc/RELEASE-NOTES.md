# MS XaoS 2.14

What has changed since
[2.11](https://github.com/michele-summo/ms-xaos/releases/tag/v2.11). Versions
2.12 and 2.13 were never released on their own, so what they brought is here
too. Everything MS XaoS adds over
[XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the notes of
[1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7),
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2),
[1.7.4](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.4),
[1.8.0](https://github.com/michele-summo/ms-xaos/releases/tag/v1.8.0),
[2.3](https://github.com/michele-summo/ms-xaos/releases/tag/v2.3),
[2.4](https://github.com/michele-summo/ms-xaos/releases/tag/v2.4),
[2.5](https://github.com/michele-summo/ms-xaos/releases/tag/v2.5),
[2.7](https://github.com/michele-summo/ms-xaos/releases/tag/v2.7) and
[2.11](https://github.com/michele-summo/ms-xaos/releases/tag/v2.11), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v2.14/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## The formula reference at hand

**Fractal → User formula** has a **Formula reference...** button beside OK,
which opens the reference at the tab it was left at. The inside and outside
**Tessellation...** dialogs have one too, which opens it at the Tilings tab.
The dialogs wait for you while the reference, a window of its own, is
scrolled and read.

## Turning a value: `rotate`

`rotate(z, angle)` turns `z` clockwise about the origin by `angle` radians,
`z*e^(-i*angle)`, whatever the view's own angle; a picture drawn from a point
put through it turns clockwise with it. A whole number of quarter turns is
exact, so a whole turn is no turn at all, to the bit, even when the angle
comes out of a division such as `180/180*pi`.

Two variable suffixes say it shorter. `_rNtD` turns by *N*/*D* of half a turn,
so `p4_r3t7` is `rotate(p4,3/7*pi)`; `_rN`, without the `t`, turns by *N*
degrees, so `z_r45` is `rotate(z,45/180*pi)`. They chain with the other
suffixes, read from left to right.

## Choosing among values: argmin and argmax

Six functions take any number of arguments and choose among them:

* `argminm(z1, z2, ...)` and `argmaxm(z1, z2, ...)` hand back whole the value
  whose modulus is the smallest or the largest.
* `argminx(z1, z2, ...)` and `argmaxx(z1, z2, ...)` take the smallest or the
  largest real part and imaginary part, each on its own, so the two parts may
  come from different values: `argminx({1,2}, {-3,5}, {4,-6})` is `{-3,-6}`.
* `argminxa(z1, z2, ...)` and `argmaxxa(z1, z2, ...)` do the same by absolute
  value and keep the sign: `argminxa({1,2}, {-3,5}, {4,-6})` is `{1,2}`.

Of two alike, every one of them keeps the one written first, `x` and `-x`
included for the absolute values, so `argminxa({1,2}, {-1,-1}, {2,1})` is
`{1,-1}`. Two moduli or absolute values a few units in the last place apart
count as alike: values equal in exact arithmetic seldom come out equal to the
bit, and a choice made by the last bits would draw as noise.
`argminm(z^2, z_r90^2)+c` draws `z^2+c`. An empty place is passed over, and
one value alone is that value.

## A delta for the random poly

The random poly of the reference's Utilities tab has a fifth field,
**Argument delta (degrees)**. Empty, as it starts, the draws are what they
were, and the same seed gives the same poly. Given a delta, the first
coefficient's argument still goes round the whole turn and every other one is
drawn within that many degrees of it.

## The formula bar

**Fractal → User formula** opens, as before, on an empty bar with the formula
used last faint in it, which OK on the empty bar takes. Once a formula has
been chosen there with OK, it opens from then on, until XaoS is closed, with
the formula in use written out. The bar, its two history arrows and its label
now line up with one line of formula, and the arrows and the label stay
beside the first line while the bar grows below.

## Translations and the guide

Every new text of the interface is translated into German, Spanish, French,
Italian, Portuguese and Russian, and the guide describes everything above in
English and in the six translations.

## Tests

435 tests, 76 more than in 2.11: the new functions, the suffixes and the
exact quarter turns in the parser's suite, the suffixes drawn against the
calls written out, the ties of the argmin and argmax functions, and the random
poly's draws with and without a delta.

## Compatibility

* Positions saved by 2.11 load and draw as 2.11 draws them; none of the
  changes moves a setting.
* Those saved by 2.14 open in 2.11 unless their user formula uses `rotate`,
  the `_r` suffixes or the argmin and argmax functions.
* The download for Windows holds XaoS.exe and XaoS-quad.exe as before; the one
  for macOS holds the ordinary build alone.
