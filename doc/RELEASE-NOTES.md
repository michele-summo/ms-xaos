# MS XaoS 2.7

What has changed since
[2.5](https://github.com/michele-summo/ms-xaos/releases/tag/v2.5). Version 2.6
was never released on its own, so what it brought is here too. Everything MS
XaoS adds over [XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the
notes of [1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7),
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2),
[1.7.4](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.4),
[1.8.0](https://github.com/michele-summo/ms-xaos/releases/tag/v1.8.0),
[2.3](https://github.com/michele-summo/ms-xaos/releases/tag/v2.3),
[2.4](https://github.com/michele-summo/ms-xaos/releases/tag/v2.4) and
[2.5](https://github.com/michele-summo/ms-xaos/releases/tag/v2.5), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v2.7/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## The formula history in the bar

Fractal → User formula reached the formulas used before only through the
History... window. The bar now steps through them itself: **Up**, from its
first line, brings in the formula used before the one shown, and **Down**,
from its last line, goes back towards the newest and then to what the bar held
before stepping began, as in a shell. Two small arrows beside the bar do the
same, each greyed when there is nothing further that way. A formula already in
the bar is passed over, so the first Up always changes something, and typing
makes the typed text the place the steps start from again. The steps go
through the bar's own editing, so they can be undone.

**The formula used last is shown faint** in the bar: all of it while the bar
is empty, and the rest of it after what you have typed for as long as what you
type is how it begins. Tab puts it in the bar whole; once what you type goes
another way the faint text goes too, and Tab moves on to the buttons as
before. The faint text is laid out on the same lines as the whole formula, so
each typed character sits over the one it completes however many lines it
takes.

The bar now opens empty, with that formula faint in it, and OK on an empty bar
uses it. Only the very first time, before any formula has been used, does it
open with the formula in use, `z^2+c`.

## Tessellation

**Fractal → Incoloring mode → Tessellation...** colours the inside by the
tiles of any of the forty-five tilings of `randsctile`, laid over the final
`z` forty to a unit as the squares of `squares` are. Its dialog has three
fields:

* **Tiling**, by the number `randsctile` knows it by and its name, `3:
  hexagons` where it starts. The Tilings tab of the formula reference has a
  picture of each.
* **Coloring**, what the colour comes from: `angle`, the angle of `z`, where
  it starts; `tile-angle`, the angle about the middle of the tile, which makes
  each tile a small wheel of colour; `zmag`; or `none`, the tiles alone.
* **Tiles**, how each tile stands apart: `turned`, moved along the palette by
  an amount of its own, where it starts; `shaded`, with a bevel from its middle
  to its edge; `edged`, its outer part drawn in the other half of the palette;
  `flat`, one colour of its own.

The frequency and function of the periodic modes size the tiles, and the
decomposition settings below cut and walk the two angles. The periodic tilings
cost about what the other inside modes do, 120 to 210 ns a pixel; the ones
that never repeat, located afresh at every point, 250 to 490. A position saves
the mode as `(incoloring 30)`, the tiling as `(intessellation 39)` only when
it is not hexagons, and the colouring as `(intessellstyle 'tile-angle 'edged)`
only when it has moved.

## Decomposition settings

**Decomposition...** in Incoloring mode and in Outcoloring mode sets how the
decompositions cut the turn: binary decomposition and color decomposition
outside, `sign(imag)`, `Decomposition-like` and `fbm + decomposition` inside,
and the angles of tessellation. Its fields are the sectors to a repetition (an
odd number taken up to the even one above in the binary ones, where two
sectors of a colour would meet), the rotation, inside the angle of `z`, `z-c`
or `z*c`, rings per unit inside or iterations per band outside, the
repetitions, the direction (right, left, or both, a mirror) and a speed: more
colours to a repetition in the colour ones, the distance between alternate
sectors' colours in the binary ones. With binary decomposition outside, the
bands and sectors together draw the grid of external rays and equipotentials.

With every field where it starts the old code runs, and the five draw as 2.5
does, pixel for pixel. They are saved as `(indecomposition ...)` and
`(outdecomposition ...)` only when moved. Moved, a decomposition is no longer
symmetric about the real axis, and XaoS stops drawing one half as the mirror
of the other.

## More periodic modes, and a function for them

Four new inside modes repeat, in the Other coloring mode submenu: `cos(|z|)`,
rings about where the orbit settled; `hexagons`, three waves a third of a turn
apart; `spiral`, five arms wound once to every factor of e in the modulus; and
`cos(|z-c|)`, rings about the parameter. With tessellation, eleven inside
modes now read **Frequency (periodic modes)**.

**Function (periodic modes)**, under it in the inside Coloring speed dialog,
bends the argument of every periodic function of the eleven by one of the
twelve colouring speed functions before the frequency multiplies it: `square`
and `cube` crowd the bands towards the edge, `log` and `atan` the other way,
and `sin`, `cos` and `triangle` fold them back on themselves. At `ident`, where
it starts, nothing changes; it is saved as `(incolorperiodicfun 6)` only when
it is not ident. Both Coloring speed dialogs are now laid out in two columns,
the speed change after ITER on the right.

## Fractional Brownian Motion sized by the count

The outside's Fractional Brownian Motion dialog has two more fields. **Noise
scaling** says what sizes the marks: `fixed`, by the frequency alone as
before, where it starts; `iteration`, by the frequency raised with the escape
count; `smooth`, with the smooth count. **Scaling factor** is how many times
the frequency doubles every ten iterations, so the marks grow finer towards
the edge of the set and follow its detail. The field is read at the two powers
of two either side of the count and the two are blended, as a texture is
between its mipmaps, so `smooth` changes the size without a seam and
`iteration` gives every band a grain of its own.

Between two doublings the colouring costs about twice the fixed one, 248 ns a
pixel against 127 at four octaves. `fixed`, or a factor of nought, draws as
before to the bit, and the setting is saved as `(outfbmscale 'smooth 1.5)`
only when it has moved.

## Reopening a dialog from a ticked entry

Clicking Tessellation... or Fractional Brownian Motion... while that mode was
already in use did nothing, so there was no way to change the settings of the
mode in use. The menu took a ticked radio entry for a checkbox to switch off, a
mistake that came from upstream. A radio entry now always opens its dialog; a
ticked checkbox, such as Julia mode or Perturbation, still turns off without
one.

## Linux and macOS

MS XaoS builds with CMake and Qt 6 on Linux and on macOS as well as on
Windows, and the qmake project links again. On macOS the build makes an
application bundle with the catalogs, the tutorial and the examples inside,
and a `deploy` target that runs macdeployqt and makes a `.dmg`; the README
says how. There is no quad build on the Mac, Apple's compiler having no
128-bit float, and on Apple Silicon `long double` is a `double`, 53 bits: the
ordinary build there goes as deep as a 64-bit double build does, not as deep
as XaoS.exe. The tests know about the 53 bits.

GitHub builds every commit on Linux and on macOS (Apple Silicon, Qt 6.10),
with CMake and with qmake, and runs the tests on both.

## Translations and the guide

Every new text of the interface is translated into German, Spanish, French,
Italian, Portuguese and Russian, and the guide describes everything above in
English and in the six translations.

## Tests

359 tests, six more than 2.5: the decomposition settings, the scaling of
Fractional Brownian Motion and tessellation, each at both precisions. The five
decompositions, the fixed fbm modes and tessellation where it starts are
pinned to the checksums they gave before the change; each setting is checked
to move what it belongs to and nothing else; the doublings of the fbm
frequency are exact, and every tiling has a tile everywhere. The frequency
test covers the eleven periodic modes and the function. The round trip of
saved positions has six more.

The dialogs were tried in the running program in Italian and German, from
their menus to the picture after OK.

## Compatibility

* Positions saved by 2.5 load and draw as 2.5 draws them, pixel for pixel.
* Those saved by 2.7 open in 2.5 unless they use something new: the five new
  inside modes, or a periodic function, decomposition setting, fbm scaling or
  tessellation setting that has moved from where it starts.
* The download for Windows holds XaoS.exe and XaoS-quad.exe as before; the
  one for macOS holds the ordinary build alone.
