# MS XaoS 2.11

What has changed since
[2.7](https://github.com/michele-summo/ms-xaos/releases/tag/v2.7). Versions
2.8, 2.9 and 2.10 were never released on their own, so what they brought is
here too. Everything MS XaoS adds over
[XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the notes of
[1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7),
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2),
[1.7.4](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.4),
[1.8.0](https://github.com/michele-summo/ms-xaos/releases/tag/v1.8.0),
[2.3](https://github.com/michele-summo/ms-xaos/releases/tag/v2.3),
[2.4](https://github.com/michele-summo/ms-xaos/releases/tag/v2.4),
[2.5](https://github.com/michele-summo/ms-xaos/releases/tag/v2.5) and
[2.7](https://github.com/michele-summo/ms-xaos/releases/tag/v2.7), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v2.11/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## Tessellation outside

**Fractal → Outcoloring mode → Tessellation...** has three modes, as
Fractional Brownian Motion has: `tessellation + smooth`, the tiles worn by the
smooth count; `tessellation + iter`, by the whole count; and `tessellation`,
the tiles alone. The tiles are laid over the pixel's place in the plane, not
over the escape point, which jumps from one band to the next and would give
only grain. The dialog has the inside's fields, with settings of their own, and
three more:

* **Frequency**, tiles to a unit of the plane, `8` where it starts. It has to
  be raised as you zoom in, as the motion's frequency does.
* **Tile scaling**: `fixed`, by the frequency alone; `iteration` or `smooth`,
  by the frequency doubled with that count, so that the tiles grow finer
  towards the edge of the set.
* **Scaling factor**, how many times the frequency doubles every ten
  iterations. A factor of a few tenths is a good start; at `1` the tiles fall
  below a pixel soon.

Between two doublings the tiles are read at both sizes and their two places on
the palette are blended, which doubles the cost of the tiling. `angle` is the
angle of the escape point and `zmag` the log of its modulus. Both are read a
few iterations past the bailout and run on across the edges of the iterations
instead of breaking there, except in `tessellation + iter` with `fixed` or
`iteration` scaling. A position saves the mode as `(outcoloring 21)` to
`(outcoloring 23)` and the settings in one command,
`(outtessellation 5 'tile-angle 'shaded 0.6 40 'smooth 0.4)`, only when one of
them has moved.

## Intensity

Both Tessellation dialogs, inside and outside, have an **Intensity** field: how
far apart along the palette the tiles are set. It multiplies every amount a
tile has of its own, its turn, the bevel of `shaded` and the edge of `edged`,
and nothing else. It takes any real value: at `0` only what the colour comes
from is left, so `angle` is a plain decomposition and `zmag` plain `zmag`; at
`1`, where it starts, the tiles are as before; above `1` they stand further
apart.

The intensity moves colours along the palette and never mixes them. Measured
with **View → Palette probe (click a point)**, at `1` a tile is moved by up to
about ten places, in every style and on any palette, so that `0.1` moves it by
up to one place and `0.01` by up to a tenth. On the default palette the
colours change within a place or two, so the range that reads as a stronger or
fainter tiling is about `0.01` to `0.5`; above that the tiles' colours are
already as far apart as the palette makes them. The inside saves it as
`(intessellintensity 0.5)`, only when it is not one.

## One scale for every tile style

A tile's own amount now spans the same part of the palette whatever its colour
comes from and whatever the size of the palette: as much as a turn of the
angle, about ten places. Before, a `flat` tile, `zmag` and `none` took the
whole palette, which on the default palette of 65534 entries made a hundredth
of intensity eighty places, and an `edged` tile was moved by half the palette.
They now take what `angle` with `turned` takes, and `edged` half of that, so
the intensity moves every style alike and `none` with `turned` is the same as
`flat`.

## Every pixel a colour of the palette

Where the outside blends two things, the two sizes of tiles between two
doublings or a value and the same value one iteration on, it blends their
places on the palette, the short way round, and looks the colour up there. It
never mixes two colours in RGB, so every pixel is a colour of the palette and
the palette probe names it. Where two places lie far apart, the blend runs
through the palette's colours between them.

## Penrose examples

`examples/MS-XaoS` holds twenty positions on the Penrose rhombs, tiling 39,
ten inside and ten outside: the same flat tiles at intensity `0.01`, `0.1`,
`1` and `141.42`, then every colouring and tile style and, outside, every
scaling and the three modes. They are installed with the other examples.

## Translations and the guide

Every new text of the interface is translated into German, Spanish, French,
Italian, Portuguese and Russian, and the guide describes everything above in
English and in the six translations.

## Tests

359 tests, as in 2.7, with more in the tessellation suites. A pixel of the
outside drawn in truecolor must be the colour of the palette entry the probe
names for it, for every colouring, every tile style, intensities from `0.01` to
`141.42` and every scaling; between two sizes the probe's place must lie on the
short way between the two sizes' places; and a tile's own amount must be the
same on a palette of 200 entries and on one of 65535. The round trip of saved
positions now also goes through the twenty examples, forty positions in all.

## Compatibility

* Positions saved by 2.7 load and draw as 2.7 draws them, except inside
  tessellation with `flat`, `zmag`, `none` or `edged` tiles: their tiles now
  span ten places of the palette instead of the whole of it, as described
  above. With `angle` or `tile-angle` and `turned` or `shaded` tiles they draw
  as before, pixel for pixel.
* Those saved by 2.11 open in 2.7 unless they use outside tessellation or an
  intensity other than one.
* The download for Windows holds XaoS.exe and XaoS-quad.exe as before; the one
  for macOS holds the ordinary build alone.
