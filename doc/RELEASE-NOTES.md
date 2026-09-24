# MS XaoS 1.7.4

What has changed since 1.7.2. Version 1.7.3 was never released on its own, so
what it brought is here too. Everything MS XaoS adds over
[XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the notes of
[1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7) and
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v1.7.4/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## No more crash on the way out

**Every exit crashed, from 1.6 on** — quitting the program and every
`-render` alike, always after the work was done, so the pictures and the
positions were never at risk and nothing showed it. The menus of numbered
entries were removed with a count written a second time, apart from the one
they were made with; when the fbm colouring modes changed the first and not
the second, the removal walked three entries past the end of the inside
colouring block.

* A block of numbered entries now ends with an empty one and knows its own
  length, so the count is written once.
* `-render` says in its exit code whether the render went through.
* The test that renders, saves and reloads positions now fails on a render
  that does not exit cleanly. It looked only at the pictures, which is how
  this went unseen.

## Noise without a grid

The fbm colouring modes, `fbm()` and `randsc` were **value noise**: a number
at each corner of a lattice, blended by a smoothstep. A smoothstep has no
slope at either end, so the field went flat along every line of the lattice,
and every octave's lattice held the first one's lines. Measured over a picture
where it was reported, the field was a third as steep on the lines as between
them: squares of eight pixels, as though the picture had been worked out at a
lower resolution than it was shown at.

All three now draw **gradient noise**, from one place: a direction at each
corner, picked by the hash from sixteen, a quintic blend, and each octave
doubled and shifted by 2 − φ and φ − 1, so that no two octaves share a line or
a corner. The slope across the lines is now 1.00 of the slope between them at
the default settings and 1.11 at those of the reported picture, where across a
line it had been nought.

What was kept, so that saved numbers mean what they meant:

* the values run from nought to the intensity — or to one — and never outside;
* **the size of the marks and blobs**: the lattice lies at three fifths of the
  frequency or of the size, which brings gradient noise to the scale value
  noise drew at — 0.56 of a cell at the default fbm settings, eleven sixteenths
  of the size for `randsc`, before and after;
* **the contrast**: gradient noise alone swings seven tenths as far from its
  middle, and a monotone curve takes every octave back onto the distribution
  value noise had — `randsc` 0.2147 from the middle against 0.2145 before, the
  default fbm 0.133 against 0.132. Monotone and smooth, it moves no outline and
  draws no crease. `tools/fbm-contrast-table.py` measures it and writes it;
* the cost, near enough, on the machine it was written on: the colouring modes
  159 ns a pixel at four octaves and 452 at eighteen against 160 and 479,
  `fbm()` 182 and 498 against 161 and 438, `randsc` 167 against 147.

What was not: a position that uses an fbm mode, `fbm()` or `randsc` keeps its
numbers, its scale and its strength and **draws different marks**.

## A picture of every tiling

The user formula reference has a **Tilings** tab: a picture of each of the
forty-five tilings `randsctile` takes, captioned as the Values tab lists them,
and a double click copies the start of the call — `randsctile(7, ` — to the
clipboard, ready for the seed.

* The pictures are drawn by `randsctile` itself, eight units square round the
  origin, by `tools/randsctile-thumbnails.cpp`, and kept in the binary: 256
  pixels a side, twice what is shown, so they stay sharp on a screen scaled up.
* Kept pictures can go stale, so each is stored with a fingerprint of what the
  function drew, and a test draws the same points again and fails, naming the
  tiling, until they are redrawn. It checks too that the Values tab lists
  exactly the tilings the parser draws.

## Tests

The same 343 test programs, with new checks: the slope across the lattice
lines, the range, the size of the marks and blobs and the spread of the values
of the noise; the tiling pictures against what `randsctile` draws; and a clean
exit from every render.

## Compatibility

* Positions saved by 1.7.2 load as they did, and draw as they did unless they
  use an fbm mode, `fbm()` or `randsc`, which draw the marks described above.
* The five mosaics — `randscq`, `randscp`, `randsch`, `randsct` and
  `randsctile` — are unchanged to the bit, and so are the palettes.
