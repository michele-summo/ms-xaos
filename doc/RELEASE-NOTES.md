# MS XaoS 1.7.2

What has changed since 1.7. Version 1.7.1 was never released on its own, so
what it brought is here too. Everything MS XaoS adds over
[XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the notes of
[1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6) and
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v1.7.2/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## Tilings: `randsctile`

**`randsctile(tiling, seed, ...)`** is the noise field of the `randsc` family
over any of **forty-five tilings**, which the first argument chooses.
Everything after it is what `randsc` takes and means the same, one place
further along, so `selfsim` is the ninth. Every tiling is scaled to a tile of
unit area on average, so the first argument changes the shape of the cells and
not their scale.

| | |
| --- | --- |
| 1–3 | the regular tilings: squares, triangles, hexagons |
| 4–11 | the eight Archimedean ones: 4.8.8, 3.6.3.6, 3.4.6.4, 3.12.12, 4.6.12, 3.3.3.4.4, 3.3.4.3.4, 3.3.3.3.6 |
| 12–19 | their duals: tetrakis square, rhombille, deltoidal trihexagonal, triakis triangular, kisrhombille, prismatic, Cairo and floret pentagons |
| 20–32 | bricks, Flemish bond, herringbone, basketweave, Pythagorean, chevrons, squares and rhombi, houses, rows of squares and triangles in two rhythms, hexagons among triangles, Greek crosses, T tetrominoes |
| 33–37 | Islamic stars: eight-pointed with crosses, six-pointed with hexagons, eight-pointed from 4.8.8, twelve-pointed from 3.12.12 and from 4.6.12 |
| 38 | Voronoi cells |
| 39–43 | tilings that never repeat: Penrose's rhombs, Penrose's kites and darts, Ammann–Beenker, rhombi in twelve and in seven directions |
| 44–45 | the pinwheel and the chair |

* The periodic ones are tables, written by `tools/randsctile-tables.py`, which
  builds each tiling and checks that twenty thousand points at random each
  fall in exactly one tile before writing it. Stars, crosses and the T are
  concave, and are placed by crossing number.
* The ones that never repeat are worked out point by point: de Bruijn's
  multigrid for the rhombus tilings, Robinson's triangles for the kites and
  darts, substitution for the pinwheel and the chair.
* The skew draws each tile's own outline shrunk, concave ones included, and
  the kaleidoscope folds them as it folds the rest of the family.
* On the machine it was written on it costs 255 to 380 ns a call for the
  periodic tilings and 410 to 680 for the ones that never repeat, against 220
  for `randscq`, and the two binaries draw all forty-five alike to the pixel.

## Palette 8: Kandinsky

The gradients of Kandinsky's paintings, **from mat to vivid**. Measured from
four vivid vignettes in his manner — black, vermilion, yellow, cerulean,
violet and pink on cream paper — and one mat picture of burnt orange, ochre,
brown, sage and petrol:

* Eighteen colours, each with its mat form, the mat picture's colour nearest
  it in lightness and hue: vermilion turns burnt orange, yellow ochre, black a
  dark brown, cerulean sage.
* Every palette draws a mood: four in five fall between mat and vivid, the
  fifth goes past the mat towards grey, and one palette in twenty or so is all
  but grey.
* What follows a colour is where his soft gradients lead it — paper into a
  wash of yellow, black into vermilion, blue into mauve into pink — counted in
  the pictures wherever the paint changes softly.
* One palette in ten has **anchors by accident**: any stop may turn to plain
  black or plain white, at a rate between one stop in twenty and one in ten,
  with no rhythm.

It stands 0.06 from the vignettes and 0.07 from the mat picture, at its two
ends, and 0.39 to 0.67 from the other seven palettes, as the Jensen–Shannon
divergence of the colours they show.

## Tests

The same 343 test programs, with new checks for every one of the forty-five
tilings — no point without a tile, near the origin or a million cells out,
values in range, flat tiles, one tile to the unit of area, no two tilings the
same field — and for the skew, the `selfsim` and the kaleidoscope over them.
And for palette 8, over four hundred palettes made as the program makes them:
dark and light where a picture looks in all but four, one colour or none in
twenty-six (the greyest, as asked, and a few single fields of colour), and an
anchor in thirty-five.

## Compatibility

* Positions saved by 1.7 load and draw as they did: palettes 1 to 7 are
  unchanged to the bit, and so are `randsc`, `randscq`, `randscp`, `randsch`
  and `randsct`.
* A position that uses palette 8 or `randsctile` names something 1.7 and the
  original XaoS do not have.
