# MS XaoS 1.7

What has changed since 1.6. Version 1.6.1 was never released on its own, so
what it brought is here too. Everything MS XaoS adds over
[XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the
[notes of 1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6), and
the [guide](https://github.com/michele-summo/ms-xaos/blob/v1.7/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## The noise functions

**`selfsim`, the eighth argument of `randsc`, `randscq`, `randscp`, `randsch`
and `randsct`.** Each pass of these functions is a field of its own with cells
smaller than the pass before, so a formula that calls one on every pass reads a
finer field each time, and once the cells are smaller than a pixel the inside of
the picture is snow. With `selfsim` the call hands back instead the average of
all the passes so far, weighed as the octaves of a fractional Brownian motion:
pass *n* weighs d^(nH), where *d* is the degradation and *H* the argument. On a
hundred passes of `randsct(672,,{0.75,0.75},,,{2,2},6)+c`, the share of
neighbouring inside pixels a sixteenth of the palette or more apart goes from
70% to 0.6%.

* `1` is the plain motion, `0.5` rougher, larger smoother.
* It may be **complex**: the imaginary part turns every pass a little further
  than the one before, a spiral of octaves.
* It is off when the argument is **left out** or its place left empty; `0` is
  the plain average of the passes, which is what every value near nought gives.
* Passes a call did not see are worked out when it is asked, so the answer does
  not depend on the order the pixels are computed in.

**`skew_mode`, the seventh argument**, says what the skew draws inside a cell,
as numbers that add together: `1` the cell's own shape, `2` a rosette with the
kaleidoscope's symmetry, `4` a turn that follows where the point stands in its
wedge of the kaleidoscope, the same in every wedge so the kaleidoscope stays
one, and `8` the modulus as well as the angle — the only one `zmag` and the
bailout can see.

## Writing formulas

**Suffixes on a variable.** `z`, `c`, `x`, `p` and `p1` to `p9999` take
suffixes that stand for the calls made on them most often:

| suffix | stands for |
| --- | --- |
| `_b` `_bi` `_br` | `bship(…)` `bshipi(…)` `bshipr(…)` |
| `_pM` | `parchment(…, M)` |
| `_paM` | `parchmenta(…, M)` |

read from left to right: `c_b_p2` is `parchment(bship(c),2)`. They draw exactly
what the calls written out draw, and the formula is saved as it was written.

## Palettes

**Palettes 4 to 7 are new.** The four that had those numbers — a spectrum, a
duotone, a triad and a complementary pair — came out one colour each on the
screen: they spread their hues along the whole palette, which in the program
is some three thousand stops long, and a picture only ever sees the first few.
What replaced them was made after measuring what 1 to 3 are made of, twenty
thousand palettes of each, and each new one was measured against 1 to 3 in turn:

| | |
| --- | --- |
| 4 | **Smog** — dark and melancholy: greys, ochre, slate and dusty violet, with oxblood and scarlet among them |
| 5 | **Warm over night** — rose, red, brown, orange and yellow over teal, petrol, blue and indigo |
| 6 | **Favoured pairs** — the pairs of colours 1 to 3 put side by side most often, set down whole between black and white |
| 7 | **Random colours** — every stop a colour at random, and nothing else |

The palette dialog now says what each algorithm is called after its number.

One finding along the way, left as it is since saved positions depend on it:
palettes 1 and 2 hold only 256 colours, and every one of them is always
followed by the same one. That is where their familiar pairs of colours come
from.

## Tests

The same 343 test programs, with new checks for `selfsim` — against the
weighted average written out by hand, at both precisions — for the suffixes,
drawn through the engine's own loop, and for the palettes, now also made as the
program makes them. The suffixes were also fuzzed, while they were written, with
four million random formulas.

## Compatibility

* Positions saved by 1.6 load unchanged, **except** that one using palette 4 to
  7 comes back in the new colours: the number and the seed are all a position
  records of its palette, and what the number means has changed. Palettes 1 to
  3 are unchanged to the bit.
* A `randsc` call that names no `selfsim` and no `skew_mode` draws what it drew
  in 1.6.
* From the repository's 1.6.1, which was never released: a written
  `selfsim` of `0` now gives the plain average instead of switching it off, and
  `skew_mode` 4 draws differently.
