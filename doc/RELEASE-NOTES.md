# MS XaoS 1.8.0

What has changed since 1.7.4. Everything MS XaoS adds over
[XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the notes of
[1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7),
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2) and
[1.7.4](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.4), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v1.8.0/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## Long formulas

**A user formula of more than about five hundred characters crashed the
program** the moment it was set, and `-render` of a position holding one
crashed too. The formula is shown as a message once it is set, and the message
was copied into a 512-byte buffer on the stack to be measured and drawn; a
longer one wrote past it. The line now goes to Qt as it is, with no buffer and
no limit.

Two more limits of the same kind went with it:

* **A position could not hold a formula of more than about a thousand
  characters.** The reader stopped at 1024 characters a token with "Token is
  too long", so such a position was written but could not be read back — and
  undo and redo, which replay one, failed the same way. The token grows with
  what it holds now: a `poly` of 254 coefficients, 4587 characters, saves,
  loads and undoes.
* The status line of the formula and the error for an unknown command were
  written with an unbounded `sprintf`.

## The formula bar

Fractal → User formula has a bar that **goes on to a new line between names,
never inside one**: `randsctile`, `z_p3` and `0.6` are kept whole, and a line
may end after a bracket, a comma or an operator. It grows a line at a time,
and once the dialog would no longer fit a 720p screen it scrolls instead. The
dialog is about half as wide as such a screen. Return is still OK, and a pasted
formula has its line breaks made spaces.

**History...** opens the ten formulas used last, newest first and each once.
One too long for the window is listed as its start, three dots and its length —
`poly(z, {0.01097,0.00699}, {-0.05... (157 characters)` — and a double click or
Use puts it in the bar whole.

The history is written **when the dialog is accepted, and only then**. It was
written by the `usrform` command, which every position and every `-render`
holding a user formula runs as well, so loading or rendering a few pictures
pushed out the formulas you had typed.

**View → Status** shows a long formula the same way: on two lines, the second
under the first, and past two cut with its length, `... (918 characters)`. It
was one line that ran off the edge of the picture.

## Utilities: a random poly

The user formula reference has a sixth tab, **Utilities**, which writes a random
`poly` from a shortest and a longest length (1 to 254 coefficients) and a
smallest and a largest magnitude (0 and 1 when left empty). Each coefficient is
complex, at an angle drawn round the whole turn, written to four significant
figures of its size. There is no `c`: a user formula starts `z` at the pixel in
either mode, so what it draws is the polynomial's own filled Julia set, and `c`
goes wherever it is wanted. **Copy** puts it on the clipboard; while a field is
empty or wrong there is no formula, just a line saying why.

## Names that say what they are

The reference and the hint under the bar named arguments `a, b, c`, which said
nothing about which place is which. Every one is named for what it is for now —
`base, exponent`, `z, low, high`, `before, after, n`. `randsc`'s `mode` is
`kaled_mirror`, and the second argument of `min`, `max` and their variants is
`threshold`. These are names in the reference only: every formula reads as it
did.

## Colouring speed functions

Fractal → Incoloring and Outcoloring mode → Coloring speed functions has three
more, after the nine there were so that every number a saved position holds
still names what it named:

* **log log**, for the deep zooms where log alone still packs the bands too
  close — at a million iterations log gives 13.8 bands and this 2.7;
* **exp**, the other way from log, for the quantities that stay small;
* **triangle**, the sine with its curves straightened, the palette walked
  forwards and back at a steady pace.

The nine give what they gave, to the bit.

### Speed change after ITER

The last entry of the same menu changes whichever function is chosen: up to an
iteration you name it goes on exactly as it was, and past it the count it is
handed runs faster or slower, so the colours speed up or slow down from there
and not before. An exponent, an ITER and one of three ways of joining the two —
**jump**, **continuous** and **smooth** — a side, saved as
`(outcolorfunafter 40 1.5 'smooth)`. It starts off, and off it does nothing at
all. The guide draws the three and measures them.

### Past the end of an int

A colour value past what an int holds — 8.4 million bands — was converted
anyway, which is undefined, and every such pixel took one and the same palette
entry. `cube` gets there at 204 iterations, so a cube-coloured picture past that
was a flat patch. It now lands where its value says, wrapped round the palette
like any smaller value.

Two smaller repairs came with it: choosing a speed function moved the tick in
the inside colouring menu rather than its own, and **initstate** did not put the
fbm's five numbers back to their defaults, so a position loaded after one that
had moved them drew with the other's.

## Aspect ratio

**View → Aspect ratio** resizes the window once so that the image takes the
shape chosen — 1:1, 4:3, 3:2, 16:10, 16:9, 2:1, 21:9 or 9:16 — keeping the
pixels it had, and leaves it free to be dragged to any shape after. From full
screen or maximised it comes back to a window first. It is for the render: a
view framed in a 16:9 window renders at 1920×1080 as it was seen, where one
framed at 4:3 shows a third more at the sides.

The selection zoom and the fixed zoom steps moved into **View → Zoom**.

## Tests

347 tests, four more than 1.7.4: the colouring speed functions against what
the nine gave and what the new ones are for, the random `poly` against the
polynomial written out term by term, the wrapping of the bar, the history and
the status, and two more positions in the round trip — one using the speed
change after ITER, one holding a formula of 4587 characters, which crashed on
the render and could not be read back.

## Compatibility

* Positions saved by 1.7.4 load as they did and draw as they did, except where
  a colour value went past what an int holds: `cube` past 204 iterations, and
  smooth colouring of Magnet and Magnet2 far from the set, show the bands
  where they had one flat colour.
* Positions saved by 1.8.0 carry the speed change after ITER — off, unless it
  was used — and 1.7.4 does not know the command, so it cannot open them.
* A formula history written by 1.7.4 can hold the same formula several times;
  it is shown once, and cleared of the others the next time a formula is
  accepted.
