# MS XaoS 2.3

What has changed since 1.8.0. Versions 2.0, 2.1 and 2.2 were never released
on their own, so what they brought is here too. Everything MS XaoS adds over
[XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the notes of
[1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7),
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2),
[1.7.4](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.4) and
[1.8.0](https://github.com/michele-summo/ms-xaos/releases/tag/v1.8.0), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v2.3/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

The first number has moved because one thing saved by 1.8.0 no longer draws to
the pixel what it drew: the speed change after ITER in **smooth** or
**continuous**. Those two modes are gone, and such a position is read into the
ones there are now, near enough to look alike. And the colour count is now
worked out whole and at an even pace, which moves the colours of a smooth
picture by a few levels and of some others by one.

## The speed change after ITER

Fractal → Incoloring and Outcoloring mode → Coloring speed functions → **Speed
change after ITER** was an ITER, an exponent and one of four modes a side. It
is now an **initial exponent** and an **initial speed**, and two groups, each
an **ITER**, an **exponent**, a **speed** and a **mode** — off, jump or smooth
— with the **radius** of smooth, which the dialog shows only when smooth is
chosen.

**Each stretch runs at its own exponent.** The count is the plain power `v^p0`
up to the first ITER of a group that is on; each stretch after it starts from
the count the one before came to, `f(ITER)`, and adds
`s·ITER·((v/ITER)^p − 1)`, whose pace at ITER is `s·p`. With `p0` = 2, `p1` =
1 at 50 and `p2` = 2 at 100, the count on `ident` is `v²` up to 50, `v + 2450`
up to 100, and `2550 + 100·((v/100)² − 1)` from there — never `p0·p1`. The
lower ITER goes first whichever group holds it, and a group that is off is left
out, its ITER with it.

**The speeds are linear where the exponents are powers**: each multiplies the
pace of its stretch, so at an exponent of one a stretch is the straight line
`f(ITER) + s·(v − ITER)`, and a speed of 3 on the first group above makes the
middle stretch `2500 + 3·(v − 50)`.

It all starts off, with the exponent and every speed at one, where nothing at
all is done to the value: every function gives to the bit what it gave before.

### Smooth, the same on both sides of ITER

**Jump** turns the pace at ITER at once, which shows as a ring where the bands
change pace. **Smooth** turns it from `ITER − r` to `ITER + r`: the pace goes in
a straight line from the one the count has at one end to the one it has at the
other, with a parabola about ITER that gives back what the line cut off the
jump's corner. It is symmetric about ITER, and it **changes the colours inside
the radius and nowhere else** — before `ITER − r` and past `ITER + r` the
picture is the jump's, to the digit. After `v²` up to 50 and an exponent of
one, a radius of 20 takes the pace from 60 at 30 to 45.5 at ITER and 1 at 70,
and goes on from 70 as `v + 2450`, as jump does.

The parabola never takes the pace below the slower of the two ends. Where it
would — a slow pace turned very fast over a short radius — it is held there,
and the count past the rounding moves on by what that keeps back. A radius
wider than the first ITER stops short of nought, and two that would overlap are
cut down alike.

**Continuous is gone.** It went on from ITER at the pace the count came with,
to join the stretches without a ring, and smooth with a wide radius does that.

### Positions

| saved as | read as |
| --- | --- |
| off, jump | the same |
| 1.8.0's smooth | smooth at a radius of three quarters of ITER: within four hundredths of ITER of what it drew up to an exponent of one and a half, a tenth up to two |
| continuous | smooth at three quarters of ITER, with the exponent `0.66 + p/4 + 0.09·p²`: within twelve hundredths of ITER of it up to a half |

Found by setting the curves side by side over four ITERs: such a picture comes
out alike, not to the pixel. A speed change that is no more than one group in
jump, with the initial exponent and every speed at one, is saved as 1.8.0 saved
it, `(outcolorfunafter 40 1.5 'jump)`, and draws the same there. Anything else
is `(outcolorspeedchange 1 40 1.5 'smooth 10 50 2 'off 10)`, and the speeds,
only where one of them is not one, `(outcolorspeedlinear 1 3 1)`. The dialog
itself is the command `outcolorspeedgroups`, its twelve fields in their order.

A dialog field can now be shown for one choice of another and hidden for the
rest, the dialog fitting itself to what it shows: that is how the radius is
there in smooth alone.

## The same picture with any number of threads

A picture drawn with several threads came out a few pixels different now and
then: `examples/Hubicka/julia3.xpf`, rendered with four threads on a busy
machine, differed in 26 of 3000 renders by up to 37 pixels along its
filaments, and never with one thread.

It came from **solid guessing**, which fills a pixel from the rows and columns
already calculated round it when they agree. Its passes ran in all the threads
at once with nothing between them, so whether a neighbour counted as calculated
depended on how far the other threads had got. Each pass is now two phases, its
rows and then its columns, and every thread finishes one before any starts the
next. The 3000 renders come out alike, and with one thread the order is the one
it always was.

Three more things changed with the threads or with the machine:

* **3D** draws a band of columns to a thread, and the last two columns of a
  band read the first two of the band to its right — which each band used to
  fill with copies of its own last column. A band now waits, row by row, for
  what the band on its right has drawn there, and the bands still run side by
  side.
* **Starfield** drew its stars from one random generator that every thread
  advanced at once. Each row now starts where one thread drawing the whole
  picture would have been. The stars still change from one start of XaoS to
  the next: the generator is seeded from the clock, as it always was.
* **Messages in a render**, such as the formula's name, went after five
  seconds of the wall's clock, so a render slowed by a busy machine drew them
  in fewer frames, or in none — which looked like antialiasing changing the
  picture. In a render they now last five seconds of the animation's own time;
  in the window nothing changes.

With this alone, of the 169 examples rendered with one thread, 167 drew to the
pixel as 1.8.0 drew them: `galaxy` has its stars from the clock, and `heart`
keeps a message that 1.8.0 had taken away before the picture was drawn. With
four threads they all draw as with one, `galaxy` apart. The colour count below
moves more of them, for reasons of its own.

## The colour count

The value a colouring mode hands to the palette lost digits on the way, and
smooth colouring worked out its fraction of an iteration at an uneven pace.
Neither showed at a colouring speed of one; both did as soon as a speed change
or a colouring speed ran the colours fast.

### Whole

Smooth colouring made its count an int of 256ths of an iteration before the
colouring speed, the function or the speed change had it, and what the
function gave back was a float, which past 65 536 entries keeps less than a
256th of one. With a speed change at 2500 entries an iteration, a view came
out in **stripes of one colour** eight pixels wide, ten entries apart, where
the palette should have run through every entry between. The count now goes on
whole, and the functions, the colouring speed and the speed change are worked
in the build's own type, as the formulas are — long double in the ordinary
build, `__float128` in quad.

Two things are kept as they were:

* **An infinity** — the middle of a potential gone to infinity, the ground
  round a formula that escapes at once — lands on the first entry, where the
  float's overflow always sent it. A finite value past what a float holds
  lands where it says.
* **A value on the edge of a 256th** stays on it: the value is rounded to a
  256th of a 256th before it is cut to a 256th. The real part of a Newton
  orbit comes to rest on a root, exactly on such an edge, and the last digits
  of the convergence put each pixel on one side or the other — a 256th of an
  entry, until an **edge filter** draws every such pixel and the contours come
  out toothed like gears. The float rounded most of them away and left
  notches; now there are none.

**Magnet and Magnet2 in smooth colouring** change the most: far from the set
their count is negative, which the int took four thousand million 256ths on
and the float drew in steps of an entry or two, flat patches with hard edges
where there is now the gradient the count describes.

### At an even pace

Smooth colouring took the fraction of the last iteration where the bailout
falls between the orbit's squared modulus before the escaping pass and after
it, linearly in its logarithm. The count ran on unbroken that way but its pace
did not: the modulus is all but squared from one pass to the next, so the pace
changed across every iteration and **turned at once at its edge** — by 4.26 in
a view of a speed change running ten thousand entries an iteration, which drew
the bands changing density at a stroke.

With the circle, the orbit is now carried on past the bailout until its
squared modulus passes a million — three or four passes more for a square —
and the fraction is taken in the logarithm of the logarithm, where a squaring
is a step of one length, with the degree the orbit's own. The count stays
within a tenth of an iteration of what it was, and the turn of 4.26 comes out
at 1.0000.

That is for the Mandelbrots, Triceratops, Mandelbar, Lambda, Manowar, Spider,
Phoenix and Octo, and for a user formula while its orbit grows as a power. The
Barnsleys multiply the orbit rather than square it, and the count at the
bailout is already the right one for them; Newton, Magnet and Catseye converge
or test their escape their own way; and the other bailout shapes measure the
orbit by something that turns with it, which carried on made the turns worse.
Those keep the count at the bailout.

### What it costs

On a picture that comes out the same, the count worked whole costs two to five
per cent in the ordinary build and seven to twenty in quad, where a logarithm
in `__float128` is worked in software. The passes and logarithms of the even
pace cost 11 to 28 per cent more in the ordinary build and 18 to 20 in quad on
the smooth pictures tried, and nothing on the others. Where stripes become a
gradient the render takes longer as well, since solid guessing fills a pixel
only from neighbours of one colour.

## Undo

Edit → **Undo** (`U`, or `Ctrl` `Z`) now goes back one step for anything that
changes the picture, and **Redo** (`Ctrl` `Y` on Windows) forward again. It
went wrong in two ways:

* **After the mouse** it went back to a picture half way through the zoom,
  and a second Undo was needed to get to where it had started: the step was
  kept at the release, while the zoom was still slowing down. It is now kept
  where the view comes to rest.
* **After a command that changes nothing** — Recalculate, hiding the messages
  — it did nothing at all. Every command kept a step, and undo went two back
  on the understanding that it had kept the last one itself, so it brought
  back the picture already on screen.

A picture is now kept only where it differs from the last one kept, and undo
and redo move one step along them whatever was kept in between. A step is a
command or a dialog — a formula, a palette, a colouring mode, a fixed zoom —
a press of the mouse with the zoom, pan or rotation it sets going, a
selection zoom, the iterations or the cycling speed set with the arrows, and
an animation played to its end. A new step after an undo cuts off what could
have been redone, as in an editor. Undo stops at the picture XaoS opened
with, after the configuration and the command line, and what moved on its
own since the last step — the autopilot, a continuous rotation — is kept when
Undo is pressed, so that Redo comes back to it. The last 256 pictures are
kept.

Each step is a saved position, so undo needs a position to come back as it
was saved. The 179 among the examples and the tests do, to the letter,
loaded alone or over one another. One more did not: a position with a
built-in palette, loaded after a custom palette, was saved again as a custom
palette, which undo took for a change and which cost the redo. A built-in
palette now clears the custom one, as a generated one always did.

## Tests

349 tests, two more than 1.8.0 and more in the others. The colouring speed test
now checks each stretch at its own exponent and speed; smooth with no corner
and never turning back, the jump's count past its radius, and at ITER the pace
half way between two steady ones; groups that are off, and speeds at one,
changing nothing to the bit; and 1.8.0's smooth and continuous against the
curves they drew. It checks the colour count too: a hundred and one colours
along a fortieth of a 256th of an iteration at a high speed, where the int gave
one; a value given back to a quarter past ten million; a value on the edge of a
256th staying on it; an infinity on the first entry; and, through the
Mandelbrot's own smooth loop, steps of the same length either side of an
iteration's edge, 8.9 entries and 8.8, where they were 8.3 and 33.5. The
fingerprints of the iteration loops were taken again in the smooth columns of
the formulas whose count changed; everywhere else they are as they were. The
round trip has two more positions: one with both groups, a radius and the
speeds, and one written as 1.8.0 and the groups without a radius saved it. The
threads were checked by rendering, as above.

The two new tests are the undo history, at each precision: a step only where
the picture differs, undo and redo one step at a time and stopping at either
end, the redo kept when the picture undo went to is kept again and cut off by
a new one, and the oldest step dropped from a full history. Undo itself was
tried in the running program — a palette, a formula, a zoom with the mouse, a
Recalculate, the iterations with the arrows, an undo while the zoom was still
slowing down — each undone and redone and set against the picture it should
come back to, pixel for pixel.

## Compatibility

* Positions saved by 1.8.0 load. A speed change after ITER that is off or in
  jump draws as it did; in smooth or continuous it draws alike, not to the
  pixel.
* Positions saved by 2.3 open in 1.8.0 when their speed change is off or one
  group in jump, with the initial exponent and the speeds at one. Any other
  speed change is written with commands 1.8.0 does not know, so it cannot open
  them.
* Positions saved by 2.0, 2.1 and 2.2 open; the smooth of 2.0's speed change
  is drawn with the rounding of 2.1.
* A position with a built-in palette, opened after a custom palette, is saved
  again with its own palette, where 1.8.0 saved it as a custom one.
* Of the 169 examples rendered one at a time at 640×480, 82 draw to the pixel
  as 1.8.0 drew them and 57 more differ by three levels at most in at most 113
  pixels. Of the other thirty, four differ by a message or by the clock, as
  above; the fifteen that colour smoothly with a power, by up to thirteen
  levels across much of the picture — `eye` by up to 39, and `spiral4`, whose
  smooth log on a true-colour palette magnifies any change, in the middle of
  one small detail; three of true colour and smooth log by up to five; `edge`
  has lost its notches; and seven differ in at most a hundred pixels, where a
  value sat on the edge of an entry and a filter or an unblended palette shows
  the step.
* With several threads a picture is now the same as with one, and the two
  builds colour as alike as they did.
