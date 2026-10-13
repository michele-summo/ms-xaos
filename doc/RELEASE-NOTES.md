# MS XaoS 2.1

What has changed since 1.8.0. Version 2.0 was never released on its own, so
what it brought is here too. Everything MS XaoS adds over
[XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the notes of
[1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7),
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2),
[1.7.4](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.4) and
[1.8.0](https://github.com/michele-summo/ms-xaos/releases/tag/v1.8.0), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v2.1/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

The first number has moved because one thing saved by 1.8.0 no longer draws to
the pixel what it drew: the speed change after ITER in **smooth** or
**continuous**. Those two modes are gone, and such a position is read into the
ones there are now, near enough to look alike.

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

Of the 169 examples rendered with one thread, 167 draw to the pixel as 1.8.0
drew them: `galaxy` has its stars from the clock, and `heart` keeps a message
that 1.8.0 had taken away before the picture was drawn. With four threads they
all draw as with one, `galaxy` apart.

## Tests

347 tests, as many as 1.8.0, with more in them. The colouring speed test now
checks each stretch at its own exponent and speed; smooth with no corner and
never turning back, the jump's count past its radius, and at ITER the pace half
way between two steady ones; groups that are off, and speeds at one, changing
nothing to the bit; and 1.8.0's smooth and continuous against the curves they
drew. The round
trip has two more positions: one with both groups, a radius and the speeds, and
one written as 1.8.0 and the groups without a radius saved it. The threads were
checked by rendering, as above.

## Compatibility

* Positions saved by 1.8.0 load. A speed change after ITER that is off or in
  jump draws as it did; in smooth or continuous it draws alike, not to the
  pixel.
* Positions saved by 2.1 open in 1.8.0 and draw the same there when their
  speed change is off or one group in jump, with the initial exponent and the
  speeds at one. Any other speed change is written with commands 1.8.0 does not
  know, so it cannot open them.
* Positions saved by 2.0 open; their smooth is drawn with the rounding of 2.1.
* Pictures drawn with one thread are 1.8.0's to the pixel, but for a message a
  render now times by its own clock; with several threads they are now the
  same as with one.
