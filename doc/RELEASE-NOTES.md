# MS XaoS 2.4

What has changed since
[2.3](https://github.com/michele-summo/ms-xaos/releases/tag/v2.3). Everything
MS XaoS adds over [XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the
notes of [1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7),
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2),
[1.7.4](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.4),
[1.8.0](https://github.com/michele-summo/ms-xaos/releases/tag/v1.8.0) and
[2.3](https://github.com/michele-summo/ms-xaos/releases/tag/v2.3), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v2.4/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## Turning and centring the view

Fractal → **View** sets the centre, radius and angle by number, and now has
two tools under the angle. They fill in its fields and do nothing else: OK
applies what they show as one step of undo, Cancel none of it, and a saved
position never hears of them.

* **Slices** cuts a full turn into that many equal ones, and **−** and **+**
  turn the angle back or on by one: 60 degrees for six, 360 for one. With
  none, the two are disabled. The turns are counted from the angle they
  started from, so six sixths are a whole turn to the last digit, and typing
  an angle or changing the slices starts the count again. The number of
  slices is kept while XaoS runs.
* **Center by axis** moves the centre straight across the picture (**R**) or
  straight up it (**I**) until the origin is in the middle that way, at the
  angle in the dialog. Unturned, that is the real or the imaginary part
  becoming nought and the other kept. Turned, it is the same along the
  picture's own axes: R at 90 degrees clears the imaginary part, and R and
  then I bring the origin to the middle at any angle. A part the turning
  leaves as rounding is taken for nought, so the field says 0 and not 1e-17.

## Turned views in the quad build

The view was turned by the sine and cosine of a double pi, which in the quad
build also meant the double functions: a pair that makes a turn only to
sixteen digits, so the centre of a turned view was drawn 5e-17 of its
distance from the origin away from where it was. Down to 1e-16 that is out
of sight. Below it, where the quad build is for, a view turned by 30 degrees
showed another place altogether. The turn is now worked in the build's own
type, pi and all, and is exact on a quarter turn.

It costs nothing that can be measured: it is worked once a frame, not once a
pixel. Rendered five times each with one thread, before and after, turned
and unturned pictures of both builds took the same time to within two and a
half per cent either way, which is the spread between runs. The deep turned
view of the quad build took 15 per cent less, but that is another picture,
since 2.3 drew another place.

## Closing after the palette editor

XaoS could crash on closing, about one time in four after the palette editor
had been used. The main window freed the context that keeps the image and
then read the image out of it, as XaoS has done since before this fork;
mostly the freed memory still held it. The image is now taken out first.
Closed forty-two times since, twenty-eight of them after the palette editor,
neither build crashed.

## Tests

351 tests, two more than 2.3: the View dialog's tools at each precision. They
check n slices of n coming round to a whole turn for every n up to 360; the
axes at every quarter turn, either way round and past a whole turn, with no
-0; the origin left in the middle of the picture as the engine turns it, at
eight other angles; the same axis twice changing nothing, and R then I giving
the origin exactly. And the turn itself: its sine and cosine making a turn to
the last digits, and a point turned and turned back being where it was, which
the quad build of 2.3 fails.

The dialog was tried in the running program: the slices, an angle typed by
hand, both axes, OK, Cancel and undo, and the picture after OK set against a
second XaoS opened on the numbers the dialog showed, pixel for pixel. Closing
was tried as above.

## Compatibility

* Positions saved by 2.3 load, and those saved by 2.4 open in 2.3.
* A turned picture in the ordinary build draws as it did but for scattered
  pixels where the detail is finer than the arithmetic. Of the 64 examples
  with an angle, 52 draw to the pixel as 2.3 drew them; `paisleypillow`,
  turned at a depth of 1.2e-16, moves a fifth of its pixels, `rpg` one in
  twenty, and the other ten less than one in a hundred.
* A position saved by the quad build of 2.3, turned and deeper than about
  1e-16, opens where its numbers say, which is not where 2.3 showed it.
