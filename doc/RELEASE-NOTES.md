# MS XaoS 2.5

What has changed since
[2.4](https://github.com/michele-summo/ms-xaos/releases/tag/v2.4). Everything
MS XaoS adds over [XaoS](https://github.com/xaos-project/XaoS) 4.3.3 is in the
notes of [1.6](https://github.com/michele-summo/ms-xaos/releases/tag/v1.6),
[1.7](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7),
[1.7.2](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.2),
[1.7.4](https://github.com/michele-summo/ms-xaos/releases/tag/v1.7.4),
[1.8.0](https://github.com/michele-summo/ms-xaos/releases/tag/v1.8.0),
[2.3](https://github.com/michele-summo/ms-xaos/releases/tag/v2.3) and
[2.4](https://github.com/michele-summo/ms-xaos/releases/tag/v2.4), and the
[guide](https://github.com/michele-summo/ms-xaos/blob/v2.5/doc/ms-xaos-guide.md)
says what each thing does and why it was made that way.

## Coloring speed in one dialog

Incoloring mode and Outcoloring mode had four entries for what takes a
colouring value to a place in the palette: Coloring speed, Coloring shift, a
submenu of the twelve speed functions and, at its foot, the speed change
after ITER in a dialog of its own. Each side now has one, **Coloring speed**:
the function chosen from a list, the speed and the shift, and under the
heading *Speed change after ITER* the fields its own dialog had. OK sets them
all as one step of undo, each as it would be set alone, so only what moved is
redrawn.

## How often the inside modes repeat

Six inside modes repeat at a pace written into them: cos(mag),
mag\*cos(real^2), sin(real^2-imag^2), squares, sin(real)\*sin(imag) and
frac(mag). The colouring speed acts on the value they hand over, not inside
them, so nothing set how close together their bands lay.

**Frequency (periodic modes)**, the first field of the inside Coloring speed
dialog, multiplies what they repeat on: at 2 the squares are an eightieth of
a unit on a side instead of a fortieth, and frac(mag) draws sixteen contours
to a unit of the modulus instead of eight. At one, where it starts and where
New puts it back, the six draw as they did to the bit. The inside modes that
do not repeat never read it.

## Fractional Brownian Motion in one dialog

It was a submenu of Other coloring mode, with its three modes as radio
entries and a Settings dialog under them for the five numbers. It is now
**Fractional Brownian Motion...** in Incoloring mode and Outcoloring mode,
beside the modes themselves, and ticked whenever one of its three is in use.
The dialog chooses the mode from a list and has the five numbers below it;
OK sets them together as one step of undo.

## Six languages complete

German, Spanish, French, Italian, Portuguese and Russian now have every text
of the interface translated. Much of what reached the screen through a table
or a variable had never been offered to the translators: the names of the
filters and of the palette algorithms, the bailout shapes, the true-colour
modes, the choices of every dialog, the descriptions and headings of the user
formula reference, and the submenus, which were named before there was a
language and so kept their English. The formula reference leaves the names
and argument lists as they are, since they are syntax, and a saved position
still names things in English. The other languages have the new texts
waiting for a translator.

OK, Cancel and the file and colour dialogs are Qt's and come from its own
translations. XaoS looked for them where Qt is installed, which a machine
running XaoS out of its folder has not got, and loaded them in the system's
language whatever language had been chosen. Those of the six languages are
now built in, from Qt 6.8.0 (Portuguese is Qt's Brazilian, there being no
other), and the language chosen is the one loaded.

The guide is in the six languages too, each saying at its top which version
of the English guide it was made from and that the English one is the
reference. Help, and What is new in About, open the guide in the language
XaoS is in, and the English one where there is none.

## A crash after New

XaoS could crash minutes after a position had been loaded and File → New
chosen, while drawing text. Loading a position, which Open, Undo and
`-loadpos` all do, freed the context that plays it and kept the pointer to
it; New wrote zeros into that freed block, and whatever the heap gave it to
next was damaged. The pointer is now cleared with the free, which puts right
the morph dialog as well, the other place that looked at it.

## Memory through a session

Every saved position left two formula parsers behind, and the undo history
saves a position before and after every command, so XaoS took more memory
the longer it was used. They are freed with the rest.

## Removed

* The OpenGL driver. `-DOPENGL=ON`, off by default, drew the picture with
  glDrawPixels from a mirrored copy of the image every frame and left out the
  full-screen handling of the menu bar. Nothing depended on it; the option,
  its CI job and its mentions in About and the guide go with it.
* Code that nothing calls: `src/formeval`, a separate formula evaluator of
  2020 that neither build compiled, and functions left with no caller, among
  them the line and rectangle drawing of the image code and the complex
  functions no formula reaches. Help no longer looks for a help file in
  `HELPPATH`; it has opened the guide on line for a long time.

## Tests

353 tests, two more than 2.4: the frequency at each precision. At one, the
six modes it reaches must give the checksums they gave before it existed;
away from one, for the four that read z alone, a frequency must draw the
same as z scaled by it, exactly; and the modes that do not repeat must never
read it. The round trip of saved positions, each drawn, saved, loaded and
drawn again, has one that sets the frequency.

The two dialogs were tried in the running program, from their menus to the
picture after OK.

## Compatibility

* Positions saved by 2.4 load. Those saved by 2.5 open in 2.4 unless they
  set the frequency, which a position saves only when it is not one.
* A position writes the colouring speed and Fractional Brownian Motion as
  before. The commands of the entries the dialogs replaced are kept, out of
  the menu, for a script: `mincolorspeed`, `mincolorshift`, `infun0` to
  `infun11`, `incolorspeedgroups` and the same outside, and `infbmset` and
  `outfbmset`.
* The radio entries of the old Fractional Brownian Motion submenu, `in23` to
  `in25` and `out18` to `out20`, are gone; `(incoloring 23)` and the like do
  what they did.
* A build configured with `-DOPENGL=ON` gets the ordinary widget.
