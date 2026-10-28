# Notes for AI agents working on MS XaoS

How the translations and the tests of this fork are kept, as they have been
done so far. The mechanics of Qt's i18n are in [i18n/README.md](i18n/README.md);
the test targets are described where they are declared, in
[tests/CMakeLists.txt](tests/CMakeLists.txt). This file says how to work with
them and what to check before calling a change done.

## Ground rules

- Agree the plan with the user before editing, and fix the cause rather than
  the symptom.
- Commit and push only when asked. Never commit `CMakeLists.txt.user`; leave
  the untracked `doc/*.xpf` files and `doc/PROBLEMI_XPF/` alone.
- Change the version only when asked.
- A position saved before a change must draw exactly as it did: every new
  setting starts at a value that leaves the old code running as it was, and a
  position saves the setting only when it has moved, so that it still opens in
  the last release.
- Commit messages are English prose: a title that says what changed, then
  paragraphs on what and why. They end with the `Co-Authored-By` line the
  session gives.

## Building on Windows

```bash
export PATH=/c/Qt/Tools/CMake_64/bin:/c/Qt/Tools/mingw1310_64/bin:/c/Qt/Tools/Ninja:/c/Qt/6.12.0/mingw_64/bin:$PATH
cmake --build build/Desktop_Qt_6_12_0_MinGW_64_bit-Release
cmake --build build/Desktop_Qt_6_12_0_MinGW_64_bit-Release --target deploy
```

`bin/XaoS.exe` (64-bit mantissa) and `bin/XaoS-quad.exe` (113-bit) are what it
builds; they need `C:\Qt\6.12.0\mingw_64\bin` and `C:\Qt\Tools\mingw1310_64\bin`
on PATH to start. The build compiles `i18n/*.ts` into `.qm` files itself, so
there is nothing to run for them. Build with no warnings.

The build does not copy Qt's DLLs into `bin/`: the `deploy` target does, with
windeployqt, together with the image format plugins (`bin/imageformats/`,
where `qwebp.dll` is what saving WebP needs). Run it after changing Qt, or
`bin/` keeps the DLLs of the Qt it had before. The kit is Qt 6.12.0 since
version 2.16; the 6.8.0 build folders point at a Qt that is gone.

## Translations

### Which languages

German, Spanish, French, Italian, Portuguese and Russian (`de es fr it pt ru`)
are complete: every text of the interface and the whole guide. Every new or
changed text is translated into these six in the same change. The other
catalogues (`cs he hi hu is ro rs sk sv tr vi`) receive the new texts from
lupdate as unfinished and are left for a translator.

### What is marked, and what is not

- A literal shown on screen is `TR("Context", "text")`. A text that reaches
  the screen through a table, a variable or a macro argument is marked with
  `QT_TRANSLATE_NOOP("Context", ...)` where it is written: lupdate sees only
  literals inside the call.
- The choices of a dialog are always shown through `TR("Dialog", ...)`, so a
  choice table is marked `QT_TRANSLATE_NOOP("Dialog", ...)` to be translated,
  as `decomp_ways` is. The same table is what a position saves as keywords,
  in English, so a translation never changes a saved file.
- Not translated: the names of the colouring modes (`incolorname`,
  `outcolorname`; they are formulas, and words such as `squares` or
  `tessellation` stay with them), choices that are formulas (`decomp_of`:
  `z`, `z-c`, `z*c`), function names and argument lists in the formula
  reference, and commands and keywords.
- One source text in one context has one translation everywhere it is used.
  Before reusing a short word (`right`, `left`, `Speed:`), look at what it is
  already translated as and where.

### After adding or changing a text

1. Update the catalogues, from the top of the tree, exactly so:

   ```bash
   lupdate -tr-function-alias QT_TRANSLATE_NOOP+=TR src -ts i18n/XaoS_*.ts
   ```

   It reports the source texts it found for each catalogue; the count only
   grows by what was added. If it is far lower than before (run without the
   alias it found 344 of about 780 and dropped the rest; over `XaoS.pro` it
   misses the sources the .pro does not list), stop and restore the
   catalogues with `git checkout -- i18n`. Do not pass `-no-obsolete`.
2. List what is new: `grep -B2 'type="unfinished"' i18n/XaoS_it.ts`.
3. Before translating, read how the neighbouring texts are translated in each
   language (`grep -A1 '<source>Frequency (periodic modes):</source>'`) and
   keep their terms and form. The terms already in use:

   | English | it | de | es | fr | pt | ru |
   | --- | --- | --- | --- | --- | --- | --- |
   | Fractal | Frattale | Fraktal | Fractal | Fractale | Fractal | Фрактал |
   | Incoloring mode | Colorazione interna | Innere Färbung | Coloreado interior | Coloriage intérieur | Coloração interior | Внутренняя раскраска |
   | Outcoloring mode | Colorazione esterna | Äußere Färbung | Coloreado exterior | Coloriage extérieur | Coloração exterior | Внешняя раскраска |
   | Coloring speed | Velocità della colorazione | Färbungsgeschwindigkeit | Velocidad del coloreado | Vitesse de coloration | Velocidade da coloração | Скорость раскраски |
   | Decomposition | Decomposizione | Zerlegung | Descomposición | Décomposition | Decomposição | Разложение |
   | periodic modes | modalità periodiche | periodische Modi | modos periódicos | modes périodiques | modos periódicos | периодические режимы |

   French puts a space before `:` and `;`. Portuguese is European
   Portuguese.
4. Fill in the six catalogues with a script that replaces
   `<source>TEXT</source>\n        <translation type="unfinished"></translation>`
   with the translation, asserting one match each. The source is XML-escaped
   in the file (`&apos;`, `&amp;`). Write the files back with `newline=''`:
   they are LF in the working tree, and git's warning that it will make them
   CRLF is expected.
5. A text that was added and removed again before any release is left by
   lupdate as `type="vanished"` where it had been translated: delete those
   `<message>` blocks. Texts that were in a release stay as vanished.
6. Check: `grep -c 'type="unfinished"'` gives 0 for the six.
7. Look at it in the running program in a language with long words, Italian
   or German (see "In the running program" below for setting the language).

### The guide

`doc/ms-xaos-guide.md` is the reference; `doc/ms-xaos-guide.LL.md` are its
translations into the same six languages. A change to the English guide is
made to all six in the same commit. The line at the top of each translation
names the version of the English guide it follows; it changes only with the
version.

- Prose is wrapped at 79 columns, counted in characters (`awk length` counts
  bytes, which is wrong for Cyrillic and accents). Table rows are one line
  each, however long.
- Rewrap only the paragraphs you change. When rewrapping by script, keep each
  inline code span on one line (replace its spaces with a placeholder before
  `textwrap.fill`, `break_long_words=False, break_on_hyphens=False`, and put
  them back), and check that no rewrapped line starts with `- `, `* `, `+ `,
  `# `, `> ` or `1. `, which would turn it into a list or a heading.
- In a translation, menu paths and dialog fields are named as that language's
  interface names them (take them from the catalogue); mode names, commands
  and keywords stay English, in backticks. Links to a heading use the
  translated heading's anchor.
- Matching paragraphs across languages by line numbers does not work: they
  wrap differently. Join each paragraph into one line, find it by a phrase,
  replace, and rewrap.
- `doc/RELEASE-NOTES.md` is English only and is rewritten whole at each
  release.

## Tests

### Running them

```bash
cd build/Desktop_Qt_6_12_0_MinGW_64_bit-Release && ctest -j4
```

The suite has 435 tests; say how many when reporting. `randsc_speed` and
`phist_speed` measure ratios of timings and can fail under a loaded machine:
run them again on their own (`ctest -R randsc_speed`) before treating a
failure as real.

The tests, renders and GUI runs write to the registry (the formula history,
the language, the last settings). Back it up before and restore it after,
every time:

```bash
reg export "HKCU\\Software\\XaoS Project" '<scratchpad>\before.reg' //y
```

```powershell
reg delete "HKCU\Software\XaoS Project" /f; reg import "<scratchpad>\before.reg"
```

Never stop a XaoS the user has open; a script stops only the process it
started.

### How the suite is laid out

- Every suite is built once per precision, `use_long_double` (XaoS.exe) and
  `use_float128` (XaoS-quad.exe), with `xaos_add_test_binary`, and without Qt.
  A new test is added to the loop in `tests/CMakeLists.txt` the same way.
- The engine tests include `formulas.cpp` rather than link it, to reach its
  static functions, and stub the few things it calls. Copy the scaffolding of
  `tests/incolor_frequency_test.cpp` or `tests/decomposition_test.cpp`.
- `engine_unchanged` compares the iteration loops' checksums with
  `engine_ref.<precision>.expected`.
- `position_roundtrip` renders every `tests/positions/*.xpf`, saves it with
  `savepos`, loads what was saved and renders it again: the two images must be
  the same and every run must exit cleanly. Every new command a position saves
  gets a position here.

### Writing a test for a change to the colouring

- **What there was.** With a new setting where it starts, the old modes must
  give what they gave, to the bit. Keep the old code running unchanged at the
  default (a `decomp_plain`-style test that chooses the old branch) instead of
  hoping a general formula reproduces it, and pin checksums over a grid taken
  from that code. Check that both precisions give the same checksum; if they
  do not, pin each.
- **Exact equivalences, not tolerances.** Frequency 2 is the picture at z
  scaled by 2 (powers of two scale exactly); a rotation of a full turn is no
  rotation; an odd number of binary sectors is the even one above; a function
  applied in the engine equals the same arithmetic done on z in the test.
- **Each setting moves what it belongs to, and nothing else**: every mode it
  reaches changes, the other modes and the other side do not.
- **Grids** a half step off the axes when a boundary lies on them, so that two
  ways of cutting the same half-planes agree on every point.
- **The context in a test** is a zero-initialised global: set every field
  whose starting value is not nought (speeds, frequency, repeats) as
  `make_fractalc` does, or the test runs a branch the program never does.

### Checks beyond the suite

- **Draw it.** `XaoS.exe -render pos.xpf -size 320x240 -basename name` writes
  `name000000.png` in the working directory, and exits 1 on a command it
  cannot read. For a visual change, render the cases side by side on one
  sheet (PIL) and look at them. When the user reports a picture, render their
  own position first.
- **Against the last release.** `tools/ms-xaos-<version>-windows-x64.zip` is
  the published build: unpack it in the scratchpad, render the same positions
  with both binaries and compare the PNGs' checksums. This is how "draws as
  2.5 does, pixel for pixel" is shown.
- **Prototype first** when the user asks for a choice of looks: draw the
  candidates outside XaoS (numpy over a Mandelbrot grid) on one sheet, then
  implement the ones chosen.

### In the running program

- The language is the registry value
  `HKCU\Software\XaoS Project\XaoS\MainWindow\language` (`en_US`, `it_IT`, ...),
  read at start-up; set it after the backup.
- Drive the window with UI Automation and posted messages. A menu entry that
  opens a modal dialog never returns from Invoke: expand the menu, post
  VK_DOWN to the popup until the entry is highlighted, then VK_RETURN. In the
  dialog, ValuePattern sets the line edits; a Qt combo box opens its list as a
  separate top-level window of the process, so find the ListItem there and
  click it. Field labels can contain `=`.
- Confirm with OK and capture the window after, not just the dialog: the
  check is that the picture changed as it should.
