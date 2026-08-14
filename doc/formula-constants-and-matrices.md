# Constant folding and constant matrices in user formulas

Design note, not yet implemented. Written 2026-08-14, from an assessment of
the SFFE parser as it stands after the precision and correctness work on this
branch.

## What is wanted

Two additions to the user-defined formula language, both concerning values
that are fixed when the formula is parsed.

**Constant folding.** A sub-expression built only from literals — `2*pi/3`,
`(1+i)^7` — is evaluated once while parsing and replaced by its result,
instead of being recomputed on every iteration.

**Constant matrices.** A matrix literal written with an `M` prefix:

    M[[a, b, c, ...], [...], ...]

The elements are constant expressions. Variables — `z`, `c`, `n`, `x`, `p`,
`p1` and the rest the engine registers — are **not** allowed inside them.

Matrices are usable only by functions written specifically to take them, which
do not exist yet. Every other interaction is an error, reported at parse time:
a product or sum of two matrices, a matrix where a scalar is expected, a scalar
where a matrix is expected. No symbolic computation of any kind is intended.

## Why the `M` prefix matters

`[` and `]` are already taken. PHASE 2 of the parser normalises `[` to `{` and
`]` to `}` (`src/sffe/sffe.cpp`, around the bracket handling), so brackets are
an alternative spelling of parentheses and `[[1,2],[3,4]]` parses today as
nested groups.

With the prefix, the tokeniser sees `M` as an identifier and can hand the rest
to a matrix scanner before that normalisation reaches the inner brackets. That
turns the main syntactic obstacle into a local change. Without it the choice
would have been different delimiters or a context-sensitive tokeniser.

## What the constraint buys

Because everything is constant and variables are refused, matrices and folded
expressions need to exist **only at parse time**. Every element and every
reference collapses to a stored value before evaluation begins.

The runtime evaluator does not change. The flat operation list, one complex
number per operand, the lazy dispatch — all of it stays as it is, which is the
hot and delicate part. A folded sub-expression leaves `sfset()` looking exactly
like a literal the user typed.

## Difficulty, piece by piece

| piece | difficulty | note |
| --- | --- | --- |
| `M[[...],...]` literal | low | contained in PHASE 2; the prefix avoids the collision |
| constant table | low | new, isolated data structure |
| refusing variables inside `M[]` | low | the registered names are already known to the parser |
| type tracking in the shunting-yard | medium-high | touches `pop_expression` |
| constant folding | medium | needs a purity flag on all 96 entries of the function table |
| a matrix type inside `sfarg` | medium | struct change, ripples to every constructor |
| semantic errors on mixed types | low | falls out of the type tracking |

## The two real risks

**`pop_expression` is the spine.** It carries the operand counting, the
variadic argument handling and the lazy dispatch, all reworked on this branch.
Adding type tracking there is where a regression would hide most comfortably.

**Folding must not touch iteration-dependent functions.** `ifiter` reads
`sffe_iteration`, which changes on every pass. Folding it at parse time yields
a formula that looks right and draws the wrong fractal — failure with no
symptom. A purity flag per table entry is mechanical to add and silent to get
wrong.

The mitigation is already in place: 119 parser cases, the engine checksum net,
the precision tests and the position round trip would all catch a regression.
That is what makes this reasonable to attempt.

## Compatibility

`M[...]` lives inside `usrform`, which is already a string, so no new `.xpf`
directive is needed and no saved file becomes unreadable. An older version
opening such a file reports a syntax error in the formula, which is the honest
outcome rather than a silently different picture.

Contrast `(precision 113)`, which had to be a command precisely because
refusing the file was the point.

## Suggested order

Do them separately, folding first.

Folding stands on its own: any formula with constant sub-expressions gains,
and the iteration loop is the hottest code in the program. It is also the
lower-risk half — no struct changes, no new types — and its correctness is
checkable with what already exists: identical results, fewer operations per
iteration, both measurable.

Matrices are the expensive half: the type inside `sfarg`, the tracking through
the parser, the ripple to every constructor. Doing them on top of a folding
that is already verified avoids having two unknowns open at once.

## One decision to make first

Folding moves when some errors surface. `1/0` inside a constant sub-expression
becomes a **parse error** rather than a NaN at run time.

That is probably wanted — it says so immediately instead of much later — but it
changes the behaviour of formulas that are accepted today, so it is a choice
rather than a detail.
