/* Generated from the sfcmplxfunc table in sffe_cmplx_gsl.cpp, whose
 * comments are the descriptions below. Plain data with no Qt in sight, so
 * that a test can link it beside the parser and check the two have not
 * drifted -- see tests/formula_help_test.cpp. Adding a function to the
 * parser without describing it here fails that test.
 *
 * Rows with no name are section headings, shown across both columns. */

#include "formulahelp.h"

const struct formula_help_row formula_help_functions[] = {
    {NULL, NULL, NULL, "operators, reached through sffe_operator rather than by name"},
    {"^", "", "a ^ b", NULL}, /* 2 */
    {"+", "", "a + b", NULL}, /* 2 */
    {"-", "", "a - b", NULL}, /* 2 */
    {"*", "", "a * b", NULL}, /* 2 */
    {"/", "", "a / b", NULL}, /* 2 */
    {"-", "", "-a", NULL}, /* 1 */
    {NULL, NULL, NULL, "trigonometry over the complex plane"},
    {"sin", "a", "nsin(a)", NULL}, /* 1 */
    {"cos", "a", "ncos(a)", NULL}, /* 1 */
    {"tan", "a", "ntan(a)", NULL}, /* 1 */
    {"cot", "a", "cot(a)", NULL}, /* 1 */
    {"asin", "a", "arcsin(a)", NULL}, /* 1 */
    {"acos", "a", "arccos(a)", NULL}, /* 1 */
    {"atan", "a", "arctan(a)", NULL}, /* 1 */
    {"acot", "a", "arccot(a)", NULL}, /* 1 */
    {"atan2", "a; b", "atan2 of the real parts, plus i times atan2 of the imaginary parts", NULL}, /* 2 */
    {NULL, NULL, NULL, "hyperbolic"},
    {"sinh", "a", "nsinh(a)", NULL}, /* 1 */
    {"cosh", "a", "ncosh(a)", NULL}, /* 1 */
    {"tanh", "a", "ntanh(a)", NULL}, /* 1 */
    {"coth", "a", "coth(a)", NULL}, /* 1 */
    {NULL, NULL, NULL, "exponential and logarithms"},
    {"exp", "a", "e ^ a", NULL}, /* 1 */
    {"log", "a", "natural logarithm of a", NULL}, /* 1 */
    {"log10", "a", "logarithm of a in base 10", NULL}, /* 1 */
    {"log2", "a", "logarithm of a in base 2", NULL}, /* 1 */
    {"logn", "a; b", "logarithm of b in base a -- the base comes first", NULL}, /* 2 */
    {NULL, NULL, NULL, "powers and roots"},
    {"pow", "a; b", "a ^ b", NULL}, /* 2 */
    {"powd", "a; b", "a ^ real(b); the imaginary part of b is ignored", NULL}, /* 2 */
    {"sqr", "a", "a * a. This used to raise a to itself: sqr(3) was 27", NULL}, /* 1 */
    {"sqrt", "a", "principal square root of a", NULL}, /* 1 */
    {"rtni", "a; b; c", "the c-th of the b b-th roots of a", NULL}, /* 3 */
    {"inv", "a", "1 / a", NULL}, /* 1 */
    {NULL, NULL, NULL, "rounding, and pulling a value apart"},
    {"ceil", "a", "ceiling of each component", NULL}, /* 1 */
    {"floor", "a", "floor of each component", NULL}, /* 1 */
    {"abs", "a", "|a|, as a real", NULL}, /* 1 */
    {"rabs", "a", "|real(a)|, as a real; the imaginary part is dropped", NULL}, /* 1 */
    {"re", "a", "real(a), as a real", NULL}, /* 1 */
    {"im", "a", "imag(a), as a real", NULL}, /* 1 */
    {"arg", "a", "angle of a, as a real", NULL}, /* 1 */
    {"mod", "a", "fractional part of each component", NULL}, /* 1 */
    {"conj", "a", "real(a) - i*imag(a); this used to swap them", NULL}, /* 1 */
    {NULL, NULL, NULL, "burning ship variants"},
    {"bship", "a", "|real(a)| + i*|imag(a)|", NULL}, /* 1 */
    {"bshipr", "a", "|real(a)| + i*imag(a)", NULL}, /* 1 */
    {"bshipi", "a", "real(a) + i*|imag(a)|", NULL}, /* 1 */
    {NULL, NULL, NULL, "assembling one value out of two"},
    {"rect", "a; b", "real(a) + i*imag(b)", NULL}, /* 2 */
    {"polar", "a; b", "|a| * e^(i*arg(b))", NULL}, /* 2 */
    {NULL, NULL, NULL, "smaller of two; the suffix says which part is compared"},
    {"min", "a; b", "smaller real part and smaller imaginary part", NULL}, /* 2 */
    {"minr", "a; b", "smaller real part; imaginary part taken from a", NULL}, /* 2 */
    {"mini", "a; b", "smaller imaginary part; real part taken from a", NULL}, /* 2 */
    {"minm", "a; b", "smaller modulus, held at the angle of a", NULL}, /* 2 */
    {NULL, NULL, NULL, "larger of two, same convention"},
    {"max", "a; b", "larger real part and larger imaginary part", NULL}, /* 2 */
    {"maxr", "a; b", "larger real part; imaginary part taken from a", NULL}, /* 2 */
    {"maxi", "a; b", "larger imaginary part; real part taken from a", NULL}, /* 2 */
    {"maxm", "a; b", "larger modulus, held at the angle of a", NULL}, /* 2 */
    {NULL, NULL, NULL, "mid(a, b, c) confines a to the range b..c."},
    {"mid", "a; b; c", "both components confined", NULL}, /* 3 */
    {"midr", "a; b; c", "real part confined; imaginary part taken from a", NULL}, /* 3 */
    {"midi", "a; b; c", "imaginary part confined; real part taken from a", NULL}, /* 3 */
    {"midm", "a; b; c", "modulus confined, held at the angle of a", NULL}, /* 3 */
    {NULL, NULL, NULL, "real trigonometry applied to each component separately"},
    {"sincos", "a", "nsin(real(a)) + i*ncos(imag(a))", NULL}, /* 1 */
    {"cossin", "a", "ncos(real(a)) + i*nsin(imag(a))", NULL}, /* 1 */
    {"sinr", "a", "nsin(real(a)); imaginary part passes through", NULL}, /* 1 */
    {"cosr", "a", "ncos(real(a)); imaginary part passes through", NULL}, /* 1 */
    {"sini", "a", "nsin(imag(a)); real part passes through", NULL}, /* 1 */
    {"cosi", "a", "ncos(imag(a)); real part passes through", NULL}, /* 1 */
    {"tancot", "a", "ntan(real(a)) + i*cot(imag(a))", NULL}, /* 1 */
    {"cottan", "a", "cot(real(a)) + i*ntan(imag(a))", NULL}, /* 1 */
    {"tanr", "a", "ntan(real(a)); imaginary part passes through", NULL}, /* 1 */
    {"cotr", "a", "cot(real(a)); imaginary part passes through", NULL}, /* 1 */
    {"tani", "a", "ntan(imag(a)); real part passes through", NULL}, /* 1 */
    {"coti", "a", "cot(imag(a)); real part passes through", NULL}, /* 1 */
    {NULL, NULL, NULL, "waveforms, again component by component"},
    {"trunc", "a", "each component truncated towards zero", NULL}, /* 1 */
    {"sawtooth", "a", "x - nfloor(x), a ramp in [0, 1)", NULL}, /* 1 */
    {"twave", "a", "triangle wave of period 2, in [-1, 1]", NULL}, /* 1 */
    {NULL, NULL, NULL, "assorted"},
    {"julian", "a; [b=1]; [c=1]", "|a|^b * e^(i*c*arg(a))", NULL}, /* 3 */
    {"inveps", "a; [b=0.01+0.01i]", "an inverse softened by b, so that it stays finite at the origin", NULL}, /* 2 */
    {"atan2s", "a; b", "atan2 of each pair of components; differs from atan2 on real arguments only", NULL}, /* 2 */
    {"ngon", "a; [b=0]; [c=3]; [d=1]", "folds a about the centre b onto a c-sided polygon, corner radius raised to d", NULL}, /* 4 */
    {"parchment", "a; b", "quantises the angle of a into |b| sectors, keeping |a|", NULL}, /* 2 */
    {"parchmenta", "a; b", "as parchment, but mirroring alternate half sectors", NULL}, /* 2 */
    {NULL, NULL, NULL, "snapping to a grid of step 1/n; n == 0 leaves the value alone"},
    {"truncv", "a; b", "both components, step 1/|b|", NULL}, /* 2 */
    {"truncc", "a; b", "real step 1/real(b), imaginary step 1/imag(b)", NULL}, /* 2 */
    {"truncvr", "a; b", "real component only, step 1/|b|", NULL}, /* 2 */
    {"truncvi", "a; b", "imaginary component only, step 1/|b|", NULL}, /* 2 */
    {"truncvm", "a; b", "the modulus, angle kept", NULL}, /* 2 */
    {"truncva", "a; b", "the angle, modulus kept", NULL}, /* 2 */
    {NULL, NULL, NULL, "special functions, none of them elementary"},
    {"erf", "a", "error function over the complex plane", NULL}, /* 1 */
    {"gamma", "a", "the gamma function over the complex plane", NULL}, /* 1 */
    {"lambertw", "a", "principal branch of the Lambert W of a", NULL}, /* 1 */
    {NULL, NULL, NULL, "choosing by iteration"},
    {"ifiter", "a; b; ...", "picks one argument by iteration number, cycling; only that one runs", NULL}, /* variadic */
    {"ifiterl", "a; b; ...", "as ifiter, but holding the last argument once the iterations run out", NULL}, /* variadic */
    {"ifiterf", "a; b", "b on the last iteration the limit allows, a on the rest; only the chosen one runs", NULL}, /* 2 */
    {"ifiterr", "a; b; n", "a while the iteration is below n, b from n on; only the chosen one runs", NULL}, /* 3 */
    {NULL, NULL, NULL, "randomness: all but the seed optional, and a new field every iteration; what the kaleidoscope and its mode do is in the Values tab"},
    {"rand", "a", "real(a) times a random number in [0, 1); depends on call order, so a redraw differs", NULL},
    {"randsc", "seed; [size=1+i]; [degradation=0.5+0.5i]; [kaleidoscope=1]; [mode=0]", "coherent noise over the point: soft blobs, size wide and size high", NULL},
    {"randscq", "seed; [size=1+i]; [degradation=0.5+0.5i]; [kaleidoscope=1]; [mode=0]", "the same field with no interpolation: a mosaic of flat square cells", NULL},
    {"randscp", "seed; [size=1+i]; [degradation=0.5+0.5i]; [kaleidoscope=1]; [mode=0]", "the same field cut into irregular flat polygons, with straight edges", NULL},
    {"randsch", "seed; [size=1+i]; [degradation=0.5+0.5i]; [kaleidoscope=1]; [mode=0]", "the same field cut into hexagons: a honeycomb of flat cells", NULL},
    {"randsct", "seed; [size=1+i]; [degradation=0.5+0.5i]; [kaleidoscope=1]; [mode=0]", "the same field cut into equilateral triangles, alternating in orientation", NULL},
    {NULL, NULL, NULL, "watching the orbit: both hand back their argument until the last iteration, and what they gathered on it, which the inside colouring modes then draw"},
    {"trap", "a; [shape=0]; [centre=0]; [size=1]", "how near the orbit ever came to a shape; the shapes are in the Values tab", NULL},
    {"stripe", "a; [density=4]", "the average of (sin(density*arg a)+1)/2 along the orbit; density a whole number, how many stripes go round a turn", NULL},
    {NULL, NULL, NULL, "polynomials"},
    {"poly", "z; k1; k2; ...", "k1*z^(m-1) + k2*z^(m-2) + ... + km: the first coefficient written multiplies the highest power, the last stands alone", NULL},
    {NULL, NULL, NULL, NULL}};

/* The numbers that appear as arguments and mean something particular.
 *
 * They were written into the description of the function that takes them,
 * which made those descriptions long enough to burst the column and left the
 * numbers themselves undocumented anywhere one would think to look. A tab of
 * their own is where one looks. Not compared against the parser -- there is
 * nothing in the parser to compare it with -- so it is kept by hand.
 */
const struct formula_help_row formula_help_values[] = {
    {NULL, NULL, NULL, "trap: which shape the orbit is measured against"},
    {"0", "", "the centre itself, a point", NULL},
    {"1", "", "a horizontal line through the centre", NULL},
    {"2", "", "a vertical line through the centre", NULL},
    {"3", "", "both of them, a cross", NULL},
    {"4", "", "a ring of radius size", NULL},
    {"5", "", "a square of half-side size", NULL},
    {"6", "", "a diamond of half-diagonal size", NULL},
    {NULL, NULL, NULL, "randsc family: how many wedges the kaleidoscope folds the plane into"},
    {"1", "", "no folding at all, which is what a call that says nothing gets", NULL},
    {"2 or more", "", "that many wedges around the origin, the field taken from one", NULL},
    {NULL, NULL, NULL, "randsc family: which mirror the kaleidoscope folds with"},
    {"0", "", "the far half of each wedge mirrors the near half, so a wedge is symmetric about its bisector", NULL},
    {"1", "", "the same the other way about, the near half mirroring the far one", NULL},
    {"anything else", "", "folds the way 0 does", NULL},
    {NULL, NULL, NULL, NULL}};

/* The variables the engine registers before parsing a user formula; see the
 * sffe_regvar calls in formulas.cpp. Anything else named in a formula is a
 * parse error, which is what keeps the two lists worth writing down. */
const struct formula_help_row formula_help_variables[] = {
    {NULL, NULL, NULL, "the iteration"},
    {"z", "", "the current value, and what the formula computes the next one from", NULL},
    {"c", "", "the point being iterated -- the pixel, in the mandelbrot sense", NULL},
    {"n", "", "the iteration number, counting from zero", NULL},
    {"x", "", "scratch value, kept across iterations", NULL},
    {NULL, NULL, NULL, "earlier values of z; Fractal -> Set p values on first iteration says what they are before there are any"},
    {"p1", "", "the value z had on the pass before this one", NULL},
    {"p2", "", "the one before that, and so on", NULL},
    {"p9999", "", "as far back as one may look; any p up to this is a variable, whether or not the picture runs that many passes", NULL},
    {"p", "", "another name for p1", NULL},
    {NULL, NULL, NULL, NULL}};

/* How values are written. These are properties of the parser rather than
 * named things it holds, so they are listed apart from the variables. */
const struct formula_help_row formula_help_notation[] = {
    {NULL, NULL, NULL, "writing a value"},
    {"{re,im}", "", "a complex number given by its two parts, as in {0,2} for 2i", NULL},
    {"{re;im}", "", "the same, with a semicolon; both separators are accepted", NULL},
    {"i", "", "the imaginary unit, as the shipped formulas use it", NULL},
    {"1.5", "", "a real number, read at the full precision of this build", NULL},
    {NULL, NULL, NULL, "calling"},
    {"f(a; b)", "", "arguments are separated by a semicolon or a comma", NULL},
    {"[b=1]", "", "an argument shown in square brackets may be left out; the value shown is what the function uses instead", NULL},
    {"f(a; ;c)", "", "one may be left out in the middle as well, by leaving its place empty; spaces make no difference", NULL},
    {"poly(z;1; ;1)", "", "an empty coefficient of poly is a term that is not there: this is z^2 + 1", NULL},
    {"ifiter(a; ;b)", "", "an empty branch of ifiter or ifiterl repeats the one before it: two iterations of a, then b", NULL},
    {"2z", "", "multiplication may be left out before a name or a bracket", NULL},
    {"(a)(b)", "", "the same, between two bracketed groups", NULL},
    {"-z^2", "", "a leading minus binds looser than a power: it means -(z^2)", NULL},
    {NULL, NULL, NULL, NULL}};
