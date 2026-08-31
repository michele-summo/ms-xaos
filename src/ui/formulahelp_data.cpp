/* Generated from the sfcmplxfunc table in sffe_cmplx_gsl.cpp, whose
 * comments are the descriptions below. Plain data with no Qt in sight, so
 * that a test can link it beside the parser and check the two have not
 * drifted -- see tests/formula_help_test.cpp. Adding a function to the
 * parser without describing it here fails that test.
 *
 * Rows with no name are section headings, shown across both columns. */

#include "formulahelp.h"

const struct formula_help_row formula_help_functions[] = {
    {NULL, NULL, "operators, reached through sffe_operator rather than by name"},
    {"^", "a ^ b", NULL}, /* 2 */
    {"+", "a + b", NULL}, /* 2 */
    {"-", "a - b", NULL}, /* 2 */
    {"*", "a * b", NULL}, /* 2 */
    {"/", "a / b", NULL}, /* 2 */
    {"-", "-a", NULL}, /* 1 */
    {NULL, NULL, "trigonometry over the complex plane"},
    {"sin", "nsin(a)", NULL}, /* 1 */
    {"cos", "ncos(a)", NULL}, /* 1 */
    {"tan", "ntan(a)", NULL}, /* 1 */
    {"cot", "cot(a)", NULL}, /* 1 */
    {"asin", "arcsin(a)", NULL}, /* 1 */
    {"acos", "arccos(a)", NULL}, /* 1 */
    {"atan", "arctan(a)", NULL}, /* 1 */
    {"acot", "arccot(a)", NULL}, /* 1 */
    {"atan2", "atan2 of the real parts, plus i times atan2 of the imaginary parts", NULL}, /* 2 */
    {NULL, NULL, "hyperbolic"},
    {"sinh", "nsinh(a)", NULL}, /* 1 */
    {"cosh", "ncosh(a)", NULL}, /* 1 */
    {"tanh", "ntanh(a)", NULL}, /* 1 */
    {"coth", "coth(a)", NULL}, /* 1 */
    {NULL, NULL, "exponential and logarithms"},
    {"exp", "e ^ a", NULL}, /* 1 */
    {"log", "natural logarithm of a", NULL}, /* 1 */
    {"log10", "logarithm of a in base 10", NULL}, /* 1 */
    {"log2", "logarithm of a in base 2", NULL}, /* 1 */
    {"logn", "logarithm of b in base a -- the base comes first", NULL}, /* 2 */
    {NULL, NULL, "powers and roots"},
    {"pow", "a ^ b", NULL}, /* 2 */
    {"powd", "a ^ real(b); the imaginary part of b is ignored", NULL}, /* 2 */
    {"sqr", "a * a. This used to raise a to itself: sqr(3) was 27", NULL}, /* 1 */
    {"sqrt", "principal square root of a", NULL}, /* 1 */
    {"rtni", "the c-th of the b b-th roots of a", NULL}, /* 3 */
    {"inv", "1 / a", NULL}, /* 1 */
    {NULL, NULL, "rounding, and pulling a value apart"},
    {"ceil", "ceiling of each component", NULL}, /* 1 */
    {"floor", "floor of each component", NULL}, /* 1 */
    {"abs", "|a|, as a real", NULL}, /* 1 */
    {"rabs", "|real(a)|, as a real; the imaginary part is dropped", NULL}, /* 1 */
    {"re", "real(a), as a real", NULL}, /* 1 */
    {"im", "imag(a), as a real", NULL}, /* 1 */
    {"arg", "angle of a, as a real", NULL}, /* 1 */
    {"mod", "fractional part of each component", NULL}, /* 1 */
    {"conj", "real(a) - i*imag(a); this used to swap them", NULL}, /* 1 */
    {NULL, NULL, "burning ship variants"},
    {"bship", "|real(a)| + i*|imag(a)|", NULL}, /* 1 */
    {"bshipr", "|real(a)| + i*imag(a)", NULL}, /* 1 */
    {"bshipi", "real(a) + i*|imag(a)|", NULL}, /* 1 */
    {NULL, NULL, "assembling one value out of two"},
    {"rect", "real(a) + i*imag(b)", NULL}, /* 2 */
    {"polar", "|a| * e^(i*arg(b))", NULL}, /* 2 */
    {NULL, NULL, "smaller of two; the suffix says which part is compared"},
    {"min", "smaller real part and smaller imaginary part", NULL}, /* 2 */
    {"minr", "smaller real part; imaginary part taken from a", NULL}, /* 2 */
    {"mini", "smaller imaginary part; real part taken from a", NULL}, /* 2 */
    {"minm", "smaller modulus, held at the angle of a", NULL}, /* 2 */
    {NULL, NULL, "larger of two, same convention"},
    {"max", "larger real part and larger imaginary part", NULL}, /* 2 */
    {"maxr", "larger real part; imaginary part taken from a", NULL}, /* 2 */
    {"maxi", "larger imaginary part; real part taken from a", NULL}, /* 2 */
    {"maxm", "larger modulus, held at the angle of a", NULL}, /* 2 */
    {NULL, NULL, "mid(a, b, c) confines a to the range b..c."},
    {"mid", "both components confined", NULL}, /* 3 */
    {"midr", "real part confined; imaginary part taken from a", NULL}, /* 3 */
    {"midi", "imaginary part confined; real part taken from a", NULL}, /* 3 */
    {"midm", "modulus confined, held at the angle of a", NULL}, /* 3 */
    {NULL, NULL, "real trigonometry applied to each component separately"},
    {"sincos", "nsin(real(a)) + i*ncos(imag(a))", NULL}, /* 1 */
    {"cossin", "ncos(real(a)) + i*nsin(imag(a))", NULL}, /* 1 */
    {"sinr", "nsin(real(a)); imaginary part passes through", NULL}, /* 1 */
    {"cosr", "ncos(real(a)); imaginary part passes through", NULL}, /* 1 */
    {"sini", "nsin(imag(a)); real part passes through", NULL}, /* 1 */
    {"cosi", "ncos(imag(a)); real part passes through", NULL}, /* 1 */
    {"tancot", "ntan(real(a)) + i*cot(imag(a))", NULL}, /* 1 */
    {"cottan", "cot(real(a)) + i*ntan(imag(a))", NULL}, /* 1 */
    {"tanr", "ntan(real(a)); imaginary part passes through", NULL}, /* 1 */
    {"cotr", "cot(real(a)); imaginary part passes through", NULL}, /* 1 */
    {"tani", "ntan(imag(a)); real part passes through", NULL}, /* 1 */
    {"coti", "cot(imag(a)); real part passes through", NULL}, /* 1 */
    {NULL, NULL, "waveforms, again component by component"},
    {"trunc", "each component truncated towards zero", NULL}, /* 1 */
    {"sawtooth", "x - nfloor(x), a ramp in [0, 1)", NULL}, /* 1 */
    {"twave", "triangle wave of period 2, in [-1, 1]", NULL}, /* 1 */
    {NULL, NULL, "assorted"},
    {"julian", "|a|^b * e^(i*c*arg(a))", NULL}, /* 3 */
    {"inveps", "an inverse softened by b, so that it stays finite at the origin", NULL}, /* 2 */
    {"atan2s", "atan2 of each pair of components; differs from atan2 on real arguments only", NULL}, /* 2 */
    {"ngon", "folds a about the centre b onto a c-sided polygon, corner radius raised to d", NULL}, /* 4 */
    {"parchment", "quantises the angle of a into |b| sectors, keeping |a|", NULL}, /* 2 */
    {"parchmenta", "as parchment, but mirroring alternate half sectors", NULL}, /* 2 */
    {NULL, NULL, "snapping to a grid of step 1/n; n == 0 leaves the value alone"},
    {"truncv", "both components, step 1/|b|", NULL}, /* 2 */
    {"truncc", "real step 1/real(b), imaginary step 1/imag(b)", NULL}, /* 2 */
    {"truncvr", "real component only, step 1/|b|", NULL}, /* 2 */
    {"truncvi", "imaginary component only, step 1/|b|", NULL}, /* 2 */
    {"truncvm", "the modulus, angle kept", NULL}, /* 2 */
    {"truncva", "the angle, modulus kept", NULL}, /* 2 */
    {"erf", "error function over the complex plane", NULL}, /* 1 */
    {"gamma", "the gamma function over the complex plane", NULL}, /* 1 */
    {"lambertw", "principal branch of the Lambert W of a", NULL}, /* 1 */
    {"ifiter", "picks one argument by iteration number, cycling; only that one is evaluated", NULL}, /* variadic */
    {"ifiterl", "as ifiter, but holding the last argument once the iterations pass the end", NULL}, /* variadic */
    {"rand", "real(a) times a random number in [0, 1)", NULL},
    {"randsc", "coherent noise over the point: randsc(seed; size; degradation), the last two optional", NULL},
    {"randscq", "the same field without interpolation: a mosaic of flat square cells", NULL}, /* 1 */
    {NULL, NULL, NULL}};

/* The variables the engine registers before parsing a user formula; see the
 * sffe_regvar calls in formulas.cpp. Anything else named in a formula is a
 * parse error, which is what keeps the two lists worth writing down. */
const struct formula_help_row formula_help_variables[] = {
    {NULL, NULL, "the iteration"},
    {"z", "the current value, and what the formula computes the next one from", NULL},
    {"c", "the point being iterated -- the pixel, in the mandelbrot sense", NULL},
    {"n", "the iteration number, counting from zero", NULL},
    {"x", "scratch value, kept across iterations", NULL},
    {NULL, NULL, "parameters, set from Fractal -> Parameters"},
    {"p", "the first parameter; the same value as p1", NULL},
    {"p1", "first parameter", NULL},
    {"p2", "second parameter", NULL},
    {"p3", "third parameter", NULL},
    {"p4", "fourth parameter", NULL},
    {"p5", "fifth parameter", NULL},
    {"p6", "sixth parameter", NULL},
    {NULL, NULL, NULL}};

/* How values are written. These are properties of the parser rather than
 * named things it holds, so they are listed apart from the variables. */
const struct formula_help_row formula_help_notation[] = {
    {NULL, NULL, "writing a value"},
    {"{re,im}", "a complex number given by its two parts, as in {0,2} for 2i", NULL},
    {"{re;im}", "the same, with a semicolon; both separators are accepted", NULL},
    {"i", "the imaginary unit, as the shipped formulas use it", NULL},
    {"1.5", "a real number, read at the full precision of this build", NULL},
    {NULL, NULL, "calling"},
    {"f(a; b)", "arguments are separated by a semicolon or a comma", NULL},
    {"2z", "multiplication may be left out before a name or a bracket", NULL},
    {"(a)(b)", "the same, between two bracketed groups", NULL},
    {"-z^2", "a leading minus binds looser than a power: it means -(z^2)", NULL},
    {NULL, NULL, NULL}};
