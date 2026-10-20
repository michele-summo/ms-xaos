/* Written from the sfcmplxfunc table in sffe_cmplx_gsl.cpp. The arguments
 * are named for what they are -- base and exponent, low and high, rest and
 * last -- and not a, b and c, since the line under the formula bar shows
 * these names while a call is being written, and a letter says nothing about
 * which place is which. The operators keep a and b, being only the sides of
 * an expression.
 *
 * Plain data with no Qt in sight, so that a test can link it beside the
 * parser and check the two have not drifted -- see
 * tests/formula_help_test.cpp. Adding a function to the parser without
 * describing it here fails that test.
 *
 * Rows with no name are section headings, shown across both columns. */

#include <cstring>

#include "formulahelp.h"

const struct formula_help_row formula_help_functions[] = {
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "operators, reached through sffe_operator rather than by name")},
    {"^", "", QT_TRANSLATE_NOOP("FormulaHelp", "a ^ b"), NULL}, /* 2 */
    {"+", "", QT_TRANSLATE_NOOP("FormulaHelp", "a + b"), NULL}, /* 2 */
    {"-", "", QT_TRANSLATE_NOOP("FormulaHelp", "a - b"), NULL}, /* 2 */
    {"*", "", QT_TRANSLATE_NOOP("FormulaHelp", "a * b"), NULL}, /* 2 */
    {"/", "", QT_TRANSLATE_NOOP("FormulaHelp", "a / b"), NULL}, /* 2 */
    {"-", "", QT_TRANSLATE_NOOP("FormulaHelp", "-a"), NULL}, /* 1 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "trigonometry over the complex plane")},
    {"sin", "z", QT_TRANSLATE_NOOP("FormulaHelp", "sine of z"), NULL}, /* 1 */
    {"cos", "z", QT_TRANSLATE_NOOP("FormulaHelp", "cosine of z"), NULL}, /* 1 */
    {"tan", "z", QT_TRANSLATE_NOOP("FormulaHelp", "tangent of z"), NULL}, /* 1 */
    {"cot", "z", QT_TRANSLATE_NOOP("FormulaHelp", "cotangent of z"), NULL}, /* 1 */
    {"asin", "z", QT_TRANSLATE_NOOP("FormulaHelp", "arcsine of z"), NULL}, /* 1 */
    {"acos", "z", QT_TRANSLATE_NOOP("FormulaHelp", "arccosine of z"), NULL}, /* 1 */
    {"atan", "z", QT_TRANSLATE_NOOP("FormulaHelp", "arctangent of z"), NULL}, /* 1 */
    {"acot", "z", QT_TRANSLATE_NOOP("FormulaHelp", "arccotangent of z"), NULL}, /* 1 */
    {"atan2", "y, x", QT_TRANSLATE_NOOP("FormulaHelp", "atan2(real(y), real(x)) + i*atan2(imag(y), imag(x))"), NULL}, /* 2 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "hyperbolic")},
    {"sinh", "z", QT_TRANSLATE_NOOP("FormulaHelp", "hyperbolic sine of z"), NULL}, /* 1 */
    {"cosh", "z", QT_TRANSLATE_NOOP("FormulaHelp", "hyperbolic cosine of z"), NULL}, /* 1 */
    {"tanh", "z", QT_TRANSLATE_NOOP("FormulaHelp", "hyperbolic tangent of z"), NULL}, /* 1 */
    {"coth", "z", QT_TRANSLATE_NOOP("FormulaHelp", "hyperbolic cotangent of z"), NULL}, /* 1 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "exponential and logarithms")},
    {"exp", "z", QT_TRANSLATE_NOOP("FormulaHelp", "e ^ z"), NULL}, /* 1 */
    {"log", "z", QT_TRANSLATE_NOOP("FormulaHelp", "natural logarithm of z"), NULL}, /* 1 */
    {"log10", "z", QT_TRANSLATE_NOOP("FormulaHelp", "logarithm of z in base 10"), NULL}, /* 1 */
    {"log2", "z", QT_TRANSLATE_NOOP("FormulaHelp", "logarithm of z in base 2"), NULL}, /* 1 */
    {"logn", "base, z", QT_TRANSLATE_NOOP("FormulaHelp", "logarithm of z in the base given -- the base comes first"), NULL}, /* 2 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "powers and roots")},
    {"pow", "base, exponent", QT_TRANSLATE_NOOP("FormulaHelp", "base ^ exponent"), NULL}, /* 2 */
    {"powd", "base, exponent", QT_TRANSLATE_NOOP("FormulaHelp", "base ^ real(exponent); the imaginary part of the exponent is ignored"), NULL}, /* 2 */
    {"sqr", "z", QT_TRANSLATE_NOOP("FormulaHelp", "z * z. This used to raise z to itself: sqr(3) was 27"), NULL}, /* 1 */
    {"sqrt", "z", QT_TRANSLATE_NOOP("FormulaHelp", "principal square root of z"), NULL}, /* 1 */
    {"rtni", "z, n, k", QT_TRANSLATE_NOOP("FormulaHelp", "the k-th of the n n-th roots of z, counting from 0 at the principal one"), NULL}, /* 3 */
    {"inv", "z", QT_TRANSLATE_NOOP("FormulaHelp", "1 / z"), NULL}, /* 1 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "rounding, and pulling a value apart")},
    {"ceil", "z", QT_TRANSLATE_NOOP("FormulaHelp", "ceiling of each component"), NULL}, /* 1 */
    {"floor", "z", QT_TRANSLATE_NOOP("FormulaHelp", "floor of each component"), NULL}, /* 1 */
    {"abs", "z", QT_TRANSLATE_NOOP("FormulaHelp", "|z|, as a real"), NULL}, /* 1 */
    {"rabs", "z", QT_TRANSLATE_NOOP("FormulaHelp", "|real(z)|, as a real; the imaginary part is dropped"), NULL}, /* 1 */
    {"re", "z", QT_TRANSLATE_NOOP("FormulaHelp", "real(z), as a real"), NULL}, /* 1 */
    {"im", "z", QT_TRANSLATE_NOOP("FormulaHelp", "imag(z), as a real"), NULL}, /* 1 */
    {"arg", "z", QT_TRANSLATE_NOOP("FormulaHelp", "angle of z, as a real"), NULL}, /* 1 */
    {"mod", "z", QT_TRANSLATE_NOOP("FormulaHelp", "fractional part of each component"), NULL}, /* 1 */
    {"conj", "z", QT_TRANSLATE_NOOP("FormulaHelp", "real(z) - i*imag(z); this used to swap them"), NULL}, /* 1 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "burning ship variants")},
    {"bship", "z", QT_TRANSLATE_NOOP("FormulaHelp", "|real(z)| + i*|imag(z)|"), NULL}, /* 1 */
    {"bshipr", "z", QT_TRANSLATE_NOOP("FormulaHelp", "|real(z)| + i*imag(z)"), NULL}, /* 1 */
    {"bshipi", "z", QT_TRANSLATE_NOOP("FormulaHelp", "real(z) + i*|imag(z)|"), NULL}, /* 1 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "assembling one value out of two")},
    {"rect", "real_from, imag_from", QT_TRANSLATE_NOOP("FormulaHelp", "real(real_from) + i*imag(imag_from): the real part of the first and the imaginary part of the second"), NULL}, /* 2 */
    {"polar", "length, direction", QT_TRANSLATE_NOOP("FormulaHelp", "|length| * e^(i*arg(direction)): as long as the first, pointing where the second points"), NULL}, /* 2 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "smaller of two; the suffix says which part is compared")},
    {"min", "z, threshold", QT_TRANSLATE_NOOP("FormulaHelp", "the smaller real part and the smaller imaginary part of the two"), NULL}, /* 2 */
    {"minr", "z, threshold", QT_TRANSLATE_NOOP("FormulaHelp", "the smaller real part; imaginary part taken from z"), NULL}, /* 2 */
    {"mini", "z, threshold", QT_TRANSLATE_NOOP("FormulaHelp", "the smaller imaginary part; real part taken from z"), NULL}, /* 2 */
    {"minm", "z, threshold", QT_TRANSLATE_NOOP("FormulaHelp", "the smaller modulus, held at the angle of z"), NULL}, /* 2 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "larger of two, same convention")},
    {"max", "z, threshold", QT_TRANSLATE_NOOP("FormulaHelp", "the larger real part and the larger imaginary part of the two"), NULL}, /* 2 */
    {"maxr", "z, threshold", QT_TRANSLATE_NOOP("FormulaHelp", "the larger real part; imaginary part taken from z"), NULL}, /* 2 */
    {"maxi", "z, threshold", QT_TRANSLATE_NOOP("FormulaHelp", "the larger imaginary part; real part taken from z"), NULL}, /* 2 */
    {"maxm", "z, threshold", QT_TRANSLATE_NOOP("FormulaHelp", "the larger modulus, held at the angle of z"), NULL}, /* 2 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "smallest or largest modulus of any number")},
    {"minmc", "z1, z2, ...", QT_TRANSLATE_NOOP("FormulaHelp", "the smallest modulus of all the arguments, held at the angle of z1; with two, minm"), NULL}, /* variadic */
    {"maxmc", "z1, z2, ...", QT_TRANSLATE_NOOP("FormulaHelp", "the largest modulus of all the arguments, held at the angle of z1; with two, maxm"), NULL}, /* variadic */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "mid(z, low, high) confines z to the range low..high.")},
    {"mid", "z, low, high", QT_TRANSLATE_NOOP("FormulaHelp", "both components confined"), NULL}, /* 3 */
    {"midr", "z, low, high", QT_TRANSLATE_NOOP("FormulaHelp", "real part confined; imaginary part taken from z"), NULL}, /* 3 */
    {"midi", "z, low, high", QT_TRANSLATE_NOOP("FormulaHelp", "imaginary part confined; real part taken from z"), NULL}, /* 3 */
    {"midm", "z, low, high", QT_TRANSLATE_NOOP("FormulaHelp", "modulus confined, held at the angle of z"), NULL}, /* 3 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "real trigonometry applied to each component separately")},
    {"sincos", "z", QT_TRANSLATE_NOOP("FormulaHelp", "sin(real(z)) + i*cos(imag(z))"), NULL}, /* 1 */
    {"cossin", "z", QT_TRANSLATE_NOOP("FormulaHelp", "cos(real(z)) + i*sin(imag(z))"), NULL}, /* 1 */
    {"sinr", "z", QT_TRANSLATE_NOOP("FormulaHelp", "sin(real(z)); imaginary part passes through"), NULL}, /* 1 */
    {"cosr", "z", QT_TRANSLATE_NOOP("FormulaHelp", "cos(real(z)); imaginary part passes through"), NULL}, /* 1 */
    {"sini", "z", QT_TRANSLATE_NOOP("FormulaHelp", "sin(imag(z)); real part passes through"), NULL}, /* 1 */
    {"cosi", "z", QT_TRANSLATE_NOOP("FormulaHelp", "cos(imag(z)); real part passes through"), NULL}, /* 1 */
    {"tancot", "z", QT_TRANSLATE_NOOP("FormulaHelp", "tan(real(z)) + i*cot(imag(z))"), NULL}, /* 1 */
    {"cottan", "z", QT_TRANSLATE_NOOP("FormulaHelp", "cot(real(z)) + i*tan(imag(z))"), NULL}, /* 1 */
    {"tanr", "z", QT_TRANSLATE_NOOP("FormulaHelp", "tan(real(z)); imaginary part passes through"), NULL}, /* 1 */
    {"cotr", "z", QT_TRANSLATE_NOOP("FormulaHelp", "cot(real(z)); imaginary part passes through"), NULL}, /* 1 */
    {"tani", "z", QT_TRANSLATE_NOOP("FormulaHelp", "tan(imag(z)); real part passes through"), NULL}, /* 1 */
    {"coti", "z", QT_TRANSLATE_NOOP("FormulaHelp", "cot(imag(z)); real part passes through"), NULL}, /* 1 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "waveforms, again component by component")},
    {"trunc", "z", QT_TRANSLATE_NOOP("FormulaHelp", "each component truncated towards zero"), NULL}, /* 1 */
    {"sawtooth", "z", QT_TRANSLATE_NOOP("FormulaHelp", "each component less its floor, a ramp in [0, 1)"), NULL}, /* 1 */
    {"twave", "z", QT_TRANSLATE_NOOP("FormulaHelp", "triangle wave of period 2, in [-1, 1], on each component"), NULL}, /* 1 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "assorted")},
    {"julian", "z, [power=1], [angle_factor=1]", QT_TRANSLATE_NOOP("FormulaHelp", "|z|^power * e^(i*angle_factor*arg(z)): the modulus raised to power, the angle multiplied by angle_factor"), NULL}, /* 3 */
    {"inveps", "z, [softening=0.01+0.01i]", QT_TRANSLATE_NOOP("FormulaHelp", "1/z softened so that it stays finite at the origin: real(z)/(|z|^2 + real(softening)) - i*imag(z)/(|z|^2 + imag(softening))"), NULL}, /* 2 */
    {"atan2s", "y, x", QT_TRANSLATE_NOOP("FormulaHelp", "atan2 of each pair of components; differs from atan2 on real arguments only"), NULL}, /* 2 */
    {"ngon", "z, [centre=0], [sides=3], [power=1]", QT_TRANSLATE_NOOP("FormulaHelp", "folds z about centre onto a polygon with that many sides, the corner radius raised to power"), NULL}, /* 4 */
    {"parchment", "z, sectors", QT_TRANSLATE_NOOP("FormulaHelp", "quantises the angle of z into |sectors| sectors, keeping |z|"), NULL}, /* 2 */
    {"parchmenta", "z, sectors", QT_TRANSLATE_NOOP("FormulaHelp", "as parchment, but mirroring alternate half sectors"), NULL}, /* 2 */
    {"rotate", "z, angle", QT_TRANSLATE_NOOP("FormulaHelp", "z turned clockwise about the origin by angle radians, the real part of it: z*e^(-i*angle). Whatever the view's own angle is; a picture drawn from a point put through it turns clockwise with it, as fbm(rotate(x,0.5),7)"), NULL}, /* 2 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "snapping to a grid of step 1/n; n == 0 leaves the value alone")},
    {"truncv", "z, n", QT_TRANSLATE_NOOP("FormulaHelp", "both components, step 1/|n|"), NULL}, /* 2 */
    {"truncc", "z, n", QT_TRANSLATE_NOOP("FormulaHelp", "real step 1/real(n), imaginary step 1/imag(n)"), NULL}, /* 2 */
    {"truncvr", "z, n", QT_TRANSLATE_NOOP("FormulaHelp", "real component only, step 1/|n|"), NULL}, /* 2 */
    {"truncvi", "z, n", QT_TRANSLATE_NOOP("FormulaHelp", "imaginary component only, step 1/|n|"), NULL}, /* 2 */
    {"truncvm", "z, n", QT_TRANSLATE_NOOP("FormulaHelp", "the modulus, step 1/|n|, angle kept"), NULL}, /* 2 */
    {"truncva", "z, n", QT_TRANSLATE_NOOP("FormulaHelp", "the angle, step 1/|n| in radians, modulus kept"), NULL}, /* 2 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "special functions, none of them elementary")},
    {"erf", "z", QT_TRANSLATE_NOOP("FormulaHelp", "error function over the complex plane"), NULL}, /* 1 */
    {"gamma", "z", QT_TRANSLATE_NOOP("FormulaHelp", "the gamma function over the complex plane"), NULL}, /* 1 */
    {"lambertw", "z", QT_TRANSLATE_NOOP("FormulaHelp", "principal branch of the Lambert W of z"), NULL}, /* 1 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "choosing by iteration")},
    {"ifiter", "first, second, ...", QT_TRANSLATE_NOOP("FormulaHelp", "picks one argument by iteration number, cycling: the first on the first iteration, the second on the next, and round again; only that one runs"), NULL}, /* variadic */
    {"ifiterl", "first, second, ...", QT_TRANSLATE_NOOP("FormulaHelp", "as ifiter, but holding the last argument once the iterations run out"), NULL}, /* variadic */
    {"ifiterf", "rest, last", QT_TRANSLATE_NOOP("FormulaHelp", "last on the last iteration the limit allows, rest on the others; only the chosen one runs"), NULL}, /* 2 */
    {"ifiterr", "before, after, n", QT_TRANSLATE_NOOP("FormulaHelp", "before while the iteration is below n, after from n on; only the chosen one runs"), NULL}, /* 3 */
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "randomness: all but the seed optional, and a new field every iteration; what the kaleidoscope and kaled_mirror do is in the Values tab")},
    {"rand", "scale", QT_TRANSLATE_NOOP("FormulaHelp", "real(scale) times a random number in [0, 1); depends on call order, so a redraw differs"), NULL},
    {"randsc", "seed, [size=1+i], [degradation=0.5+0.5i], [kaleidoscope=1], [kaled_mirror=0], [skew=0], [skew_mode=0], [selfsim=off]", QT_TRANSLATE_NOOP("FormulaHelp", "coherent noise over the point: soft blobs, size wide and size high -- gradient noise, one octave of what fbm sums"), NULL},
    {"randscq", "seed, [size=1+i], [degradation=0.5+0.5i], [kaleidoscope=1], [kaled_mirror=0], [skew=0], [skew_mode=0], [selfsim=off]", QT_TRANSLATE_NOOP("FormulaHelp", "the same field with no interpolation: a mosaic of flat square cells"), NULL},
    {"randscp", "seed, [size=1+i], [degradation=0.5+0.5i], [kaleidoscope=1], [kaled_mirror=0], [skew=0], [skew_mode=0], [selfsim=off]", QT_TRANSLATE_NOOP("FormulaHelp", "the same field cut into irregular flat polygons, with straight edges"), NULL},
    {"randsch", "seed, [size=1+i], [degradation=0.5+0.5i], [kaleidoscope=1], [kaled_mirror=0], [skew=0], [skew_mode=0], [selfsim=off]", QT_TRANSLATE_NOOP("FormulaHelp", "the same field cut into hexagons: a honeycomb of flat cells"), NULL},
    {"randsct", "seed, [size=1+i], [degradation=0.5+0.5i], [kaleidoscope=1], [kaled_mirror=0], [skew=0], [skew_mode=0], [selfsim=off]", QT_TRANSLATE_NOOP("FormulaHelp", "the same field cut into equilateral triangles, alternating in orientation"), NULL},
    {"randsctile", "tiling, seed, [size=1+i], [degradation=0.5+0.5i], [kaleidoscope=1], [kaled_mirror=0], [skew=0], [skew_mode=0], [selfsim=off]", QT_TRANSLATE_NOOP("FormulaHelp", "the same field over one of forty-five tilings, chosen by the first argument (see the values, and the Tilings tab for a picture of each): regular, Archimedean and their duals, bricks and other patterns, Islamic stars, Voronoi, Penrose, Ammann-Beenker and more that never repeat. The rest as randsc"), NULL},
    {"fbm", "value, seed, [intensity=4], [frequency=8], [octaves=4], [roughness=0.5]", QT_TRANSLATE_NOOP("FormulaHelp", "a fractional Brownian motion over the value written in front of it, the same motion the fbm colouring modes draw: octaves of gradient noise, each at twice the frequency of the one before and keeping roughness of its height. Runs from nought to intensity and never below. fbm(z,7) moves with the orbit, fbm(x,7) stands still on the plane; parchmenta(z,6) inside it folds the plane into sectors"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "watching the orbit: both hand back their argument until the last iteration, and what they gathered on it, which the inside colouring modes then draw")},
    {"trap", "orbit, [shape=0], [centre=0], [size=1]", QT_TRANSLATE_NOOP("FormulaHelp", "how near the orbit ever came to a shape; the shapes are in the Values tab"), NULL},
    {"stripe", "orbit, [density=4]", QT_TRANSLATE_NOOP("FormulaHelp", "the average of (sin(density*arg orbit)+1)/2 along the orbit; density a whole number, how many stripes go round a turn"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "figures, drawn the same in mandelbrot and julia mode as the noise is -- but fractals of their own rather than fields: each carries the point to the parent of the part it stands in, so the pass it leaves on is the level of that part and writing one alone draws the figure. The topmost part, having no parent, is thrown out of the bailout rather than declared gone, so that the outside colouring modes that read z have a point of the plane to read there. radius is read the way bailout is, as the square of the distance, and the figure is inscribed in the shape a bailout of that number draws: the gasket in the triangle at -90 degrees, the carpet in the square, the snowflake with its points on the corners of the hexagon at 0. Each takes the kaleidoscope the randsc family takes, written last: a figure folds z, the point it is carrying, where a field folds the position, and what comes of it is the same -- the figure drawn in one wedge and repeated round the origin. A bailout larger than the radius leaves a margin between figure and shape: that is no part of the figure and gets no iteration, being turned by the figure's own symmetry so that it never leaves and takes the inside colour")},
    {"sierpinskyt", "[radius=4], [kaleidoscope=1], [kaled_mirror=0]", QT_TRANSLATE_NOOP("FormulaHelp", "the Sierpinski gasket: doubling away from the nearest corner carries a point to the parent of its hole, so it leaves on the pass numbered by the cut that took it"), NULL},
    {"sierpinskyc", "[radius=4], [squares=3], [kaleidoscope=1], [kaled_mirror=0]", QT_TRANSLATE_NOOP("FormulaHelp", "the Sierpinski carpet: cut into squares by squares, the border ring kept and all it encloses thrown away, so one square stands in the middle and 4*squares-4 around it. 3 is the usual carpet, 5 or 7 a lacier one, 2 a gasket"), NULL},
    {"snowflake", "[radius=4], [kaleidoscope=1], [kaled_mirror=0]", QT_TRANSLATE_NOOP("FormulaHelp", "the Koch snowflake, read from its middle out and banded by generation: a regular hexagon, six triangles of its own side standing on its six sides, twelve of a third that on their free edges, and so on. Each step onto a parent expands -- by the root of three onto the hexagon, by three below that -- so that the figure stays sensitive when it is written inside a larger formula. Whether a point is ground is decided by walking down the Koch curve one level a pass, so what a picture shows is what its iteration count paid for. The ground it leaves in the corners of the hexagon is no part of the figure: it is turned half about, which lands ground on ground, so it never leaves and takes the inside colour -- turned rather than left where it stood, so that the point keeps moving as it does in the other two"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "polynomials")},
    {"poly", "z, k1, k2, ...", QT_TRANSLATE_NOOP("FormulaHelp", "k1*z^(m-1) + k2*z^(m-2) + ... + km: the first coefficient written multiplies the highest power, the last stands alone"), NULL},
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
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "trap: which shape the orbit is measured against")},
    {"0", "", QT_TRANSLATE_NOOP("FormulaHelp", "the centre itself, a point"), NULL},
    {"1", "", QT_TRANSLATE_NOOP("FormulaHelp", "a horizontal line through the centre"), NULL},
    {"2", "", QT_TRANSLATE_NOOP("FormulaHelp", "a vertical line through the centre"), NULL},
    {"3", "", QT_TRANSLATE_NOOP("FormulaHelp", "both of them, a cross"), NULL},
    {"4", "", QT_TRANSLATE_NOOP("FormulaHelp", "a ring of radius size"), NULL},
    {"5", "", QT_TRANSLATE_NOOP("FormulaHelp", "a square of half-side size"), NULL},
    {"6", "", QT_TRANSLATE_NOOP("FormulaHelp", "a diamond of half-diagonal size"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "randsctile: the tiling, its first argument. Every one scaled to a tile of unit area on average, so it changes the shape of the cells and not their scale. The ones that never repeat are worked out point by point and cost the most; they go flat past four thousand million cells from the origin. The Tilings tab has a picture of each")},
    {"1", "", QT_TRANSLATE_NOOP("FormulaHelp", "squares"), NULL},
    {"2", "", QT_TRANSLATE_NOOP("FormulaHelp", "equilateral triangles"), NULL},
    {"3", "", QT_TRANSLATE_NOOP("FormulaHelp", "hexagons"), NULL},
    {"4", "", QT_TRANSLATE_NOOP("FormulaHelp", "octagons and squares (4.8.8)"), NULL},
    {"5", "", QT_TRANSLATE_NOOP("FormulaHelp", "hexagons and triangles, the trihexagonal (3.6.3.6)"), NULL},
    {"6", "", QT_TRANSLATE_NOOP("FormulaHelp", "hexagons, squares and triangles (3.4.6.4)"), NULL},
    {"7", "", QT_TRANSLATE_NOOP("FormulaHelp", "dodecagons and triangles (3.12.12)"), NULL},
    {"8", "", QT_TRANSLATE_NOOP("FormulaHelp", "dodecagons, hexagons and squares (4.6.12)"), NULL},
    {"9", "", QT_TRANSLATE_NOOP("FormulaHelp", "rows of squares and of triangles (3.3.3.4.4)"), NULL},
    {"10", "", QT_TRANSLATE_NOOP("FormulaHelp", "the snub square (3.3.4.3.4)"), NULL},
    {"11", "", QT_TRANSLATE_NOOP("FormulaHelp", "the snub hexagonal, hexagons in a sea of triangles (3.3.3.3.6)"), NULL},
    {"12", "", QT_TRANSLATE_NOOP("FormulaHelp", "tetrakis square: right triangles, eight round a point"), NULL},
    {"13", "", QT_TRANSLATE_NOOP("FormulaHelp", "rhombille: rhombi that read as cubes"), NULL},
    {"14", "", QT_TRANSLATE_NOOP("FormulaHelp", "deltoidal trihexagonal: kites"), NULL},
    {"15", "", QT_TRANSLATE_NOOP("FormulaHelp", "triakis triangular: obtuse triangles"), NULL},
    {"16", "", QT_TRANSLATE_NOOP("FormulaHelp", "kisrhombille: half-equilateral triangles, twelve round a point"), NULL},
    {"17", "", QT_TRANSLATE_NOOP("FormulaHelp", "prismatic pentagons"), NULL},
    {"18", "", QT_TRANSLATE_NOOP("FormulaHelp", "Cairo pentagons"), NULL},
    {"19", "", QT_TRANSLATE_NOOP("FormulaHelp", "floret pentagons, six to a flower"), NULL},
    {"20", "", QT_TRANSLATE_NOOP("FormulaHelp", "bricks, running bond"), NULL},
    {"21", "", QT_TRANSLATE_NOOP("FormulaHelp", "Flemish bond: a brick and a half brick by turns"), NULL},
    {"22", "", QT_TRANSLATE_NOOP("FormulaHelp", "herringbone"), NULL},
    {"23", "", QT_TRANSLATE_NOOP("FormulaHelp", "basketweave"), NULL},
    {"24", "", QT_TRANSLATE_NOOP("FormulaHelp", "Pythagorean: squares of two sizes"), NULL},
    {"25", "", QT_TRANSLATE_NOOP("FormulaHelp", "chevrons"), NULL},
    {"26", "", QT_TRANSLATE_NOOP("FormulaHelp", "squares and rhombi in alternate rows"), NULL},
    {"27", "", QT_TRANSLATE_NOOP("FormulaHelp", "houses: a pentagon, up and hanging by turns"), NULL},
    {"28", "", QT_TRANSLATE_NOOP("FormulaHelp", "rows: two of squares, one of triangles"), NULL},
    {"29", "", QT_TRANSLATE_NOOP("FormulaHelp", "rows: one of squares, two of triangles"), NULL},
    {"30", "", QT_TRANSLATE_NOOP("FormulaHelp", "hexagons among triangles"), NULL},
    {"31", "", QT_TRANSLATE_NOOP("FormulaHelp", "Greek crosses"), NULL},
    {"32", "", QT_TRANSLATE_NOOP("FormulaHelp", "T tetrominoes"), NULL},
    {"33", "", QT_TRANSLATE_NOOP("FormulaHelp", "eight-pointed stars and crosses"), NULL},
    {"34", "", QT_TRANSLATE_NOOP("FormulaHelp", "six-pointed stars and hexagons"), NULL},
    {"35", "", QT_TRANSLATE_NOOP("FormulaHelp", "eight-pointed stars with four-pointed ones and rhombi"), NULL},
    {"36", "", QT_TRANSLATE_NOOP("FormulaHelp", "twelve-pointed stars with three-pointed ones and rhombi"), NULL},
    {"37", "", QT_TRANSLATE_NOOP("FormulaHelp", "twelve-pointed stars with notched hexagons and squares"), NULL},
    {"38", "", QT_TRANSLATE_NOOP("FormulaHelp", "Voronoi cells, irregular -- randscp's"), NULL},
    {"39", "", QT_TRANSLATE_NOOP("FormulaHelp", "Penrose rhombs: never repeats, fivefold"), NULL},
    {"40", "", QT_TRANSLATE_NOOP("FormulaHelp", "Penrose kites and darts"), NULL},
    {"41", "", QT_TRANSLATE_NOOP("FormulaHelp", "Ammann-Beenker: squares and rhombi, eightfold, never repeats"), NULL},
    {"42", "", QT_TRANSLATE_NOOP("FormulaHelp", "rhombi in twelve directions, never repeats"), NULL},
    {"43", "", QT_TRANSLATE_NOOP("FormulaHelp", "rhombi in seven directions, never repeats"), NULL},
    {"44", "", QT_TRANSLATE_NOOP("FormulaHelp", "the pinwheel: triangles facing every way there is"), NULL},
    {"45", "", QT_TRANSLATE_NOOP("FormulaHelp", "the chair: Ls cut into Ls"), NULL},
    {"anything else", "", QT_TRANSLATE_NOOP("FormulaHelp", "nought, as a zero size gives"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "randsc family and the three figures: how many wedges the kaleidoscope folds the plane into. The fields fold the position; a figure has no use for it and reads z, the point it carries down towards the part it stands in, so that is what a figure folds")},
    {"1", "", QT_TRANSLATE_NOOP("FormulaHelp", "no folding at all, which is what a call that says nothing gets"), NULL},
    {"2 or more", "", QT_TRANSLATE_NOOP("FormulaHelp", "that many wedges around the origin, the field taken from one"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "randsc family and the three figures: kaled_mirror, which mirror the kaleidoscope folds with")},
    {"0", "", QT_TRANSLATE_NOOP("FormulaHelp", "the far half of each wedge mirrors the near half, so a wedge is symmetric about its bisector"), NULL},
    {"1", "", QT_TRANSLATE_NOOP("FormulaHelp", "the same the other way about, the near half mirroring the far one"), NULL},
    {"anything else", "", QT_TRANSLATE_NOOP("FormulaHelp", "folds the way 0 does"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "randsc family: the skew, which turns a cell's value by how far out of the middle of its cell the point stands, so that the colouring modes have something to read inside a cell. It is measured in the cell's own geometry, so the colour follows the shape the field is cut into: squares in randscq, hexagons in randsch, triangles in randsct, each cell's own polygon in randscp, and its own blobs in randsc. The real part is that gradient, the imaginary part a flat turn the whole cell shares. A turn and not a scaling, so the modulus is untouched and a circular bailout cannot see it at all, at any skew; a polygon reads the components instead and does see it, the more so the further its corners stand out past its sides. zmag reads the modulus and so stays flat; colour with real, imag, angle or real/imag")},
    {"0", "", QT_TRANSLATE_NOOP("FormulaHelp", "the default: no turn at all, the value stays on the real axis and is the number it always was, to the bit"), NULL},
    {"0.05", "", QT_TRANSLATE_NOOP("FormulaHelp", "a seventh of a colour band across a cell: shows at a colour speed of about 8, and costs nothing under a circular, square or hexagonal bailout"), NULL},
    {"0.3 to 0.6", "", QT_TRANSLATE_NOOP("FormulaHelp", "about one band across a cell, which is what shows at a colour speed of 1. Free under a circular bailout; one to three per cent of the picture changes which cells leave under a polygon"), NULL},
    {"larger", "", QT_TRANSLATE_NOOP("FormulaHelp", "no more colour to speak of past 1, and more of the picture moving under a polygonal bailout. Pair a large skew with a circular bailout"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "randsc family: what the skew does. The numbers add together rather than choosing one another out, so 3 is the shape and the rosette at once and 11 is those two with the modulus moving as well. A skew of nought leaves every one of them doing nothing")},
    {"0 or 1", "", QT_TRANSLATE_NOOP("FormulaHelp", "the shape: the cell's own outline, shrunk step by step -- squares in randscq, hexagons in randsch, triangles in randsct, its polygon in randscp, its blobs in randsc. What a call that gives no skew_mode gets, and what is used when neither 1 nor 2 is named"), NULL},
    {"2", "", QT_TRANSLATE_NOOP("FormulaHelp", "the rosette: the angle round the middle of the cell, folded into as many turns as the kaleidoscope has wedges, so each cell carries the picture's own symmetry. Plain spokes when nothing is folded, and a spiral when added to 1"), NULL},
    {"4", "", QT_TRANSLATE_NOOP("FormulaHelp", "across the wedge: a turn that follows where the point stands in its wedge of the kaleidoscope -- nothing on the edges of the wedge, the whole turn the imaginary part of the skew gives on its own down the middle (127 degrees for 2) -- the same in every wedge, so the picture stays a kaleidoscope, mirrors and all. Does nothing when nothing is folded, or when the imaginary part is nought"), NULL},
    {"8", "", QT_TRANSLATE_NOOP("FormulaHelp", "radial: the modulus moves as well as the angle. The only one zmag and the bailout can see, since those read the modulus and nothing else -- and for that reason the only one that cuts the figure where the bailout falls inside a cell"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "randsc family: selfsim, every pass so far rather than this one alone. Each pass is a field of its own, its cells the degradation times those of the pass before, so a formula calling one on every pass reads a finer field each time -- and once the cells are smaller than a pixel, snow. With selfsim the call hands back the passes so far averaged, pass n weighed by d^(nH): d is the degradation, the geometric mean of its two components taken without their signs, and H is this number. The fine passes weigh little, so their snow lies under the coarse ones rather than over them. It is the skewed value that is averaged, so the skew and its modes go on doing what they do. An average draws the values together, so a higher colour speed brings back the contrast. With a degradation of one every pass weighs the same, and the answer is their plain average, which settles toward a flat middle")},
    {"left out", "", QT_TRANSLATE_NOOP("FormulaHelp", "off: the one pass, to the bit -- what a call that says nothing gets, and a place left empty"), NULL},
    {"0", "", QT_TRANSLATE_NOOP("FormulaHelp", "the plain average: every pass weighed alike, as a degradation of one weighs them, and what every H near nought gives. Not off"), NULL},
    {"1", "", QT_TRANSLATE_NOOP("FormulaHelp", "the plain motion: each pass weighs d times the one before"), NULL},
    {"0.5", "", QT_TRANSLATE_NOOP("FormulaHelp", "rougher: the fine passes keep more of their weight, so more of the detail shows"), NULL},
    {"larger", "", QT_TRANSLATE_NOOP("FormulaHelp", "smoother: the first few passes take nearly all of it"), NULL},
    {"negative", "", QT_TRANSLATE_NOOP("FormulaHelp", "the fine passes outweigh the coarse ones, which is the snow back again"), NULL},
    {"complex, as {1,2}", "", QT_TRANSLATE_NOOP("FormulaHelp", "the real part weighs the passes as above; the imaginary part turns each pass by that times the logarithm of d further than the one before, before they are averaged -- a spiral of octaves, and a value with an imaginary part even without a skew. The weights are divided by the sum of their sizes, so the answer stays within the values averaged. Turns nothing at a degradation of one. Off only when left out"), NULL},
    {NULL, NULL, NULL, NULL}};

/* The variables the engine registers before parsing a user formula; see the
 * sffe_regvar calls in formulas.cpp. Anything else named in a formula is a
 * parse error, which is what keeps the two lists worth writing down. */
const struct formula_help_row formula_help_variables[] = {
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "the iteration")},
    {"z", "", QT_TRANSLATE_NOOP("FormulaHelp", "the current value, and what the formula computes the next one from; in the initialization, where z starts"), NULL},
    {"c", "", QT_TRANSLATE_NOOP("FormulaHelp", "the point being iterated -- the pixel in mandelbrot mode, the constant in julia mode"), NULL},
    {"n", "", QT_TRANSLATE_NOOP("FormulaHelp", "the iteration number"), NULL},
    {"x", "", QT_TRANSLATE_NOOP("FormulaHelp", "the plain coordinate: like c, but the pixel in julia mode too, and unchanged from pass to pass"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "earlier values of z; Fractal -> Set p values on first iteration says what they are before there are any")},
    {"p1", "", QT_TRANSLATE_NOOP("FormulaHelp", "the value z had on the pass before this one"), NULL},
    {"p2", "", QT_TRANSLATE_NOOP("FormulaHelp", "the one before that, and so on"), NULL},
    {"p9999", "", QT_TRANSLATE_NOOP("FormulaHelp", "as far back as one may look; any p up to this is a variable, whether or not the picture runs that many passes"), NULL},
    {"p", "", QT_TRANSLATE_NOOP("FormulaHelp", "another name for p1"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "suffixes on z, c, x and the p's, a shorthand for the calls made most often on one. Read from left to right, each wrapping what the ones before made: c_b_p2 is parchment(bship(c), 2), p12_p2_p3 is parchment(parchment(p12, 2), 3). M is a whole number of one or more, written in figures. Not on n, which is a count")},
    {"_b", "", QT_TRANSLATE_NOOP("FormulaHelp", "bship(...)"), NULL},
    {"_bi", "", QT_TRANSLATE_NOOP("FormulaHelp", "bshipi(...)"), NULL},
    {"_br", "", QT_TRANSLATE_NOOP("FormulaHelp", "bshipr(...)"), NULL},
    {"_pM", "", QT_TRANSLATE_NOOP("FormulaHelp", "parchment(..., M)"), NULL},
    {"_paM", "", QT_TRANSLATE_NOOP("FormulaHelp", "parchmenta(..., M)"), NULL},
    {"_rNtD", "", QT_TRANSLATE_NOOP("FormulaHelp", "rotate(..., N/D*pi): turned clockwise by N/D of half a turn, p4_r3t7 being rotate(p4,3/7*pi). D is one or more, N may be 0"), NULL},
    {"_rN", "", QT_TRANSLATE_NOOP("FormulaHelp", "rotate(..., N/180*pi): turned clockwise by N degrees, z_r45 being rotate(z,45/180*pi)"), NULL},
    {NULL, NULL, NULL, NULL}};

/* How values are written. These are properties of the parser rather than
 * named things it holds, so they are listed apart from the variables. */
const struct formula_help_row formula_help_notation[] = {
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "writing a value")},
    {"{re,im}", "", QT_TRANSLATE_NOOP("FormulaHelp", "a complex number given by its two parts, as in {0,2} for 2i"), NULL},
    {"i", "", QT_TRANSLATE_NOOP("FormulaHelp", "the imaginary unit, as the shipped formulas use it"), NULL},
    {"1.5", "", QT_TRANSLATE_NOOP("FormulaHelp", "a real number, read at the full precision of this build"), NULL},
    {NULL, NULL, NULL, QT_TRANSLATE_NOOP("FormulaHelp", "calling")},
    {"f(a, b)", "", QT_TRANSLATE_NOOP("FormulaHelp", "arguments are separated by a comma, and by nothing else: a semicolon used to serve as well and is now refused where it stands"), NULL},
    {"[b=1]", "", QT_TRANSLATE_NOOP("FormulaHelp", "an argument shown in square brackets may be left out; the value shown is what the function uses instead"), NULL},
    {"f(a, ,c)", "", QT_TRANSLATE_NOOP("FormulaHelp", "one may be left out in the middle as well, by leaving its place empty; spaces make no difference"), NULL},
    {"poly(z,1, ,1)", "", QT_TRANSLATE_NOOP("FormulaHelp", "an empty coefficient of poly is a term that is not there: this is z^2 + 1"), NULL},
    {"ifiter(a, ,b)", "", QT_TRANSLATE_NOOP("FormulaHelp", "an empty branch of ifiter or ifiterl repeats the one before it: two iterations of a, then b"), NULL},
    {"2z", "", QT_TRANSLATE_NOOP("FormulaHelp", "multiplication may be left out before a name or a bracket"), NULL},
    {"(a)(b)", "", QT_TRANSLATE_NOOP("FormulaHelp", "the same, between two bracketed groups"), NULL},
    {"-z^2", "", QT_TRANSLATE_NOOP("FormulaHelp", "a leading minus binds looser than a power: it means -(z^2)"), NULL},
    {NULL, NULL, NULL, NULL}};

/* The numbered rows under the randsctile heading of the Values tab, which the
 * Tilings tab captions its pictures with. Found here rather than written out a
 * second time, so that the two tabs cannot come to say different things; a
 * test checks that what is found is every tiling the parser draws. */
static int whole_number(const char *text)
{
    if (!*text)
        return -1;
    int n = 0;
    for (; *text; text++) {
        if (*text < '0' || *text > '9')
            return -1;
        n = n * 10 + (*text - '0');
    }
    return n;
}

int formula_help_tilings(const struct formula_help_row **rows)
{
    static const char heading[] = "randsctile:";
    const struct formula_help_row *r = formula_help_values;
    for (; r->name || r->section; r++)
        if (r->section && !strncmp(r->section, heading, sizeof heading - 1))
            break;
    *rows = NULL;
    if (!r->section)
        return 0;
    *rows = ++r;
    int count = 0;
    while (r[count].name && whole_number(r[count].name) == count + 1)
        count++;
    return count;
}
