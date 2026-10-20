/* The Utilities tab of Help -> User formula reference: a random polynomial,
 * written as a poly call.
 *
 * poly(z, k1, ..., km) is k1*z^(m-1) + ... + km, so m -- how many coefficients
 * the call has -- is its length, and the polynomial is of degree m-1.
 *
 * Every coefficient is a complex number drawn in polar form: its magnitude
 * evenly between the two the fields give, its argument evenly round the whole
 * turn. Given a delta, only the first coefficient's argument goes round the
 * whole turn, and every other is drawn evenly within delta degrees across,
 * centred on it: the same draws in the same order, so that the length and
 * the magnitudes are the ones the same seed gives without one. Written out, each part keeps four significant figures of the
 * magnitude, however large or small that is, so that a range of 0.001 to
 * 0.002 is not rounded to nought and one of 100 to 1000 does not trail twelve
 * decimals. A coefficient whose imaginary part rounds to nought is written as
 * the real number it is.
 *
 * There is no c in it. A user formula starts z at the pixel in mandelbrot mode
 * as well as in julia mode (calculate.h: it is not one of the formulas that
 * start at nought), so with nothing else changing from pixel to pixel the
 * picture is the polynomial's own filled Julia set, and the same in both.
 *
 * Without Qt, so that formula_poly_test can parse every formula this writes
 * and check it against the polynomial written out term by term.
 */

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>

#include "formulahelp.h"

/* The field with the spaces at either end left off, and whether it read as
 * a number: FORMULA_POLY_OK and the number in *value, or FORMULA_POLY_EMPTY
 * or FORMULA_POLY_NOT_NUMBER.
 *
 * strtod rather than a walk over the digits, so that "3.5" can be told from
 * "abc": the one is a number of the wrong kind for a length, the other none at
 * all. XaoS reads numbers in the C locale throughout, so the point is a full
 * stop here as it is in a formula. */
static int number_in(const char *text, double *value)
{
    while (*text == ' ' || *text == '\t')
        text++;
    size_t n = strlen(text);
    while (n && (text[n - 1] == ' ' || text[n - 1] == '\t'))
        n--;
    if (!n)
        return FORMULA_POLY_EMPTY;
    const std::string trimmed(text, n);
    char *end = NULL;
    *value = strtod(trimmed.c_str(), &end);
    if (end != trimmed.c_str() + n || !std::isfinite(*value))
        return FORMULA_POLY_NOT_NUMBER;
    return FORMULA_POLY_OK;
}

int formula_poly_length(const char *text, int *length)
{
    double value = 0;
    const int reading = number_in(text, &value);
    if (reading != FORMULA_POLY_OK)
        return reading;
    if (value != std::floor(value))
        return FORMULA_POLY_NOT_WHOLE;
    if (value < 1)
        return FORMULA_POLY_TOO_SHORT;
    if (value > FORMULA_POLY_LONGEST)
        return FORMULA_POLY_TOO_LONG;
    *length = (int)value;
    return FORMULA_POLY_OK;
}

int formula_poly_magnitude(const char *text, double fallback,
                           double *magnitude)
{
    double value = 0;
    const int reading = number_in(text, &value);
    if (reading == FORMULA_POLY_EMPTY) {
        *magnitude = fallback;
        return FORMULA_POLY_OK;
    }
    if (reading != FORMULA_POLY_OK)
        return reading;
    if (value < 0)
        return FORMULA_POLY_NEGATIVE;
    if (value > FORMULA_POLY_LARGEST)
        return FORMULA_POLY_TOO_LARGE;
    *magnitude = value;
    return FORMULA_POLY_OK;
}

int formula_poly_delta(const char *text, double *delta)
{
    double value = 0;
    const int reading = number_in(text, &value);
    if (reading == FORMULA_POLY_EMPTY) {
        *delta = -1;
        return FORMULA_POLY_OK;
    }
    if (reading != FORMULA_POLY_OK)
        return reading;
    if (value < 0)
        return FORMULA_POLY_NEGATIVE;
    if (value > FORMULA_POLY_WIDEST)
        return FORMULA_POLY_TOO_WIDE;
    *delta = value;
    return FORMULA_POLY_OK;
}

/* splitmix64: small, and enough for picking coefficients. */
static uint64_t next(uint64_t *state)
{
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* from lo to hi, both included */
static int between(uint64_t *state, int lo, int hi)
{
    return lo + (int)(next(state) % (uint64_t)(hi - lo + 1));
}

/* from nought up to but not including one, in steps of 2^-53 */
static double fraction(uint64_t *state)
{
    return (double)(next(state) >> 11) * (1.0 / 9007199254740992.0);
}

/* A number rounded to a count of decimals, written as the parser reads it:
 * by hand rather than through printf, whose decimal point is whatever the
 * locale says. Trailing noughts are left off, and a nought is never signed. */
static std::string decimal(double x, int decimals)
{
    long long scale = 1;
    for (int k = 0; k < decimals; k++)
        scale *= 10;
    long long v = llround(x * (double)scale);
    std::string out;
    if (v < 0) {
        out += '-';
        v = -v;
    }
    out += std::to_string(v / scale);
    long long fraction = v % scale;
    int shown = decimals;
    while (shown > 0 && fraction % 10 == 0) {
        fraction /= 10;
        shown--;
    }
    if (shown > 0) {
        std::string digits = std::to_string(fraction);
        out += '.';
        out.append(shown - digits.size(), '0');
        out += digits;
    }
    return out;
}

std::string formula_poly_random(int shortest, int longest, double smallest,
                                double largest, unsigned long long seed,
                                int *length, double delta)
{
    uint64_t state = seed;
    const int m = between(&state, shortest, longest);
    if (length)
        *length = m;

    const double turn = 6.283185307179586477;
    double first = 0;
    std::string out = "poly(z";
    for (int k = 1; k <= m; k++) {
        const double r = smallest + (largest - smallest) * fraction(&state);
        double angle = turn * fraction(&state);
        /* The same draw, spread over delta across rather than the turn. */
        if (k == 1)
            first = angle;
        else if (delta >= 0)
            angle = first + (angle / turn - 0.5) * delta * (turn / 360);
        /* Four significant figures of the magnitude: three decimals at one,
         * one more for each tenth smaller, one fewer for each ten times
         * larger. Twelve at most, which is past anything a picture shows. */
        int decimals = 0;
        if (r > 0)
            decimals = 3 - (int)std::floor(std::log10(r));
        decimals = decimals < 0 ? 0 : decimals > 12 ? 12 : decimals;
        const std::string re = decimal(r * std::cos(angle), decimals);
        const std::string im = decimal(r * std::sin(angle), decimals);
        out += ", ";
        if (im == "0")
            out += re;
        else
            out += "{" + re + "," + im + "}";
    }
    out += ')';
    return out;
}
