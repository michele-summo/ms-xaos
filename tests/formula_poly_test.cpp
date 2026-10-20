/* The random poly of the formula reference's Utilities tab.
 *
 * What has to hold is that it never hands out a formula the parser refuses,
 * nor one that is not the polynomial it looks like: each call it writes is
 * parsed, and evaluated against the same polynomial written out term by term
 * -- k1*z^(m-1) + ... + km -- which is a second reading of the call that
 * shares nothing with poly's own Horner loop. Around that, the fields: read as
 * the tab reads them, and kept to -- the length, and the magnitude of every
 * coefficient, whose argument has to go all the way round, or keep within
 * the delta of the first's when one is given.
 */

#include <cmath>
#include <cstdio>
#include <cstring>
#include <set>
#include <string>
#include <vector>

#include "sffe.h"
#include "formulahelp.h"

const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

static int failures = 0;

static void check(int condition, const char *what)
{
    if (condition) {
        printf("ok     %s\n", what);
    } else {
        printf("FAIL   %s\n", what);
        failures++;
    }
}

/* The arguments of a call as the generator writes them: "poly(" and ")"
 * around arguments parted by ", ". A complex coefficient holds a comma of its
 * own but no space, so it is not cut. */
static bool arguments(const std::string &call, std::vector<std::string> &out)
{
    out.clear();
    if (call.compare(0, 5, "poly(") != 0 || call.back() != ')')
        return false;
    const std::string inside = call.substr(5, call.size() - 6);
    size_t from = 0;
    for (;;) {
        const size_t at = inside.find(", ", from);
        if (at == std::string::npos) {
            out.push_back(inside.substr(from));
            return true;
        }
        out.push_back(inside.substr(from, at - from));
        from = at + 2;
    }
}

static sfNumber vz, vc;

static sffe *parser(const char *text)
{
    sffe *p = sffe_alloc();
    sffe_regvar(&p, &vz, "z");
    sffe_regvar(&p, &vc, "c");
    if (sffe_parse(&p, text) != 0) {
        sffe_free(&p);
        return NULL;
    }
    return p;
}

/* A coefficient read by the parser on its own, as a value; false if it does
 * not parse. */
static bool value_of(const std::string &text, double *re, double *im)
{
    sffe *k = parser(text.c_str());
    if (!k)
        return false;
    const sfNumber v = sffe_eval(k);
    sffe_free(&k);
    *re = (double)GSL_REAL(v);
    *im = (double)GSL_IMAG(v);
    return true;
}

/* Everything a call is made of: no c, no locale's decimal comma, no exponent,
 * no stray character. */
static bool only_expected_characters(const std::string &call)
{
    for (char ch : call)
        if (!strchr("poly(z, ){}-.0123456789", ch))
            return false;
    return true;
}

/* What the calls examined so far were made of, for the checks on their
 * spread. */
static std::vector<double> moduli, angles;

/* One call checked through: its shape, its length, the magnitude of each
 * coefficient, that the parser takes it, and that it is the polynomial
 * written out. Returns what went wrong, or NULL. */
static const char *examine(const std::string &call, int length, int lo, int hi,
                           double smallest, double largest)
{
    std::vector<std::string> args;
    if (!arguments(call, args) || args.empty() || args[0] != "z")
        return "not poly(z, ...)";
    if (!only_expected_characters(call))
        return "a character that is not part of a number";
    const int m = (int)args.size() - 1;
    if (m != length)
        return "the length it says is not the number of coefficients";
    if (m < lo || m > hi)
        return "a length outside the two asked for";

    /* Each part keeps four significant figures of the magnitude, so the
     * magnitude written is the one drawn to within a thousandth of it -- and
     * twelve decimals are the most, a millionth of a millionth. */
    double size = 0;
    for (int j = 1; j <= m; j++) {
        double re, im;
        if (args[j].empty())
            return "an empty coefficient";
        if (!value_of(args[j], &re, &im))
            return "a coefficient does not parse alone";
        const double r = std::hypot(re, im);
        if (r < smallest * (1 - 1e-3) - 1e-12 ||
            r > largest * (1 + 1e-3) + 1e-12)
            return "a coefficient's magnitude outside the two asked for";
        moduli.push_back(r);
        if (r > 0)
            angles.push_back(std::atan2(im, re));
        size += r;
    }

    std::string expanded;
    for (int j = 1; j <= m; j++) {
        if (!expanded.empty())
            expanded += "+";
        expanded += "(" + args[j] + ")";
        if (m - j)
            expanded += "*z^" + std::to_string(m - j);
    }

    sffe *poly = parser(call.c_str());
    if (!poly)
        return "the parser refuses the call";
    sffe *sum = parser(expanded.c_str());
    if (!sum) {
        sffe_free(&poly);
        return "the parser refuses the call written out";
    }
    const char *wrong = NULL;
    /* z on the unit circle, where no power grows: the two readings then
     * differ by rounding, a little of it per term and in proportion to the
     * coefficients. */
    for (int t = 0; t < 4 && !wrong; t++) {
        const double angle = 0.7 + 1.9 * t;
        GSL_SET_COMPLEX(&vz, cos(angle), sin(angle));
        GSL_SET_COMPLEX(&vc, -1.3 + 0.9 * t, 1.1 - 0.6 * t);
        const sfNumber a = sffe_eval(poly);
        const sfNumber b = sffe_eval(sum);
        const double apart =
            std::hypot((double)(GSL_REAL(a) - GSL_REAL(b)),
                       (double)(GSL_IMAG(a) - GSL_IMAG(b)));
        if (!(apart <= 1e-12 * (1 + m) * (1 + size)))
            wrong = "the call is not the polynomial it looks like";
    }
    sffe_free(&poly);
    sffe_free(&sum);
    return wrong;
}

int main()
{
    char what[256];

    /* --- the fields ------------------------------------------------------ */
    struct {
        const char *text;
        int reading, length;
    } lengths[] = {
        {"", FORMULA_POLY_EMPTY, 0},
        {"   ", FORMULA_POLY_EMPTY, 0},
        {"abc", FORMULA_POLY_NOT_NUMBER, 0},
        {"3x", FORMULA_POLY_NOT_NUMBER, 0},
        {"3,5", FORMULA_POLY_NOT_NUMBER, 0},
        {"inf", FORMULA_POLY_NOT_NUMBER, 0},
        {"nan", FORMULA_POLY_NOT_NUMBER, 0},
        {"3.5", FORMULA_POLY_NOT_WHOLE, 0},
        {"0", FORMULA_POLY_TOO_SHORT, 0},
        {"-2", FORMULA_POLY_TOO_SHORT, 0},
        {"255", FORMULA_POLY_TOO_LONG, 0},
        {"1", FORMULA_POLY_OK, 1},
        {" 7 ", FORMULA_POLY_OK, 7},
        {"3.0", FORMULA_POLY_OK, 3},
        {"254", FORMULA_POLY_OK, 254},
    };
    for (const auto &f : lengths) {
        int length = -1;
        const int reading = formula_poly_length(f.text, &length);
        sprintf(what, "the length \"%s\" reads as %d%s", f.text, f.reading,
                f.reading == FORMULA_POLY_OK ? ", with its value" : "");
        check(reading == f.reading &&
                  (reading != FORMULA_POLY_OK || length == f.length),
              what);
    }

    /* An empty magnitude is its default, which the tab gives: 0 for the
     * minimum and 1 for the maximum. */
    struct {
        const char *text;
        double fallback;
        int reading;
        double magnitude;
    } magnitudes[] = {
        {"", 0, FORMULA_POLY_OK, 0},
        {"  ", 1, FORMULA_POLY_OK, 1},
        {"0.25", 1, FORMULA_POLY_OK, 0.25},
        {" 2 ", 0, FORMULA_POLY_OK, 2},
        {"0", 1, FORMULA_POLY_OK, 0},
        {"1e3", 0, FORMULA_POLY_OK, 1000},
        {"1000000", 0, FORMULA_POLY_OK, 1000000},
        {"abc", 0, FORMULA_POLY_NOT_NUMBER, 0},
        {"0,5", 0, FORMULA_POLY_NOT_NUMBER, 0},
        {"inf", 0, FORMULA_POLY_NOT_NUMBER, 0},
        {"-0.5", 0, FORMULA_POLY_NEGATIVE, 0},
        {"1000001", 0, FORMULA_POLY_TOO_LARGE, 0},
    };
    for (const auto &f : magnitudes) {
        double magnitude = -1;
        const int reading =
            formula_poly_magnitude(f.text, f.fallback, &magnitude);
        sprintf(what, "the magnitude \"%s\" (empty is %g) reads as %d%s",
                f.text, f.fallback, f.reading,
                f.reading == FORMULA_POLY_OK ? ", with its value" : "");
        check(reading == f.reading &&
                  (reading != FORMULA_POLY_OK || magnitude == f.magnitude),
              what);
    }

    /* An empty delta is none, which reads as a negative one. */
    struct {
        const char *text;
        int reading;
        double delta;
    } deltas[] = {
        {"", FORMULA_POLY_OK, -1},
        {"  ", FORMULA_POLY_OK, -1},
        {"0", FORMULA_POLY_OK, 0},
        {" 30 ", FORMULA_POLY_OK, 30},
        {"12.5", FORMULA_POLY_OK, 12.5},
        {"360", FORMULA_POLY_OK, 360},
        {"abc", FORMULA_POLY_NOT_NUMBER, 0},
        {"30,5", FORMULA_POLY_NOT_NUMBER, 0},
        {"-1", FORMULA_POLY_NEGATIVE, 0},
        {"360.5", FORMULA_POLY_TOO_WIDE, 0},
    };
    for (const auto &f : deltas) {
        double delta = -2;
        const int reading = formula_poly_delta(f.text, &delta);
        sprintf(what, "the delta \"%s\" reads as %d%s", f.text, f.reading,
                f.reading == FORMULA_POLY_OK ? ", with its value" : "");
        check(reading == f.reading &&
                  (reading != FORMULA_POLY_OK || delta == f.delta),
              what);
    }

    /* --- with no delta, the calls drawn before there was one --------------
     *
     * Written down from the generator as it was before the delta, so that the
     * same seed still gives the same call, to the character. */
    {
        int m = 0;
        check(formula_poly_random(3, 5, 0.5, 2, 1, &m) ==
                      "poly(z, {1.592,-0.293}, {-1.096,0.4}, {1.18,-1.145}, "
                      "{-0.284,1.253}, {-1.393,0.958})" &&
                  formula_poly_random(3, 5, 0.5, 2, 42, &m, -1) ==
                      "poly(z, {-0.1322,0.728}, {0.987,0.241}, "
                      "{0.355,1.767}, {-0.911,1.437})",
              "with no delta, seeds 1 and 42 give the calls they gave before "
              "there was a delta");
    }

    /* --- every call parses, keeps to the fields and is its polynomial -----
     *
     * The ends are there on purpose: one coefficient, the longest call the
     * parser takes, magnitudes far below one and far above it -- where the
     * figures written have to follow the size -- and one magnitude only. */
    struct {
        int lo, hi;
        double smallest, largest;
    } ranges[] = {
        {1, 1, 0, 1},      {1, 3, 0, 1},       {2, 2, 0.5, 0.5},
        {3, 8, 0.5, 2},    {2, 20, 0, 1},      {1, 254, 0.001, 0.002},
        {3, 6, 100, 1000}, {254, 254, 0, 1},   {4, 4, 0, 1000000},
        {3, 3, 1e-9, 1e-8},
    };
    for (const auto &r : ranges) {
        const char *wrong = NULL;
        std::string example;
        for (unsigned long long seed = 1; seed <= 200 && !wrong; seed++) {
            int m = -1;
            example = formula_poly_random(r.lo, r.hi, r.smallest, r.largest,
                                          seed, &m);
            wrong = examine(example, m, r.lo, r.hi, r.smallest, r.largest);
        }
        if (wrong)
            sprintf(what, "lengths %d to %d, magnitudes %g to %g: %s, in %.90s",
                    r.lo, r.hi, r.smallest, r.largest, wrong, example.c_str());
        else
            sprintf(what,
                    "lengths %d to %d, magnitudes %g to %g: 200 calls parse, "
                    "keep to both and are their polynomials",
                    r.lo, r.hi, r.smallest, r.largest);
        check(!wrong, what);
    }

    /* --- and they are random ------------------------------------------- */
    {
        std::set<std::string> seen;
        std::set<int> lengths_seen;
        int same = 1;
        moduli.clear();
        angles.clear();
        for (unsigned long long seed = 1; seed <= 600; seed++) {
            int m = 0;
            const std::string call =
                formula_poly_random(3, 8, 0.5, 2, seed, &m);
            seen.insert(call);
            lengths_seen.insert(m);
            if (formula_poly_random(3, 8, 0.5, 2, seed, NULL) != call)
                same = 0;
            examine(call, m, 3, 8, 0.5, 2);
        }
        sprintf(what,
                "600 seeds give %d different calls and every length from 3 to "
                "8 (%d of them)",
                (int)seen.size(), (int)lengths_seen.size());
        check(seen.size() == 600 && lengths_seen.size() == 6, what);
        check(same, "the same seed gives the same call");

        /* the magnitudes fill the range evenly: each quarter of it holds
         * about a quarter of them, and the ends are reached */
        int quarter[4] = {0, 0, 0, 0};
        double least = 1e9, most = 0;
        for (double r : moduli) {
            const int q = (int)((r - 0.5) / 1.5 * 4);
            quarter[q < 0 ? 0 : q > 3 ? 3 : q]++;
            least = std::fmin(least, r);
            most = std::fmax(most, r);
        }
        const int n = (int)moduli.size();
        sprintf(what,
                "%d magnitudes from 0.5 to 2 fill the range evenly (%d %d %d "
                "%d a quarter) and reach its ends (%.3f to %.3f)",
                n, quarter[0], quarter[1], quarter[2], quarter[3], least, most);
        check(quarter[0] > n / 5 && quarter[1] > n / 5 && quarter[2] > n / 5 &&
                  quarter[3] > n / 5 && least < 0.52 && most > 1.98,
              what);

        /* and the arguments go all the way round: each of eight sectors of
         * the turn holds about an eighth of them */
        int sector[8] = {0};
        for (double a : angles) {
            int s = (int)std::floor((a + M_PI) / (2 * M_PI) * 8);
            sector[s < 0 ? 0 : s > 7 ? 7 : s]++;
        }
        const int k = (int)angles.size();
        int fewest = k;
        for (int s : sector)
            fewest = s < fewest ? s : fewest;
        sprintf(what,
                "their arguments go all the way round: the emptiest of eight "
                "sectors holds %d of %d",
                fewest, k);
        check(fewest > k / 10, what);
    }

    /* --- a delta keeps the arguments near the first's ---------------------
     *
     * Every argument after the first within half the delta of it, measured
     * the short way round; and spread over the whole of that width, not
     * bunched: the farthest of all of them reaches nearly half the delta.
     * The length and the magnitudes are the ones the seed gives without a
     * delta, since the draws are the same. The magnitudes are kept from 0.5
     * up, where four significant figures hold an argument to a thousandth of
     * a radian or better. A call of one coefficient has no other argument to
     * keep near, and is the call drawn without a delta. */
    {
        const double deltas_tried[] = {0, 1, 30, 90, 180, 360};
        for (double delta : deltas_tried) {
            const char *wrong = NULL;
            double farthest = 0;
            for (unsigned long long seed = 1; seed <= 300 && !wrong; seed++) {
                int m = 0, m_free = 0;
                const std::string call =
                    formula_poly_random(2, 12, 0.5, 2, seed, &m, delta);
                const std::string free =
                    formula_poly_random(2, 12, 0.5, 2, seed, &m_free);
                wrong = examine(call, m, 2, 12, 0.5, 2);
                if (!wrong && m != m_free)
                    wrong = "the length is not the one drawn without a delta";
                std::vector<std::string> args, args_free;
                arguments(call, args);
                arguments(free, args_free);
                double first = 0;
                for (int j = 1; j <= m && !wrong; j++) {
                    double re, im, re_free, im_free;
                    value_of(args[j], &re, &im);
                    value_of(args_free[j], &re_free, &im_free);
                    const double r = std::hypot(re, im);
                    if (std::fabs(r - std::hypot(re_free, im_free)) > 2e-3 * r)
                        wrong = "a magnitude is not the one drawn without a "
                                "delta";
                    const double a = std::atan2(im, re);
                    if (j == 1) {
                        first = a;
                        if (args[1] != args_free[1])
                            wrong = "the first coefficient is not the one "
                                    "drawn without a delta";
                        continue;
                    }
                    double apart = std::fabs(a - first);
                    apart = std::fmin(apart, 2 * M_PI - apart);
                    farthest = std::fmax(farthest, apart);
                    if (apart > delta / 2 * M_PI / 180 + 2e-3)
                        wrong = "an argument farther from the first's than "
                                "half the delta";
                }
            }
            const double half = std::fmin(delta / 2, 180) * M_PI / 180;
            if (!wrong && farthest < 0.95 * half - 2e-3)
                wrong = "the arguments keep well within the delta";
            if (wrong)
                sprintf(what, "delta %g: %s", delta, wrong);
            else
                sprintf(what,
                        "delta %g: 300 calls keep every argument within %g "
                        "degrees of the first's, and reach %.2f",
                        delta, delta / 2, farthest * 180 / M_PI);
            check(!wrong, what);
        }
        int m = 0, m_free = 0;
        check(formula_poly_random(1, 1, 0.5, 2, 7, &m, 30) ==
                      formula_poly_random(1, 1, 0.5, 2, 7, &m_free) &&
                  m == 1,
              "a call of one coefficient is the same with a delta and "
              "without");
    }

    printf("%s\n", failures ? "FAILED" : "all passed");
    return failures ? 1 : 0;
}
