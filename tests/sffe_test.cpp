/*
 * Regression tests for the SFFE formula parser used by XaoS user formulas.
 *
 * Each case runs in its own process (see tests/CMakeLists.txt) so that a
 * segfault in the parser is reported as a single failing case instead of
 * taking the whole suite down with it.
 *
 *   sffe_test                 run every case in-process (stops at a crash)
 *   sffe_test --case N        run case N only
 *   sffe_test --list          print all cases
 *   sffe_test --expect-count N  fail unless there are exactly N cases
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include "sffe.h"

/* XaoS provides this out of src/ui; the parser only uses it for error text. */
const char *qt_gettext(const char * /*context*/, const char *text)
{
    return text;
}

enum expectation {
    Ok,    /* must parse and evaluate to a finite value */
    OkVal, /* ...and match the expected value */
    Err    /* must report this error without crashing */
};

struct testcase {
    const char *expr;
    expectation expect;
    int error;     /* enum sffe_error, when expect == Err */
    double re;     /* expected value, when expect == OkVal */
    double im;
    unsigned iter; /* value of sffe_iteration during evaluation */
    /* How many times c0/c1/c2 must each have run after LAZY_ITERS iterations,
     * as "a,b,c". NULL for cases that do not measure evaluation. */
    const char *counts;
    const char *note;
};

#define T_OK(e, n) {e, Ok, 0, 0, 0, 0, NULL, n}
#define T_VAL(e, r, i, n) {e, OkVal, 0, r, i, 0, NULL, n}
#define T_ERR(e, c, n) {e, Err, c, 0, 0, 0, NULL, n}
/* ...evaluated on iteration k, for the iteration-dependent functions */
#define T_VAL_AT(e, k, r, i, n) {e, OkVal, 0, r, i, k, NULL, n}
/* ...counting which branches actually ran */
#define T_LAZY(e, c, n) {e, Ok, 0, 0, 0, 0, c, n}

#define LAZY_ITERS 6

/* Counting stand-ins used as ifiter branches: an argument that never runs
 * leaves its counter alone, which is the only way to observe laziness from
 * outside (the value the formula produces is the same either way). */
static int hits[3];
#define COUNTER(n)                                                             \
    static sfarg *counter##n(sfarg *const p)                                   \
    {                                                                          \
        hits[n] += 1;                                                          \
        sfvalue(p) = sfvalue(sfaram1(p));                                      \
        return p;                                                              \
    }
COUNTER(0)
COUNTER(1)
COUNTER(2)

/*
 * Variables are bound to z=2, c=3, p=1, n=2 so that the expected values below
 * stay small integers wherever the maths allows it.
 */
static const testcase cases[] = {
    /* --- baseline arithmetic ------------------------------------------- */
    T_VAL("z^2+c", 7, 0, "baseline"),
    T_VAL("z^9+c", 515, 0, "baseline"),
    T_VAL("z*(-c)", -6, 0, "baseline"),

    /* --- unary minus ---------------------------------------------------
     * Some of these only pass today because sffe_parse rewrites the source
     * text in PHASE 1 (f(-x) -> f(0-x)). They must keep passing once that
     * workaround is replaced by a real unary operator.
     */
    T_VAL("-z", -2, 0, "unary minus at start of formula"),
    T_VAL("+z", 2, 0, "unary plus at start of formula"),
    T_VAL("--z", 2, 0, "stacked unary minus"),
    T_VAL("-(z+c)", -5, 0, "unary minus applied to a parenthesised group"),
    T_VAL("-z*c", -6, 0, "unary minus binds tighter than *"),
    T_VAL("-z^2", -4, 0, "unary minus binds looser than ^: -(z^2)"),
    T_VAL("-2^2", -4, 0, "same, with a literal base"),
    T_VAL("z^-2", 0.25, 0, "unary minus after ^, numeric operand"),
    T_VAL("2*-3", -6, 0, "unary minus after *, numeric operand"),
    T_VAL("z--c", 5, 0, "binary minus followed by unary minus"),
    T_VAL("z-+c", -1, 0, "binary minus followed by unary plus"),
    T_VAL("-{1,2}", -1, -2, "unary minus applied to a complex literal"),
    T_VAL("2-{1,2}", 1, -2, "binary minus before a complex literal"),
    T_VAL("{-1,-2}", -1, -2, "signs inside a complex literal"),
    T_VAL("sin(-z)", -0.909297426825681695, 0, "unary minus as sole argument"),
    T_VAL("-sin(z)", -0.909297426825681695, 0, "unary minus applied to a call"),
    T_VAL("sin(-(z))", -0.909297426825681695, 0, "unary minus before a group"),
    /* atan2s works component-wise; with real arguments the imaginary part
     * comes out of atan2(-0,-0), so use literals to pin an exact value. */
    T_OK("atan2s(-z,-c)", "unary minus in both arguments"),
    T_VAL("atan2s({1,2},{3,4})", 0.321750554396642194,
          0.463647609000806116, "component-wise atan2"),

    /* atan2(y, x) is the angle of the real parts plus i times the angle of the
     * imaginary parts, and must collapse to the plain atan2 on real input. */
    T_VAL("atan2(1,1)", 0.785398163397448310, 0, "atan2 of two reals"),
    T_VAL("atan2(z,c)", 0.588002603547567532, 0, "atan2(2,3)"),
    T_VAL("atan2(-z,-c)", -2.553590050042225713, 0,
          "negated reals stay real: no -pi from atan2(-0,-0)"),
    T_VAL("atan2(0,1)", 0, 0, "atan2 of a zero numerator"),
    T_VAL("atan2({1,2},{3,4})", 0.321750554396642194, 0.463647609000806116,
          "each component gets its own angle"),
    T_VAL("atan2({0,2},{0,4})", 0, 0.463647609000806116,
          "a zero horizontal pair contributes nothing"),

    T_ERR("rad(z)", UnknownFunction, "listed with no implementation"),
    T_ERR("deg(z)", UnknownFunction, "listed with no implementation"),
    T_ERR("sign(z)", UnknownFunction, "listed with no implementation"),
    T_OK("logn(z,-c)", "unary minus after a comma"),
    T_OK("logn(-z,c)", "unary minus as first argument"),

    /* These three are the gap the PHASE 1 workaround does not cover: it only
     * injects a 0 after '(' ',' or at the start, so a sign in front of a
     * non-numeric operand after * / ^ still reaches xstrtonum and fails.
     */
    T_VAL("z*-c", -6, 0, "unary minus after *, symbolic operand"),
    T_VAL("z^-c", 0.125, 0, "unary minus after ^, symbolic operand"),
    T_VAL("2^-z", 0.25, 0, "unary minus after ^, symbolic exponent"),

    /* --- wrong parameter counts must report, not crash ------------------ */
    T_ERR("sin(z,c)", InvalidParameters, "too many arguments"),
    T_ERR("sin(z,c,p)", InvalidParameters, "far too many arguments"),
    T_ERR("logn(z)", InvalidParameters, "too few arguments"),
    T_ERR("rtni(z,c)", InvalidParameters, "too few arguments, arity 3"),
    T_ERR("sin()", InvalidParameters, "empty argument list"),
    T_ERR("sin(,z)", InvalidParameters, "missing first argument"),
    /* rtni does not return its root: it stores it over its own first argument
     * and evaluates to -1. Saved formulas depend on that, so it is kept.
     * See the note on sfrtni. */
    T_VAL("rtni(2,12,6)", -1, 0, "arity 3 satisfied"),
    T_VAL("rtni(z,12,6)*(1-z)+c", 0.940536905640704735, 0,
          "the shipped heart.xpf formula keeps its value"),
    /* rtni2 is the same root returned the ordinary way: it leaves z alone, so
     * adding z back gives root + 2 rather than root + root. */
    T_VAL("rtni2(2,12,6)", -1.059463094359295265, 0, "rtni2 returns the root"),
    T_VAL("rtni2(z,12,6)+z", 0.940536905640704735, 0,
          "rtni2 does not overwrite its argument"),
    T_VAL("rtni(z,12,6)+z", -2.059463094359295265, 0,
          "...whereas rtni does, and both z read the root"),

    /* --- behaviour pinned by the table's documentation ------------------
     * These names do not mean what they look like. The comments in
     * sfcmplxfunc say so; these cases make sure the descriptions stay true.
     */
    T_VAL("sqr(3)", 9, 0, "sqr(a) is a squared"),
    T_VAL("sqr({1,2})", -3, 4, "...on complex arguments too: (1+2i)^2"),
    T_VAL("inv(4)", 0.25, 0, "inv(a) is 1/a, reachable now its name is fixed"),
    T_VAL("inv({0,2})", 0, -0.5, "1/(2i) = -i/2"),
    T_VAL("conj({3,4})", 4, 3, "conj swaps the components, it does not negate"),
    T_VAL("rect({1,2},{3,4})", 3, 2, "rect takes real from b, imaginary from a"),
    T_VAL("polar({8,0},{0,3})", 3, 0, "polar takes the modulus from b"),
    T_VAL("logn(2,8)", 3, 0, "logn takes the base first"),
    T_VAL("mod({3.7,-2.3})", 0.7, -0.3, "mod is the fractional part"),
    T_VAL("truncv({3.77,-2.22},10)", 3.7, -2.2, "truncv snaps to a 1/n grid"),
    T_VAL("truncv({3.77,-2.22},0)", 3.77, -2.22, "...and n = 0 changes nothing"),
    T_VAL("mid(0,10,1)", 10, 0, "an inverted mid range sends values to the far end"),
    T_VAL("mid(99,10,1)", 1, 0, "...at both ends"),
    T_VAL("mid(5,1,10)", 5, 0, "an ordinary mid range just clamps"),
    T_VAL("sawtooth({3.7,-2.3})", 0.7, 0.7, "sawtooth is x - floor(x)"),
    T_VAL("bship({2.7,-3.9})", 2.7, 3.9, "bship keeps the fractional parts"),

    /* --- gamma and Lambert W ---------------------------------------------
     * Both were wrong before: gamma was scaled by log(sqrt(2*pi)) instead of
     * sqrt(2*pi), and lambertw returned its own starting guess whenever the
     * Householder step cancelled to zero, which is what happens at z = 1.
     */
    T_VAL("gamma(1)", 1, 0, "gamma(n) is (n-1)!"),
    T_VAL("gamma(5)", 24, 0, "...so gamma(5) is 24, not 8.798"),
    T_VAL("gamma(10)", 362880, 0, "and gamma(10) is 9!"),
    T_VAL("gamma(0.5)", 1.772453850905516027, 0, "gamma(1/2) is sqrt(pi)"),
    T_VAL("gamma({1,1})", 0.498015668118356043, -0.154949828301810685,
          "gamma(1+i)"),
    T_VAL("gamma(z+1)-z*gamma(z)", 0, 0, "the recurrence gamma(z+1) = z*gamma(z)"),

    T_VAL("lambertw(0)", 0, 0, "W(0) is 0"),
    T_VAL("lambertw(1)", 0.567143290409783873, 0, "W(1) is the omega constant"),
    T_VAL("lambertw(e)", 1, 0, "W(e) is 1"),
    T_VAL("lambertw(2)", 0.852605502013725491, 0, "W(2)"),
    /* the defining identity, which pins the answer without a reference table */
    T_VAL("lambertw(z)*exp(lambertw(z))", 2, 0, "W(z)*e^W(z) is z again"),
    T_VAL("lambertw({-3,2})*exp(lambertw({-3,2}))", -3, 2,
          "...off the real axis too"),

    /* --- malformed input must report, not crash ------------------------- */
    T_ERR("z,c", InvalidOperators, "comma outside a function call"),
    T_ERR(",", InvalidOperators, "bare comma"),
    T_ERR("()", EmptyFormula, "empty parenthesised expression"),
    T_ERR("z+", InvalidOperators, "binary operator missing its right operand"),
    T_ERR("*z", InvalidOperators, "binary operator missing its left operand"),

    /* --- error paths that already work ---------------------------------- */
    T_ERR("", EmptyFormula, "empty formula"),
    T_ERR("sin(", UnbalancedBrackets, "unbalanced parentheses"),
    T_ERR("foo(z)", UnknownFunction, "unknown function"),
    T_ERR("q+1", UnknownVariable, "unknown variable"),

    /* --- variadic calls (ifiter/ifiterl) --------------------------------
     * A variadic callee cannot tell where its own arguments end, so the parser
     * passes it the count as an extra argument. Without that, anything around
     * the call used to be swallowed: "10+ifiter(1,2,3)" evaluated to 5.
     */
    T_VAL_AT("ifiter(1,2,3)", 0, 1, 0, "picks the first argument"),
    T_VAL_AT("ifiter(1,2,3)", 1, 2, 0, "picks the second argument"),
    T_VAL_AT("ifiter(1,2,3)", 2, 3, 0, "picks the third argument"),
    T_VAL_AT("ifiter(1,2,3)", 3, 1, 0, "wraps around to the first"),
    T_VAL_AT("ifiterl(1,2,3)", 2, 3, 0, "reaches the last argument"),
    T_VAL_AT("ifiterl(1,2,3)", 9, 3, 0, "holds on the last argument"),
    T_VAL_AT("ifiter(5)", 4, 5, 0, "a single argument is always picked"),
    T_VAL_AT("10+ifiter(1,2,3)", 1, 12, 0, "term before the call survives"),
    T_VAL_AT("99*ifiter(7,8)", 0, 693, 0, "factor before the call survives"),
    T_VAL_AT("ifiter(1,2,3)+10", 1, 12, 0, "term after the call survives"),
    T_VAL_AT("ifiter(1,2)+ifiter(3,4,5)", 2, 6, 0, "two calls in one formula"),
    T_VAL_AT("ifiter(1,2,3)*ifiterl(10,20)", 3, 20, 0, "both variants at once"),
    T_VAL_AT("ifiter(z^2,z^3,z^4)+c", 1, 11, 0, "expressions as arguments"),
    T_VAL_AT("-ifiter(1,2)", 0, -1, 0, "unary minus applied to the call"),
    T_VAL_AT("ifiter(-z,-c)", 1, -3, 0, "unary minus inside the arguments"),
    T_VAL_AT("ifiter(ifiter(1,2),3)", 0, 1, 0, "nested variadic calls"),
    T_ERR("ifiter()", InvalidParameters, "variadic call with no arguments"),
    T_ERR("ifiter(1,)", InvalidParameters, "variadic call with a trailing comma"),

    /* --- lazy evaluation -------------------------------------------------
     * Only the selected branch may run. Over 6 iterations an eagerly
     * evaluated ifiter(a,b,c) would run all three 6 times each; lazily it
     * runs 6 branches in total.
     */
    T_LAZY("ifiter(c0(z),c1(z),c2(z))", "2,2,2",
           "cycling runs each branch in turn, one per iteration"),
    T_LAZY("ifiterl(c0(z),c1(z),c2(z))", "1,1,4",
           "holding settles on the last branch"),
    T_LAZY("ifiter(c0(z),c1(z))+c2(z)", "3,3,6",
           "only the branches are skipped, not the rest of the formula"),
    T_LAZY("c0(z)*ifiter(c1(z),c2(z))", "6,3,3",
           "a term before the call still runs every iteration"),
    T_LAZY("ifiter(ifiter(c0(z),c1(z)),c2(z))", "3,0,3",
           "a skipped branch does not run the calls nested inside it"),

    /* --- formulas shipped in examples/ and catalogs/ --------------------
     * These are the compatibility corpus: whatever we change in the parser,
     * every formula XaoS ships must still parse and evaluate.
     */
    T_VAL("(z^5)+((0.8+0.4i)*(z^4))+z", 46.8, 6.4, "shipped: implicit i suffix"),
    T_OK("C*Z*LOG(Z)", "shipped"),
    T_OK("COSH(Z^3+C)SINH(Z^1.2)+C", "shipped: implicit multiplication"),
    T_OK("EXP(Z)+C", "shipped"),
    T_OK("I*LOGN(10;Z^6)^3+C", "shipped: ';' argument separator"),
    T_OK("LOGN(5;Z^2)+C", "shipped"),
    T_OK("POWD(SINH(POWD(Z;1.2));2.8)+C", "shipped: nested calls"),
    T_OK("POWD(SINH(POWD(Z;1.2));2.8)+C-0.2P^2", "shipped: implicit 0.2*P"),
    T_VAL("POWI(RABS(Z)+I*RABS(IM(Z));2)+C", 7, 0, "shipped"),
    T_OK("RTNI(Z;12;6)(1-Z)+C", "shipped: implicit multiplication after ')'"),
    T_VAL("z+(1/z)*(-1)^n", 2.5, 0, "shipped: (-1)^n"),
    T_VAL("(abs(re(z))+i*abs(im(z)))^2+c", 7, 0, "shipped: USER_FORMULA default"),
};

static const int case_count = (int)(sizeof(cases) / sizeof(cases[0]));

/* Returns 0 when the case behaves as expected, 1 otherwise. */
static int run_case(const testcase &tc)
{
    char errbuf[256];
    memset(errbuf, 0, sizeof(errbuf));

    sffe *parser = sffe_alloc();
    if (!parser) {
        printf("FAIL  cannot allocate parser\n");
        return 1;
    }

    sfNumber vz, vc, vp, vn;
    GSL_SET_COMPLEX(&vz, 2, 0);
    GSL_SET_COMPLEX(&vc, 3, 0);
    GSL_SET_COMPLEX(&vp, 1, 0);
    GSL_SET_COMPLEX(&vn, 2, 0);
    sffe_regvar(&parser, &vz, "z");
    sffe_regvar(&parser, &vc, "c");
    sffe_regvar(&parser, &vp, "p");
    sffe_regvar(&parser, &vn, "n");
    sffe_regfunc(&parser, "c0", 1, counter0);
    sffe_regfunc(&parser, "c1", 1, counter1);
    sffe_regfunc(&parser, "c2", 1, counter2);

    parser->errormsg = errbuf;
    int err = sffe_parse(&parser, tc.expr);
    parser->errormsg = NULL;

    int failed = 0;
    if (tc.expect == Err) {
        if (err == 0) {
            printf("FAIL  expected error %d, formula parsed successfully\n",
                   tc.error);
            failed = 1;
        } else {
            if (err != tc.error) {
                printf("FAIL  expected error %d, got %d (%s)\n", tc.error, err,
                       errbuf);
                failed = 1;
            }
            /* A reported error is only useful if it names what went wrong. */
            if (errbuf[0] == '\0') {
                printf("FAIL  error %d reported with an empty message\n", err);
                failed = 1;
            } else if (strstr(errbuf, "(null)")) {
                printf("FAIL  error message has no context: \"%s\"\n", errbuf);
                failed = 1;
            }
            /*
             * uih_sffeset() reuses the same parser to restore the previous
             * formula after a failed parse, so reporting an error is only
             * useful if the parser survives it.
             */
            if (sffe_parse(&parser, "z^2+c") != 0) {
                printf("FAIL  parser unusable after the error was reported\n");
                failed = 1;
            } else {
                sfNumber back = sffe_eval(parser);
                if (fabs((double)GSL_REAL(back) - 7.0) > 1e-9) {
                    printf("FAIL  recovery formula evaluated to %g, want 7\n",
                           (double)GSL_REAL(back));
                    failed = 1;
                }
            }
        }
    } else {
        if (err != 0) {
            printf("FAIL  expected success, got error %d (%s)\n", err, errbuf);
            failed = 1;
        } else if (tc.counts) {
            memset(hits, 0, sizeof(hits));
            for (unsigned i = 0; i < LAZY_ITERS; i++) {
                sffe_iteration = i;
                sffe_eval(parser);
            }
            char got[32];
            snprintf(got, sizeof(got), "%d,%d,%d", hits[0], hits[1], hits[2]);
            if (strcmp(got, tc.counts) != 0) {
                printf("FAIL  branches ran %s times, expected %s\n", got,
                       tc.counts);
                failed = 1;
            }
        } else {
            sffe_iteration = tc.iter;
            sfNumber r = sffe_eval(parser);
            double re = (double)GSL_REAL(r);
            double im = (double)GSL_IMAG(r);
            if (!std::isfinite(re) || !std::isfinite(im)) {
                printf("FAIL  result is not finite: %g%+gi\n", re, im);
                failed = 1;
            } else if (tc.expect == OkVal) {
                double tol = 1e-9 * (1.0 + fabs(tc.re) + fabs(tc.im));
                if (fabs(re - tc.re) > tol || fabs(im - tc.im) > tol) {
                    printf("FAIL  expected %g%+gi, got %g%+gi\n", tc.re, tc.im,
                           re, im);
                    failed = 1;
                }
            }
        }
    }

    sffe_free(&parser);
    return failed;
}

static void describe(int i)
{
    printf("case %d: \"%s\"  (%s)\n", i, cases[i].expr, cases[i].note);
    fflush(stdout);
}

int main(int argc, char *argv[])
{
    if (argc == 3 && !strcmp(argv[1], "--case")) {
        int i = atoi(argv[2]);
        if (i < 0 || i >= case_count) {
            printf("no such case: %d (have %d)\n", i, case_count);
            return 2;
        }
        describe(i);
        return run_case(cases[i]);
    }

    if (argc == 3 && !strcmp(argv[1], "--expect-count")) {
        int n = atoi(argv[2]);
        if (n != case_count) {
            printf("case count changed: build expects %d, source has %d.\n"
                   "Update SFFE_TEST_CASE_COUNT in tests/CMakeLists.txt.\n",
                   n, case_count);
            return 1;
        }
        return 0;
    }

    if (argc == 2 && !strcmp(argv[1], "--list")) {
        for (int i = 0; i < case_count; i++)
            describe(i);
        return 0;
    }

    int failures = 0;
    for (int i = 0; i < case_count; i++) {
        describe(i);
        failures += run_case(cases[i]);
    }
    printf("\n%d case(s), %d failure(s)\n", case_count, failures);
    return failures ? 1 : 0;
}
