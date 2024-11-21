/*/////////////////////////////////////////////////////////////////////////////////////
// project : sFFe ( SegFault (or Segmentation Fault :) ) formula evalutaor )
// author  : Mateusz Malczak ( mateusz@malczak.info )
// wpage   : www.segfaultlabs.com/projects/sffe
///////////////////////////////////////////////////////////////////////////////////////
// special build for XaoS, for more info visit
// http://www.segfaultlabs.com/projects/sfXaos
/////////////////////////////////////////////////////////////////////////////////////*/

#ifdef SFFE_CMPLX_GSL

#include "sffe.h"
#include "sffe_cmplx_gsl.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>

const sffunction sfcmplxfunc[sffnctscount] = {
    /* operators */
    {sfpow, 2, "^\0"},
    {sfadd, 2, "+\0"},
    {sfsub, 2, "-\0"},
    {sfmul, 2, "*\0"},
    {sfdiv, 2, "/\0"},

    /* functions */
    {sfsin, 1, "sin\0"},
    {sfcos, 1, "cos\0"},
    {sftan, 1, "tan\0"},
    {sfcot, 1, "cot\0"},
    {sfasin, 1, "asin\0"},
    {sfacos, 1, "acos\0"},
    {sfatan, 1, "atan\0"},
    {sfacot, 1, "acot\0"},
    {sfatan2, 2, "atan2\0"},
    {sfsinh, 1, "sinh\0"},
    {sfcosh, 1, "cosh\0"},
    {sftanh, 1, "tanh\0"},
    {sfcoth, 1, "coth\0"},
    {sfexp, 1, "exp\0"},
    {sflog, 1, "log\0"},
    {sflog10, 1, "log10\0"},
    {sflog2, 1, "log2\0"},
    {sflogN, 2, "logn\0"},
    {sflogN, 2, "logcn\0"},
    /*power functions */
    {sfpow, 2, "pow\0"},
    {sfpowd, 2, "powd\0"},
    {sfpow, 2, "powi\0"},
    {sfpow, 2, "powdc\0"},
    {sfsqr, 1, "sqr\0"},
    {sfsqrt, 1, "sqrt\0"},
    {sfrtni, 3, "rtni"},
    {sfinv, 1, "inv\n"},
    {sfceil, 1, "ceil\0"},
    {sffloor, 1, "floor\0"},
    {sfabs, 1, "abs\0"},
    {sfrabs, 1, "rabs\0"},
    {sfre, 1, "re\0"},
    {sfim, 1, "im\0"},
    {sfcarg, 1, "arg\0"},
    {sfmod, 1, "mod\0"},
    {sfconj, 1, "conj\0"},

    {sfbship, 1, "bship\0"},
    {sfbshipr, 1, "bshipr\0"},
    {sfbshipi, 1, "bshipi\0"},

    {sfrect, 2, "rect\0"},
    {sfpolar, 2, "polar\0"},

    {sfmin, 2, "min\0"},
    {sfminr, 2, "minr\0"},
    {sfmini, 2, "mini\0"},
    {sfminm, 2, "minm\0"},

    {sfmax, 2, "max\0"},
    {sfmaxr, 2, "maxr\0"},
    {sfmaxi, 2, "maxi\0"},
    {sfmaxm, 2, "maxm\0"},

    {sfmid, 3, "mid\0"},
    {sfmidr, 3, "midr\0"},
    {sfmidi, 3, "midi\0"},
    {sfmidm, 3, "midm\0"},

    {sfsincos, 1, "sincos\0"},
    {sfcossin, 1, "cossin\0"},
    {sfsinr, 1, "sinr\0"},
    {sfcosr, 1, "cosr\0"},
    {sfsini, 1, "sini\0"},
    {sfcosi, 1, "cosi\0"},

    {sftancot, 1, "tancot\0"},
    {sfcottan, 1, "cottan\0"},
    {sftanr, 1, "tanr\0"},
    {sfcotr, 1, "cotr\0"},
    {sftani, 1, "tani\0"},
    {sfcoti, 1, "coti\0"},

    {sftrunc, 1, "trunc\0"},
    {sfsawtooth, 1, "sawtooth\0"},
    {sftwave, 1, "twave\0"},

    {sfjulian, 3, "julian\0"},
    {sfinveps, 2, "inveps\0"},
    {sfatan2s, 2, "atan2s\0"},

    {sfngon, 4, "ngon\0"}, //ngon(z, center, n, pow)
    {sfparchment, 2, "parchment\0"}, //z, n
    {sfparchmenta, 2, "parchmenta\0"}, //z, n

    {NULL, 1, "rad\0"},
    {NULL, 1, "deg\0"},
    {NULL, 1, "sign\0"},
    {NULL, 1, "trunc\0"},
    {sfrand, 1, "rand\0"}};

const char sfcnames[sfvarscount][6] = {"pi\0", "pi_2\0", "pi2\0",
                                       "e\0",  "i\0",    "rnd\0"};

const cfptr sfcvals[sfvarscount] = {sfcPI, sfcPI2, sfc2PI, sfcE, sfcI, sfcRND};

sfarg *sfadd(sfarg *const p)
{ /* + */
    sfvalue(p) = gsl_complex_add(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfsub(sfarg *const p)
{ /* - */
    sfvalue(p) = gsl_complex_sub(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfmul(sfarg *const p)
{ /* *  */
    sfvalue(p) = gsl_complex_mul(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfdiv(sfarg *const p)
{ /*  /   */
    sfvalue(p) = gsl_complex_div(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfsin(sfarg *const p)
{ /* sin */
    sfvalue(p) = gsl_complex_sin(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcos(sfarg *const p)
{ /* cos */
    sfvalue(p) = gsl_complex_cos(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sftan(sfarg *const p)
{ /* tan */
    sfvalue(p) = gsl_complex_tan(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcot(sfarg *const p)
{ /* ctan */
    sfvalue(p) = gsl_complex_cot(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfasin(sfarg *const p)
{ /* asin */
    sfvalue(p) = gsl_complex_arcsin(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfacos(sfarg *const p)
{ /* acos */
    sfvalue(p) = gsl_complex_arccos(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfatan(sfarg *const p)
{ /* atan */
    sfvalue(p) = gsl_complex_arctan(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfacot(sfarg *const p)
{ /* actan */
    sfvalue(p) = gsl_complex_arccot(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfatan2(sfarg *const p)
{ /* atan2 */
    return sfaram2(p);
}

sfarg *sfsinh(sfarg *const p)
{ /* sinh */
    sfvalue(p) = gsl_complex_sinh(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcosh(sfarg *const p)
{ /* cosh */
    sfvalue(p) = gsl_complex_cosh(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sftanh(sfarg *const p)
{ /* tanh */
    sfvalue(p) = gsl_complex_tanh(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcoth(sfarg *const p)
{ /* ctanh */
    sfvalue(p) = gsl_complex_coth(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfexp(sfarg *const p)
{ /* exp */
    sfvalue(p) = gsl_complex_exp(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sflog(sfarg *const p)
{ /* log */
    sfvalue(p) = gsl_complex_log(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sflog10(sfarg *const p)
{ /* log10 */
    sfvalue(p) = gsl_complex_log10(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sflog2(sfarg *const p)
{ /* log2 */
    sfNumber base;
    real(base) = 2;
    imag(base) = 0;
    sfvalue(p) = gsl_complex_log_b(sfvalue(sfaram1(p)), base);
    return sfaram1(p);
}

sfarg *sflogN(sfarg *const p)
{ /* logN */
    sfvalue(p) = gsl_complex_log_b(sfvalue(sfaram1(p)), sfvalue(sfaram2(p)));
    return sfaram2(p);
}

sfarg *sfpow(sfarg *const p)
{ /* cmplx pow */
    sfvalue(p) = gsl_complex_pow(sfvalue(sfaram2(p)), sfvalue(sfaram1(p)));
    return sfaram2(p);
}

sfarg *sfpowd(sfarg *const p)
{ /* int pow */
    sfvalue(p) = gsl_complex_pow_real(sfvalue(sfaram2(p)),
                                      GSL_REAL(sfvalue(sfaram1(p))));
    return sfaram2(p);
}

sfarg *sfsqr(sfarg *const p)
{ /* sqr */
    sfvalue(p) = gsl_complex_pow(sfvalue(sfaram1(p)), sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfsqrt(sfarg *const p)
{ /* sqrt */
    sfvalue(p) = gsl_complex_sqrt(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfrtni(sfarg *const p)
{ /* rtni */
    double nrz = pow(gsl_complex_abs(sfvalue(sfaram3(p))),
                     1.0 / (double)(int)real(sfvalue(sfaram2(p))));
    double alfi = (gsl_complex_arg(sfvalue(sfaram3(p))) +
                   8 * atan(1.0) * (double)(int)real(sfvalue(sfaram1(p)))) /
                  (double)(int)real(sfvalue(sfaram2(p)));

    cmplxset(sfvalue(sfaram3(p)), nrz * cos(alfi), nrz * sin(alfi));
    return sfaram3(p);
}

sfarg *sfinv(sfarg *const p)
{ /* cinv */
    sfvalue(p) = gsl_complex_inverse(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfceil(sfarg *const p)
{ /* ceil */
    // sfvalue(p) = ceil( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = ceil(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = ceil(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sffloor(sfarg *const p)
{ /* floor */
    // sfvalue(p) = floor( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = floor(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = floor(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcarg(sfarg *const p)
{ /* floor */
    // sfvalue(p) = floor( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = gsl_complex_arg(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = 0.0;
    return sfaram1(p);
}

sfarg *sfmod(sfarg *const p)
{ /* floor */
    // sfvalue(p) = floor( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = fmod(GSL_REAL(sfvalue(sfaram1(p))), 1);
    GSL_IMAG(sfvalue(p)) = fmod(GSL_IMAG(sfvalue(sfaram1(p))), 1);
    return sfaram1(p);
}

sfarg *sfconj(sfarg *const p)
{ /* floor */
    // sfvalue(p) = floor( sfvalue( sfaram1(p) ) );
    GSL_REAL(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfabs(sfarg *const p)
{ /* abs - |z| */
    GSL_REAL(sfvalue(p)) = gsl_complex_abs(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = 0.0;
    return sfaram1(p);
}

sfarg *sfrabs(sfarg *const p)
{ /* abs - real numbers */
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    if (GSL_REAL(sfvalue(p)) < 0)
        GSL_REAL(sfvalue(p)) = -GSL_REAL(sfvalue(p));
    GSL_IMAG(sfvalue(p)) = 0;
    return sfaram1(p);
}

sfarg *sfre(sfarg *const p)
{ /* RE */
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = 0.0;
    return sfaram1(p);
}

sfarg *sfim(sfarg *const p)
{ /* IM */
    GSL_REAL(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = 0.0;
    return sfaram1(p);
}

sfarg *sfrand(sfarg *const p)
{ /* rand */
    GSL_REAL(sfvalue(p)) =
        GSL_REAL(sfvalue(sfaram1(p))) * (double)rand() / (double)RAND_MAX;
    GSL_IMAG(sfvalue(p)) = 0;
    return sfaram1(p);
}

sfarg *sfbship(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = abs(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = abs(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfbshipr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = abs(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfbshipi(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = abs(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfrect(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    return sfaram1(p);
}

sfarg *sfpolar(sfarg *const p)
{
    double radius = gsl_complex_abs(sfvalue(sfaram1(p)));
    double theta = gsl_complex_arg(sfvalue(sfaram2(p)));
    sfvalue(p) = gsl_complex_polar(radius, theta);
    return sfaram1(p);
}

sfarg *sfmin(sfarg *const p)
{
    double r1 = GSL_REAL(sfvalue(sfaram2(p)));
    double r2 = GSL_REAL(sfvalue(sfaram1(p)));

    double i1 = GSL_IMAG(sfvalue(sfaram2(p)));
    double i2 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = r1 < r2 ? r2 : r1;
    GSL_IMAG(sfvalue(p)) = i1 < i2 ? i2 : i1;
    return sfaram2(p);
}

sfarg *sfminr(sfarg *const p)
{
    double r1 = GSL_REAL(sfvalue(sfaram2(p)));
    double r2 = GSL_REAL(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = r1 < r2 ? r2 : r1;
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    return sfaram2(p);
}

sfarg *sfmini(sfarg *const p)
{
    double i1 = GSL_IMAG(sfvalue(sfaram2(p)));
    double i2 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
    GSL_IMAG(sfvalue(p)) = i1 < i2 ? i2 : i1;
    return sfaram2(p);
}

sfarg *sfminm(sfarg *const p)
{
    double r1 = gsl_complex_abs(sfvalue(sfaram2(p)));
    double r2 = gsl_complex_abs(sfvalue(sfaram1(p)));
    double theta = gsl_complex_arg(sfvalue(sfaram2(p)));

    sfvalue(p) = gsl_complex_polar(r1 < r2 ? r2 : r1, theta);
    return sfaram2(p);
}

sfarg *sfmax(sfarg *const p)
{
    double r1 = GSL_REAL(sfvalue(sfaram2(p)));
    double r2 = GSL_REAL(sfvalue(sfaram1(p)));

    double i1 = GSL_IMAG(sfvalue(sfaram2(p)));
    double i2 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = r1 < r2 ? r1 : r2;
    GSL_IMAG(sfvalue(p)) = i1 < i2 ? i1 : i2;
    return sfaram2(p);
}

sfarg *sfmaxr(sfarg *const p)
{
    double r1 = GSL_REAL(sfvalue(sfaram2(p)));
    double r2 = GSL_REAL(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = r1 < r2 ? r1 : r2;
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram2(p)));
    return sfaram2(p);
}

sfarg *sfmaxi(sfarg *const p)
{
    double i1 = GSL_IMAG(sfvalue(sfaram2(p)));
    double i2 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram2(p)));
    GSL_IMAG(sfvalue(p)) = i1 < i2 ? i1 : i2;
    return sfaram2(p);
}

sfarg *sfmaxm(sfarg *const p)
{
    double r1 = gsl_complex_abs(sfvalue(sfaram2(p)));
    double r2 = gsl_complex_abs(sfvalue(sfaram1(p)));
    double theta = gsl_complex_arg(sfvalue(sfaram2(p)));

    sfvalue(p) = gsl_complex_polar(r1 < r2 ? r1 : r2, theta);
    return sfaram2(p);
}

double calc_mid(double v1, double v2, double v3) {
    if (v2 < v3) {
        if (v1 < v2) {
            return v2;
        } else if (v1 > v3) {
            return v3;
        } else {
            return v1;
        }
    } else {
        if (v1 < v2 && v1 < v3) {
            return v2;
        } else if (v1 > v3 && v1 > v2) {
            return v3;
        } else {
            return v1;
        }
    }
}

sfarg *sfmid(sfarg *const p)
{
    double r1 = GSL_REAL(sfvalue(sfaram3(p)));
    double r2 = GSL_REAL(sfvalue(sfaram2(p)));
    double r3 = GSL_REAL(sfvalue(sfaram1(p)));

    double i1 = GSL_IMAG(sfvalue(sfaram3(p)));
    double i2 = GSL_IMAG(sfvalue(sfaram2(p)));
    double i3 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = calc_mid(r1, r2, r3);
    GSL_IMAG(sfvalue(p)) = calc_mid(i1, i2, i3);
    return sfaram3(p);
}

sfarg *sfmidr(sfarg *const p)
{
    double r1 = GSL_REAL(sfvalue(sfaram3(p)));
    double r2 = GSL_REAL(sfvalue(sfaram2(p)));
    double r3 = GSL_REAL(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = calc_mid(r1, r2, r3);
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram3(p)));
    return sfaram3(p);
}

sfarg *sfmidi(sfarg *const p)
{
    double i1 = GSL_IMAG(sfvalue(sfaram3(p)));
    double i2 = GSL_IMAG(sfvalue(sfaram2(p)));
    double i3 = GSL_IMAG(sfvalue(sfaram1(p)));

    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram3(p)));
    GSL_IMAG(sfvalue(p)) = calc_mid(i1, i2, i3);
    return sfaram3(p);
}

sfarg *sfmidm(sfarg *const p)
{
    double r1 = gsl_complex_abs(sfvalue(sfaram3(p)));
    double r2 = gsl_complex_abs(sfvalue(sfaram2(p)));
    double r3 = gsl_complex_abs(sfvalue(sfaram1(p)));
    double theta = gsl_complex_arg(sfvalue(sfaram3(p)));

    sfvalue(p) = gsl_complex_polar(calc_mid(r1, r2, r3), theta);
    return sfaram3(p);
}

sfarg *sfsincos(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = sin(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = cos(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcossin(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = cos(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = sin(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfsinr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = sin(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcosr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = cos(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfsini(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = sin(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcosi(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = cos(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

double cot(double x) {
    return cos(x)/sin(x);
}

sfarg *sftancot(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = tan(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = cot(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcottan(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = cot(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = tan(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sftanr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = tan(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sfcotr(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = cot(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = GSL_IMAG(sfvalue(sfaram1(p)));
    return sfaram1(p);
}

sfarg *sftani(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = tan(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfcoti(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = GSL_REAL(sfvalue(sfaram1(p)));
    GSL_IMAG(sfvalue(p)) = cot(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sftrunc(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = trunc(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = trunc(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

double sawtooth(double x) {
    return x - floor(x);
}

sfarg *sfsawtooth(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = sawtooth(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = sawtooth(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

double twave(double x) {
    double xf = x/2.0;
    return 2.0*abs(2.0*(xf-floor(xf+0.5)))-1.0;
}

sfarg *sftwave(sfarg *const p)
{
    GSL_REAL(sfvalue(p)) = twave(GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = twave(GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram1(p);
}

sfarg *sfjulian(sfarg *const p)
{
    gsl_complex z = sfvalue(sfaram3(p));
    gsl_complex m;
    GSL_SET_COMPLEX(&m, gsl_complex_abs(z), 0);
    m = gsl_complex_pow(m, sfvalue(sfaram2(p)));
    gsl_complex b = sfvalue(sfaram1(p));
    double mx = GSL_REAL(m);
    double my = GSL_IMAG(m);
    double arg = gsl_complex_arg(z);
    double byg = exp(-GSL_IMAG(b)*arg);
    double bxg = arg * GSL_REAL(b);
    double cosbxg = cos(bxg);
    double sinbxg = sin(bxg);

    GSL_REAL(sfvalue(p)) = byg*(mx*cosbxg - my*sinbxg);
    GSL_IMAG(sfvalue(p)) = byg*(my*cosbxg + mx*sinbxg);
    return sfaram3(p);
}

sfarg *sfinveps(sfarg *const p)
{ /* cinv */
    double x = GSL_REAL(sfvalue(sfaram2(p)));
    double y = GSL_IMAG(sfvalue(sfaram2(p)));
    double delta = (x*x + y*y);
    GSL_REAL(sfvalue(p)) = x/(delta + GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = -y/(delta + GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram2(p);
}

sfarg *sfatan2s(sfarg *const p)
{ /* cinv */
    GSL_REAL(sfvalue(p)) = atan2(GSL_REAL(sfvalue(sfaram2(p))), GSL_REAL(sfvalue(sfaram1(p))));
    GSL_IMAG(sfvalue(p)) = atan2(GSL_IMAG(sfvalue(sfaram2(p))), GSL_IMAG(sfvalue(sfaram1(p))));
    return sfaram2(p);
}

#define M_1_2PI     0.15915494309189533576888376337251
#define M_2PI       6.283185307179586476925286766559

sfarg *sfngon(sfarg *const p)
{
    gsl_complex i;
    GSL_SET_COMPLEX(&i, 0.0, 1.0);

    gsl_complex n = sfvalue(sfaram2(p));
    gsl_complex zc = gsl_complex_sub(sfvalue(sfaram4(p)), sfvalue(sfaram3(p)));
    double t = gsl_complex_arg(zc);
    gsl_complex tn = gsl_complex_mul_real(n, t * M_1_2PI);
    tn = gsl_complex_add_real(tn, 0.5);
    GSL_REAL(tn) = floor(GSL_REAL(tn));
    GSL_IMAG(tn) = floor(GSL_IMAG(tn));
    tn = gsl_complex_mul_real(tn, M_2PI);
    tn = gsl_complex_div(tn, n);
    double cr = cos(t);
    double sr = sin(t);
    gsl_complex ccn = gsl_complex_cos(tn);
    gsl_complex scn = gsl_complex_sin(tn);
    gsl_complex rn = gsl_complex_add(gsl_complex_mul_real(ccn, cr),
                                     gsl_complex_mul_real(scn, sr));
    rn = gsl_complex_mul_real(gsl_complex_pow(rn, sfvalue(sfaram1(p))),
                              gsl_complex_abs(zc));
    gsl_complex argexp = gsl_complex_exp(gsl_complex_mul_real(i, t));
    sfvalue(p) =  gsl_complex_add(gsl_complex_mul(rn, argexp), sfvalue(sfaram3(p)));

    return sfaram4(p);
}

sfarg *sfparchment(sfarg *const p)
{
    gsl_complex z = sfvalue(sfaram2(p));
    double n = gsl_complex_abs(sfvalue(sfaram1(p)));
    //if (n == 2 && abs(GSL_REAL(z) - (-1.5)) < 0.1  && abs(GSL_IMAG(z) - (-1)) < 0.1) {
    //    int vb = 1;
    //}

    double t = gsl_complex_arg(z);
    double dN = n * M_1_2PI;
    double nN = 1/dN;

    double trc = ceil(t * dN) * nN;

    double trm = t - trc + nN;

    sfvalue(p) = gsl_complex_polar(gsl_complex_abs(z), trm);
    return sfaram2(p);
}

sfarg *sfparchmenta(sfarg *const p)
{
    gsl_complex z = sfvalue(sfaram2(p));
    double n = gsl_complex_abs(sfvalue(sfaram1(p)));
    //if (n == 5 && abs(GSL_REAL(z) - (-1.5)) < 0.1  && abs(GSL_IMAG(z) - (-1)) < 0.1) {
    //    int vb = 1;
    //}

    double t = gsl_complex_arg(z);
    double dN = n * M_1_2PI;
    double nN = 1/dN;
    double trc = ceil(t * dN) * nN;

    dN = dN*2;
    nN = 1/dN;
    double trc2 = ceil(t * dN) * nN;

    double trm = trc < trc2 + 0.1 / n ? trc2 - t : t + nN - trc2;

    sfvalue(p) = gsl_complex_polar(gsl_complex_abs(z), trm);
    return sfaram2(p);
}

// const eval
void sfcPI(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, 4 * atan(1), 0); }

void sfcPI2(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, 2 * atan(1), 0); }

void sfc2PI(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, 8 * atan(1), 0); }

void sfcE(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, exp(1), 0); }

void sfcI(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, 0, 1); }

void sfcRND(sfNumber *cnst) { GSL_SET_COMPLEX(cnst, rand(), 0); }

#endif
