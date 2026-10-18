#define __USE_MINGW_ANSI_STDIO 1 // for long double support on Windows
#include <cstdio>
#include <climits>
#include <cstring>
#include <cerrno>
#include <cstdlib>
#include "filter.h"
#include "fractal.h"
#include "formulas.h"
#include "ui_helper.h"
#include "config.h"
#include "xmenu.h"
#include "play.h"
#include "i18n.h"

#ifdef USE_FLOAT128
#include <quadmath.h>
#endif

#define myputs(s)                                                              \
    ((xio_puts(s, uih->savec->file) == XIO_EOF) ? outputerror(uih), 1 : 0)
#define myputc(s)                                                              \
    ((xio_putc(s, uih->savec->file) == XIO_EOF) ? outputerror(uih), 1 : 0)
static int first;
static int changed;
static int last;
/* Keywords of a saved position, and the choices of their dialogs, which show
 * them translated in the "Dialog" context. */
const char *const save_fastmode[] = {QT_TRANSLATE_NOOP("Dialog", "zero"),
                                     QT_TRANSLATE_NOOP("Dialog", "never"),
                                     QT_TRANSLATE_NOOP("Dialog", "animation"),
                                     QT_TRANSLATE_NOOP("Dialog", "new"),
                                     QT_TRANSLATE_NOOP("Dialog", "always"),
                                     NULL};

const char *const xtextposnames[] = {QT_TRANSLATE_NOOP("Dialog", "left"),
                                     QT_TRANSLATE_NOOP("Dialog", "center"),
                                     QT_TRANSLATE_NOOP("Dialog", "right"), NULL};
const char *const ytextposnames[] = {QT_TRANSLATE_NOOP("Dialog", "top"),
                                     QT_TRANSLATE_NOOP("Dialog", "middle"),
                                     QT_TRANSLATE_NOOP("Dialog", "bottom"), NULL};

static void outputerror(struct uih_context *uih)
{
    static char error[245];
    if (uih->savec->writefailed)
        return;
    sprintf(error, "Write failed:%s", xio_errorstring());
    uih_error(uih, error);
    uih->savec->writefailed = 1;
}

static void start_save(struct uih_context *uih, const char *name)
{
    if (!changed && !uih->savec->firsttime) {
        char s[256];
        sprintf(s, "\n(usleep %i)\n", tl_lookup_timer(uih->savec->timer));
        myputs(s);
        tl_reset_timer(uih->savec->timer);
    }
    changed = 1;
    myputc('(');
    myputs(name);
    first = 0;
}

static void stop_save(struct uih_context *uih)
{
    myputc(')');
    myputc('\n');
}

#ifdef SAVEKEYWORDUSED

static void save_keyword(struct uih_context *uih, const char *name)
{
    if (!first)
        myputc(' ');
    else
        first = 0;
    myputs(name);
}
#endif

static void save_keystring(struct uih_context *uih, const char *name)
{
    if (!first)
        myputc(' ');
    else
        first = 0;
    myputc('\'');
    myputs(name);
}

static void save_float(struct uih_context *uih, number_t number)
{
    if (!first)
        myputc(' ');
    else
        first = 0;
    char s[256];
    char fs[10];
    /* Full precision for the build, so that reloading the file puts the view
     * back exactly where it was saved. */
#ifdef USE_FLOAT128
    snprintf(fs, sizeof(fs), "%%.%iQG", NUMBER_DIGITS);
    quadmath_snprintf(s, 256, fs, (__float128)number);
#else
#ifdef USE_LONG_DOUBLE
    snprintf(fs, sizeof(fs), "%%.%iLG", NUMBER_DIGITS);
    snprintf(s, 256, fs, (long double)number);
#else
    snprintf(fs, sizeof(fs), "%%.%iG", NUMBER_DIGITS);
    snprintf(s, 256, fs, (double)number);
#endif
#endif
    myputs(s);
}

static void save_float2(struct uih_context *uih, number_t number, int places)
{
    char fs[10];
    if (!first)
        myputc(' ');
    else
        first = 0;
    if (places < 0)
        places = 0;
    if (places > NUMBER_DIGITS)
        places = NUMBER_DIGITS;
    char s[256];
#ifdef USE_FLOAT128
    /* The format string is built from places and then used, which it was not:
     * this branch composed fs and went on to print a hardcoded 34 digits, so
     * asking for fewer had no effect. */
    snprintf(fs, sizeof(fs), "%%.%iQG", places);
    quadmath_snprintf(s, 256, fs, (__float128)number);
#else
#ifdef USE_LONG_DOUBLE
    snprintf(fs, sizeof(fs), "%%.%iLG", places);
    snprintf(s, 256, fs, (long double)number);
#else
    /* ...and this one passed 10 and 256 as the format strings. */
    snprintf(fs, sizeof(fs), "%%.%iG", places);
    snprintf(s, 256, fs, (double)number);
#endif
#endif
    myputs(s);
}

static void save_int(struct uih_context *uih, int number)
{
    char s[256];
    if (!first)
        myputc(' ');
    else
        first = 0;
    sprintf(s, "%i", number);
    myputs(s);
}

static void save_onoff(struct uih_context *uih, int number)
{
    if (!first)
        myputc(' ');
    else
        first = 0;
    myputs(number ? "#t" : "#f");
}

static void save_string(struct uih_context *uih, const char *text)
{
    int i = 0;
    if (!first)
        myputc(' ');
    else
        first = 0;
    myputc('"');
    while (text[i] != 0) {
        if (text[i] == '"')
            myputc('\\');
        myputc(text[i]);
        i++;
    }
    myputc('"');
}

static void save_intc(struct uih_context *uih, const char *name, int number)
{
    start_save(uih, name);
    save_int(uih, number);
    stop_save(uih);
}

static void save_onoffc(struct uih_context *uih, const char *name, int number)
{
    start_save(uih, name);
    save_onoff(uih, number);
    stop_save(uih);
}

static void save_floatc(struct uih_context *uih, const char *name,
                        number_t number)
{
    start_save(uih, name);
    save_float(uih, number);
    stop_save(uih);
}

static void save_float2c(struct uih_context *uih, const char *name,
                         number_t number, int places)
{
    start_save(uih, name);
    save_float2(uih, number, places);
    stop_save(uih);
}

static void save_coordc(struct uih_context *uih, const char *name,
                        number_t number, number_t number2)
{
    start_save(uih, name);
    save_float(uih, number);
    save_float(uih, number2);
    stop_save(uih);
}

static void save_keystringc(struct uih_context *uih, const char *name,
                            const char *param)
{
    start_save(uih, name);
    save_keystring(uih, param);
    stop_save(uih);
}

static void save_stringc(struct uih_context *uih, const char *name,
                         const char *param)
{
    start_save(uih, name);
    save_string(uih, param);
    stop_save(uih);
}

static void save_noparam(struct uih_context *uih, const char *name)
{
    start_save(uih, name);
    stop_save(uih);
}

static void save_nstring(struct uih_context *uih, int number,
                         const char *const *const texts)
{
    save_keystring(uih, texts[number]);
}

static void save_nstringc(struct uih_context *uih, const char *name, int number,
                          const char *const *const texts)
{
    save_keystringc(uih, name, texts[number]);
}

/* The view used to be written through ndecimals(), which counted how many
 * decimal places the size needed and passed that as the digit count. Two
 * things were wrong with it. The count went to a %.*G format, where it means
 * significant digits rather than places after the point, so it under-reported
 * for any centre of modulus above one. And it was capped well below what
 * number_t holds: at a zoom of 1e-13 it asked for 17 digits of a coordinate
 * carrying 21, so a deep position came back a fraction of a pixel away from
 * where it was saved, which at a high iteration count is a visibly different
 * picture.
 *
 * The view now goes out at full precision, like the juliaseed and the
 * perturbation always have. It costs a few characters per line and buys a
 * position that reloads exactly. */
static void savepos(struct uih_context *uih);

static void savepos(struct uih_context *uih)
{
    start_save(uih, "view");
    save_float(uih, uih->fcontext->s.cr);
    save_float(uih, uih->fcontext->s.ci);
    save_float(uih, uih->fcontext->s.rr);
    save_float(uih, uih->fcontext->s.ri);
    stop_save(uih);
    uih->savec->fcontext->s = uih->fcontext->s;
}

static void savepos2(struct uih_context *uih);

static void savepos2(struct uih_context *uih)
{
    start_save(uih, "animateview");
    save_float(uih, uih->fcontext->s.cr);
    save_float(uih, uih->fcontext->s.ci);
    save_float(uih, uih->fcontext->s.rr);
    save_float(uih, uih->fcontext->s.ri);
    stop_save(uih);
    uih->savec->fcontext->s = uih->fcontext->s;
}

static void savepos3(struct uih_context *uih);

static void savepos3(struct uih_context *uih)
{
    start_save(uih, "morphview");
    save_float(uih, uih->fcontext->s.cr);
    save_float(uih, uih->fcontext->s.ci);
    save_float(uih, uih->fcontext->s.rr);
    save_float(uih, uih->fcontext->s.ri);
    stop_save(uih);
    uih->savec->fcontext->s = uih->fcontext->s;
}

/* The speed change after ITER of one side, as pointers into a context, so
 * that the side being saved and the side last saved read alike. */
struct colorfun_fields {
    number_t *initial, *iter, *exp, *radius, *iter2, *exp2, *radius2;
    int *mode, *mode2;
    number_t *speed0, *speed, *speed2;
};

static colorfun_fields colorfun_side(fractal_context *c, int inset)
{
    if (inset)
        return {&c->incolorfuninitial,      &c->incolorfuniter,
                &c->incolorfunexp,          &c->incolorfunradius,
                &c->incolorfun2iter,        &c->incolorfun2exp,
                &c->incolorfun2radius,      &c->incolorfunmode,
                &c->incolorfun2mode,        &c->incolorfuninitialspeed,
                &c->incolorfunspeed,        &c->incolorfun2speed};
    return {&c->outcolorfuninitial,      &c->outcolorfuniter,
            &c->outcolorfunexp,          &c->outcolorfunradius,
            &c->outcolorfun2iter,        &c->outcolorfun2exp,
            &c->outcolorfun2radius,      &c->outcolorfunmode,
            &c->outcolorfun2mode,        &c->outcolorfuninitialspeed,
            &c->outcolorfunspeed,        &c->outcolorfun2speed};
}

static void save_colorfungroups(struct uih_context *uih, int inset,
                                colorfun_fields now, int unit);

/* The speed change after ITER of one side, inset 1 for the inside.
 *
 * Written as 1.8.0 wrote it, (incolorfunafter iter exponent 'jump) or 'off,
 * while that is all there is -- one group, off or in jump, the initial
 * exponent at one and the rest as it starts -- so that such a position opens
 * in 1.8.0 and draws there as here. Otherwise as (incolorspeedchange initial
 * iter exponent 'mode radius iter2 exponent2 'mode2 radius2), in the order of
 * the dialog, which 1.8.0 does not know. Either only when something in it has
 * moved, or for a whole position, as the fbm's numbers are; and the old one
 * puts the rest back to how it starts when it is read, so that an animation
 * that turns it off again replays as it was recorded.
 *
 * The speeds go apart, as (incolorspeedlinear initialspeed speed speed2), and
 * only where one of them is not one or they have moved: a position that does
 * not use them opens in 2.0, which knows no speeds, as it did. */
static void save_colorfunchange(struct uih_context *uih, int inset)
{
    colorfun_fields now = colorfun_side(uih->fcontext, inset);
    colorfun_fields was = colorfun_side(uih->savec->fcontext, inset);
    const int moved = *now.initial != *was.initial ||
                      *now.iter != *was.iter || *now.exp != *was.exp ||
                      *now.mode != *was.mode || *now.radius != *was.radius ||
                      *now.iter2 != *was.iter2 || *now.exp2 != *was.exp2 ||
                      *now.mode2 != *was.mode2 || *now.radius2 != *was.radius2;
    const int speedmoved = *now.speed0 != *was.speed0 ||
                           *now.speed != *was.speed ||
                           *now.speed2 != *was.speed2;
    const int unit = *now.speed0 == 1 && *now.speed == 1 && *now.speed2 == 1;
    if (!uih->save && !moved && !speedmoved)
        return;
    if (uih->save || moved)
        save_colorfungroups(uih, inset, now, unit);
    if (speedmoved || (uih->save && !unit)) {
        start_save(uih, inset ? "incolorspeedlinear" : "outcolorspeedlinear");
        save_float(uih, *now.speed0);
        save_float(uih, *now.speed);
        save_float(uih, *now.speed2);
        stop_save(uih);
    }
    *was.initial = *now.initial;
    *was.iter = *now.iter;
    *was.exp = *now.exp;
    *was.mode = *now.mode;
    *was.radius = *now.radius;
    *was.iter2 = *now.iter2;
    *was.exp2 = *now.exp2;
    *was.mode2 = *now.mode2;
    *was.radius2 = *now.radius2;
    *was.speed0 = *now.speed0;
    *was.speed = *now.speed;
    *was.speed2 = *now.speed2;
}

/* The groups themselves: as 1.8.0 wrote them while that is all there is --
 * which takes the speeds at one, since reading it puts them back there -- or
 * as (incolorspeedchange ...). */
static void save_colorfungroups(struct uih_context *uih, int inset,
                                colorfun_fields now, int unit)
{
    /* off and jump are 1.8.0's first two modes, by the same numbers */
    const int plain = unit && *now.initial == 1 && *now.mode != 2 &&
                      *now.radius == COLORFUN_RADIUS &&
                      *now.iter2 == COLORFUN2_ITER &&
                      *now.exp2 == COLORFUN2_EXP && *now.mode2 == 0 &&
                      *now.radius2 == COLORFUN2_RADIUS;
    if (plain) {
        start_save(uih, inset ? "incolorfunafter" : "outcolorfunafter");
        save_float(uih, *now.iter);
        save_float(uih, *now.exp);
        save_nstring(uih, *now.mode, colorfun_after_modes);
    } else {
        start_save(uih, inset ? "incolorspeedchange" : "outcolorspeedchange");
        save_float(uih, *now.initial);
        save_float(uih, *now.iter);
        save_float(uih, *now.exp);
        save_nstring(uih, *now.mode, colorfun_modes);
        save_float(uih, *now.radius);
        save_float(uih, *now.iter2);
        save_float(uih, *now.exp2);
        save_nstring(uih, *now.mode2, colorfun_modes);
        save_float(uih, *now.radius2);
    }
    stop_save(uih);
}

void uih_saveframe(struct uih_context *uih)
{
    struct uih_savedcontext *s = uih->savec;
    int i;
    int resetsync = 0;
    if (uih->save) {
        changed = 0;
        if (s->firsttime)
            save_noparam(uih, "initstate");
        if (s->nonfractalscreen && !uih->nonfractalscreen)
            save_noparam(uih, "display"), s->nonfractalscreen = 0;
        /* The palette goes out before the filters, not after.
         *
         * Enabling the palette filter and then setting the palette does not
         * leave the same state as doing it the other way round, and the file
         * has to describe the order the state was actually built in. Written
         * filters-first, examples/Hubicka/spiral.xpf -- default palette 34
         * with the palette filter on -- came back with different colours than
         * it was saved with, which is how this was found. */
        if (uih->palettechanged) {
            switch (uih->palettetype) {
                case 0:
                    save_intc(uih, "defaultpalette", uih->paletteshift);
                    break;
                default:
                    start_save(uih, "palette");
                    save_int(uih, uih->palettetype);
                    save_int(uih, uih->paletteseed);
                    save_int(uih, uih->paletteshift);
                    stop_save(uih);
                    break;
            }
            uih->palettechanged = 0;
            s->manualpaletteshift = 0;
        }
        if (uih->palettepickerenabled) {
            start_save(uih, "palettecolors");
            unsigned char colors[31][3];
            getDEFSEGMENTColor(colors);
            for (int i=0; i<31; i++) {
                char currcolor[6];
                rgbtohex(colors[i][0], colors[i][1],
                         colors[i][2], currcolor);
                save_string(uih, currcolor);
            }
            stop_save(uih);
        }
        if (s->manualpaletteshift != uih->manualpaletteshift)
            save_intc(uih, "shiftpalette",
                      uih->manualpaletteshift - s->manualpaletteshift),
                s->manualpaletteshift = uih->manualpaletteshift;
        for (i = uih_nfilters; i >= 0; i--) {
            if (uih->filter[i] != NULL) {
                if (s->filter[i] != 1) {
                    start_save(uih, "filter");
                    save_keystring(uih, uih->filter[i]->action->shortname);
                    save_onoff(uih, 1);
                    s->filter[i] = 1;
                    stop_save(uih);
                }
            } else if (s->filter[i] != 0) {
                s->filter[i] = 0;
                start_save(uih, "filter");
                save_keystring(uih, uih_filters[i]->shortname);
                save_onoff(uih, 0);
                stop_save(uih);
            }
        }
        if (s->fcontext->currentformula != uih->fcontext->currentformula) {
            save_keystringc(uih, "formula",
                            uih->fcontext->currentformula->shortname),
                s->fcontext->currentformula = uih->fcontext->currentformula;
#ifdef USE_SFFE
            /*SFFE : malczak */
            if (uih->fcontext->currentformula->flags & SFFE_FRACTAL) {
                if (uih->fcontext->userformula->expression)
                    save_stringc(uih, "usrform",
                                 uih->fcontext->userformula->expression);
                if (uih->fcontext->userinitial->expression)
                    save_stringc(uih, "usrformInit",
                                 uih->fcontext->userinitial->expression);
            };
/*SFFE : malczak */
#endif
            set_formula(s->fcontext,
                        (int)(uih->fcontext->currentformula - formulas));
        }
        if (s->mode >= UIH_SAVEALL)
            save_intc(uih, "letterspersec", uih->letterspersec);
        if (s->mode > UIH_SAVEPOS) {
            if (s->speedup != uih->speedup)
                save_floatc(uih, "speedup", uih->speedup),
                    s->speedup = uih->speedup;
            if (s->maxstep != uih->maxstep)
                save_floatc(uih, "maxstep", uih->maxstep),
                    s->maxstep = uih->maxstep;
            if (s->fastmode != uih->fastmode)
                save_nstringc(uih, "fastmode", uih->fastmode, save_fastmode),
                    s->fastmode = uih->fastmode;
        }
        if (s->juliamode != uih->juliamode)
            save_onoffc(uih, "fastjulia", uih->juliamode),
                s->juliamode = uih->juliamode;
        if (s->cycling != uih->cycling)
            save_onoffc(uih, "cycling", uih->cycling),
                s->cycling = uih->cycling;
        if (s->mode >= UIH_SAVEPOS &&
            s->fcontext->periodicity != uih->fcontext->periodicity)
            save_onoffc(uih, "periodicity", uih->fcontext->periodicity),
                s->fcontext->periodicity = uih->fcontext->periodicity;
        if ((uih->cycling || s->mode >= UIH_SAVEALL) &&
            (s->cyclingspeed != uih->cyclingspeed ||
             s->direction != uih->direction * uih->cyclingdirection))
            save_intc(uih, "cyclingspeed",
                      uih->cyclingspeed * uih->direction *
                          uih->cyclingdirection),
                s->cyclingspeed = uih->cyclingspeed,
                s->direction = uih->direction * uih->cyclingdirection;
        if ((s->mode > UIH_SAVEPOS && (uih->step || uih->zoomactive)) &&
            (s->xcenter != uih->xcenter || s->ycenter != uih->ycenter))
            save_coordc(uih, "zoomcenter", uih->xcenter, uih->ycenter),
                s->xcenter = uih->xcenter, s->ycenter = uih->ycenter;
        if ((!uih->fcontext->mandelbrot || uih->juliamode) &&
            (s->fcontext->pre != uih->fcontext->pre ||
             s->fcontext->pim != uih->fcontext->pim)) {
            if (uih->juliamode && uih->pressed)
                save_coordc(uih, "morphjulia", uih->fcontext->pre,
                            uih->fcontext->pim),
                    s->fcontext->pre = uih->fcontext->pre,
                    s->fcontext->pim = uih->fcontext->pim;
            else
                save_coordc(uih, "juliaseed", uih->fcontext->pre,
                            uih->fcontext->pim),
                    s->fcontext->pre = uih->fcontext->pre,
                    s->fcontext->pim = uih->fcontext->pim;
        }
        if (uih->fcontext->bre != s->fcontext->bre ||
            uih->fcontext->bim != s->fcontext->bim) {
            save_coordc(uih, "perturbation", uih->fcontext->bre,
                        uih->fcontext->bim),
                s->fcontext->bre = uih->fcontext->bre,
                s->fcontext->bim = uih->fcontext->bim;
        }
        if (uih->fastrotate != s->fastrotate && s->mode > UIH_SAVEPOS) {
            save_onoffc(uih, "fastrotate", uih->fastrotate);
            s->fastrotate = uih->fastrotate;
        }
        if (uih->fcontext->angle != s->fcontext->angle && s->autorotate != 1) {
            if (s->rotatepressed && s->mode == UIH_SAVEANIMATION)
                save_float2c(uih, "morphangle", uih->fcontext->angle, 5);
            else
                save_float2c(uih, "angle", uih->fcontext->angle, 5);
            s->rotatepressed = uih->rotatepressed;
            s->fcontext->angle = uih->fcontext->angle;
        }
        if (uih->rotationspeed != s->rotationspeed &&
            ((s->mode > UIH_SAVEPOS && uih->rotatemode == ROTATE_CONTINUOUS) ||
             s->mode >= UIH_SAVEALL)) {
            save_float2c(uih, "rotationspeed", uih->rotationspeed, 6);
            s->rotationspeed = uih->rotationspeed;
        }
        if (s->autorotate != (uih->rotatemode == ROTATE_CONTINUOUS)) {
            s->autorotate = (uih->rotatemode == ROTATE_CONTINUOUS);
            save_onoffc(uih, "autorotate", s->autorotate);
        }
        if (s->fcontext->maxiter != uih->fcontext->maxiter)
            save_intc(uih, "maxiter", uih->fcontext->maxiter),
                s->fcontext->maxiter = uih->fcontext->maxiter;
        if (s->fcontext->bailout != uih->fcontext->bailout)
            save_floatc(uih, "bailout", uih->fcontext->bailout),
                s->fcontext->bailout = uih->fcontext->bailout;
        /* Written only when it is not the circle every earlier version used,
         * so a position that does not ask for a shape stays loadable by them.
         * A file that does ask carries a command they do not know, and they
         * refuse it -- which is right, since they would draw it differently. */
        if (s->fcontext->bailoutmode != uih->fcontext->bailoutmode)
            save_intc(uih, "bailoutmode", uih->fcontext->bailoutmode),
                s->fcontext->bailoutmode = uih->fcontext->bailoutmode;
        if (uih->save || s->fcontext->coloringmode != uih->fcontext->coloringmode)
            save_intc(uih, "outcoloring", uih->fcontext->coloringmode.AsInt()),
                s->fcontext->coloringmode = uih->fcontext->coloringmode;
        if (s->fcontext->incoloringmode != uih->fcontext->incoloringmode)
            save_intc(uih, "incoloring", uih->fcontext->incoloringmode),
                s->fcontext->incoloringmode = uih->fcontext->incoloringmode;
        if (s->fcontext->incoloringmode != uih->fcontext->incoloringmode)
            save_intc(uih, "incoloring", uih->fcontext->incoloringmode),
                s->fcontext->incoloringmode = uih->fcontext->incoloringmode;
        /* The fractional Brownian motion the fbm colouring modes are drawn
         * from: five numbers a side, written as one command because they are
         * read as one, and only when one of them has moved. */
        if (uih->save ||
            s->fcontext->infbmintensity != uih->fcontext->infbmintensity ||
            s->fcontext->infbmfrequency != uih->fcontext->infbmfrequency ||
            s->fcontext->infbmroughness != uih->fcontext->infbmroughness ||
            s->fcontext->infbmoctaves != uih->fcontext->infbmoctaves ||
            s->fcontext->infbmseed != uih->fcontext->infbmseed) {
            start_save(uih, "infbmset");
            save_float(uih, uih->fcontext->infbmintensity);
            save_float(uih, uih->fcontext->infbmfrequency);
            save_int(uih, uih->fcontext->infbmoctaves);
            save_float(uih, uih->fcontext->infbmroughness);
            save_int(uih, uih->fcontext->infbmseed);
            stop_save(uih);
            s->fcontext->infbmintensity = uih->fcontext->infbmintensity;
            s->fcontext->infbmfrequency = uih->fcontext->infbmfrequency;
            s->fcontext->infbmroughness = uih->fcontext->infbmroughness;
            s->fcontext->infbmoctaves = uih->fcontext->infbmoctaves;
            s->fcontext->infbmseed = uih->fcontext->infbmseed;
        }
        if (uih->save ||
            s->fcontext->outfbmintensity != uih->fcontext->outfbmintensity ||
            s->fcontext->outfbmfrequency != uih->fcontext->outfbmfrequency ||
            s->fcontext->outfbmroughness != uih->fcontext->outfbmroughness ||
            s->fcontext->outfbmoctaves != uih->fcontext->outfbmoctaves ||
            s->fcontext->outfbmseed != uih->fcontext->outfbmseed) {
            start_save(uih, "outfbmset");
            save_float(uih, uih->fcontext->outfbmintensity);
            save_float(uih, uih->fcontext->outfbmfrequency);
            save_int(uih, uih->fcontext->outfbmoctaves);
            save_float(uih, uih->fcontext->outfbmroughness);
            save_int(uih, uih->fcontext->outfbmseed);
            stop_save(uih);
            s->fcontext->outfbmintensity = uih->fcontext->outfbmintensity;
            s->fcontext->outfbmfrequency = uih->fcontext->outfbmfrequency;
            s->fcontext->outfbmroughness = uih->fcontext->outfbmroughness;
            s->fcontext->outfbmoctaves = uih->fcontext->outfbmoctaves;
            s->fcontext->outfbmseed = uih->fcontext->outfbmseed;
        }
        /* How the outside's fbm marks are sized, and the inside's tiling and
         * how it is coloured: written only when moved, so that a position
         * that leaves them alone opens in 2.6, which knows nothing of them. */
        if (s->fcontext->outfbmscaling != uih->fcontext->outfbmscaling ||
            s->fcontext->outfbmfactor != uih->fcontext->outfbmfactor) {
            start_save(uih, "outfbmscale");
            save_nstring(uih, uih->fcontext->outfbmscaling, fbm_scalings);
            save_float(uih, uih->fcontext->outfbmfactor);
            stop_save(uih);
            s->fcontext->outfbmscaling = uih->fcontext->outfbmscaling;
            s->fcontext->outfbmfactor = uih->fcontext->outfbmfactor;
        }
        if (s->fcontext->intessellation != uih->fcontext->intessellation)
            save_intc(uih, "intessellation", uih->fcontext->intessellation),
                s->fcontext->intessellation = uih->fcontext->intessellation;
        if (s->fcontext->intessellcolor != uih->fcontext->intessellcolor ||
            s->fcontext->intessellfill != uih->fcontext->intessellfill) {
            start_save(uih, "intessellstyle");
            save_nstring(uih, uih->fcontext->intessellcolor, tessell_colors);
            save_nstring(uih, uih->fcontext->intessellfill, tessell_fills);
            stop_save(uih);
            s->fcontext->intessellcolor = uih->fcontext->intessellcolor;
            s->fcontext->intessellfill = uih->fcontext->intessellfill;
        }
        /* The decompositions' settings, one command a side, written only
         * when one of them has moved: a position that leaves them alone
         * opens in 2.5, which knows nothing of them. */
        if (s->fcontext->indecompsectors != uih->fcontext->indecompsectors ||
            s->fcontext->indecomprotation != uih->fcontext->indecomprotation ||
            s->fcontext->indecompof != uih->fcontext->indecompof ||
            s->fcontext->indecompradial != uih->fcontext->indecompradial ||
            s->fcontext->indecomprepeats != uih->fcontext->indecomprepeats ||
            s->fcontext->indecompway != uih->fcontext->indecompway ||
            s->fcontext->indecompspeed != uih->fcontext->indecompspeed) {
            start_save(uih, "indecomposition");
            save_int(uih, uih->fcontext->indecompsectors);
            save_float(uih, uih->fcontext->indecomprotation);
            save_nstring(uih, uih->fcontext->indecompof, decomp_of);
            save_int(uih, uih->fcontext->indecompradial);
            save_float(uih, uih->fcontext->indecomprepeats);
            save_nstring(uih, uih->fcontext->indecompway, decomp_ways);
            save_float(uih, uih->fcontext->indecompspeed);
            stop_save(uih);
            s->fcontext->indecomprepeats = uih->fcontext->indecomprepeats;
            s->fcontext->indecompway = uih->fcontext->indecompway;
            s->fcontext->indecompspeed = uih->fcontext->indecompspeed;
            s->fcontext->indecompsectors = uih->fcontext->indecompsectors;
            s->fcontext->indecomprotation = uih->fcontext->indecomprotation;
            s->fcontext->indecompof = uih->fcontext->indecompof;
            s->fcontext->indecompradial = uih->fcontext->indecompradial;
        }
        if (s->fcontext->outdecompsectors != uih->fcontext->outdecompsectors ||
            s->fcontext->outdecomprotation !=
                uih->fcontext->outdecomprotation ||
            s->fcontext->outdecompradial != uih->fcontext->outdecompradial ||
            s->fcontext->outdecomprepeats != uih->fcontext->outdecomprepeats ||
            s->fcontext->outdecompway != uih->fcontext->outdecompway ||
            s->fcontext->outdecompspeed != uih->fcontext->outdecompspeed) {
            start_save(uih, "outdecomposition");
            save_int(uih, uih->fcontext->outdecompsectors);
            save_float(uih, uih->fcontext->outdecomprotation);
            save_int(uih, uih->fcontext->outdecompradial);
            save_float(uih, uih->fcontext->outdecomprepeats);
            save_nstring(uih, uih->fcontext->outdecompway, decomp_ways);
            save_float(uih, uih->fcontext->outdecompspeed);
            stop_save(uih);
            s->fcontext->outdecomprepeats = uih->fcontext->outdecomprepeats;
            s->fcontext->outdecompway = uih->fcontext->outdecompway;
            s->fcontext->outdecompspeed = uih->fcontext->outdecompspeed;
            s->fcontext->outdecompsectors = uih->fcontext->outdecompsectors;
            s->fcontext->outdecomprotation = uih->fcontext->outdecomprotation;
            s->fcontext->outdecompradial = uih->fcontext->outdecompradial;
        }
        if ((s->fcontext->incoloringmode == INCOLORING_TRUECOLOR || s->mode >= UIH_SAVEALL) &&
            s->fcontext->intcolor != uih->fcontext->intcolor)
            save_intc(uih, "intcoloring", uih->fcontext->intcolor),
                s->fcontext->intcolor = uih->fcontext->intcolor;
        if ((s->fcontext->coloringmode == OutColormodeType::ColOut_True_color || s->mode >= UIH_SAVEALL) &&
            s->fcontext->outtcolor != uih->fcontext->outtcolor)
            save_intc(uih, "outtcoloring", uih->fcontext->outtcolor),
                s->fcontext->outtcolor = uih->fcontext->outtcolor;
        if (s->fcontext->incolorfun != uih->fcontext->incolorfun)
            save_intc(uih, "incolorfun", uih->fcontext->incolorfun),
                s->fcontext->incolorfun = uih->fcontext->incolorfun;
        if (s->fcontext->incolorspeed != uih->fcontext->incolorspeed)
            save_floatc(uih, "incolorspeed", uih->fcontext->incolorspeed),
                s->fcontext->incolorspeed = uih->fcontext->incolorspeed;
        if (s->fcontext->incolorshift != uih->fcontext->incolorshift)
            save_intc(uih, "incolorshift", uih->fcontext->incolorshift),
                s->fcontext->incolorshift = uih->fcontext->incolorshift;
        /* Written only when it is not the one every earlier version drew
         * with, so a position that does not change it stays loadable by
         * them. */
        if (s->fcontext->incolorfrequency != uih->fcontext->incolorfrequency)
            save_floatc(uih, "incolorfrequency", uih->fcontext->incolorfrequency),
                s->fcontext->incolorfrequency = uih->fcontext->incolorfrequency;
        if (s->fcontext->incolorperiodicfun != uih->fcontext->incolorperiodicfun)
            save_intc(uih, "incolorperiodicfun",
                      uih->fcontext->incolorperiodicfun),
                s->fcontext->incolorperiodicfun =
                    uih->fcontext->incolorperiodicfun;
        if (s->fcontext->outcolorfun != uih->fcontext->outcolorfun)
            save_intc(uih, "outcolorfun", uih->fcontext->outcolorfun),
                s->fcontext->outcolorfun = uih->fcontext->outcolorfun;
        if (s->fcontext->outcolorspeed != uih->fcontext->outcolorspeed)
            save_floatc(uih, "outcolorspeed", uih->fcontext->outcolorspeed),
                s->fcontext->outcolorspeed = uih->fcontext->outcolorspeed;
        if (s->fcontext->outcolorshift != uih->fcontext->outcolorshift)
            save_intc(uih, "outcolorshift", uih->fcontext->outcolorshift),
                s->fcontext->outcolorshift = uih->fcontext->outcolorshift;
        save_colorfunchange(uih, 1);
        save_colorfunchange(uih, 0);
        if (s->fcontext->pndefault != uih->fcontext->pndefault)
            save_onoffc(uih, "pndefault", uih->fcontext->pndefault),
                s->fcontext->pndefault = uih->fcontext->pndefault;
        if (s->fcontext->newtonmodesffe != uih->fcontext->newtonmodesffe)
            save_onoffc(uih, "newtonmodesffe", uih->fcontext->newtonmodesffe),
                s->fcontext->newtonmodesffe = uih->fcontext->newtonmodesffe;
        if (s->fcontext->newtonconvergence != uih->fcontext->newtonconvergence)
            save_floatc(uih, "newtonconvergence", uih->fcontext->newtonconvergence),
                s->fcontext->newtonconvergence = uih->fcontext->newtonconvergence;
        if (s->fcontext->mandelbrot != uih->fcontext->mandelbrot)
            save_onoffc(uih, "julia", !uih->fcontext->mandelbrot),
                s->fcontext->mandelbrot = uih->fcontext->mandelbrot;
        if (s->mode > UIH_SAVEPOS && s->fcontext->range != uih->fcontext->range)
            save_intc(uih, "range", uih->fcontext->range),
                s->fcontext->range = uih->fcontext->range;
        if (s->fcontext->plane != uih->fcontext->plane)
            save_intc(uih, "plane", uih->fcontext->plane),
                s->fcontext->plane = uih->fcontext->plane;
        if (s->zoomactive != uih->zoomactive && s->mode > UIH_SAVEPOS) {
            switch (uih->zoomactive) {
                case -1:
                    save_noparam(uih, "unzoom");
                    break;
                case 1:
                    save_noparam(uih, "zoom");
                    break;
                default:
                    save_noparam(uih, "stop");
                    break;
            }
            s->zoomactive = uih->zoomactive;
        }
        if ((s->mode >= UIH_SAVEPOS || uih->displaytext) &&
            s->color != uih->color) {
            start_save(uih, "color");
            save_nstring(uih, uih->color, uih_colornames);
            stop_save(uih);
            s->color = uih->color;
        }
        if (s->clearscreen) {
            save_noparam(uih, "clearscreen");
            s->clearscreen = 0;
            s->nonfractalscreen = 1;
        }
        if (uih->displaytext) {
            for (i = 0; i < 3; i++) {
                if (uih->displaytext & (1 << i)) {
                    if (s->ytextpos != i || s->xtextpos != uih->textpos[i]) {
                        start_save(uih, "textposition");
                        save_nstring(uih, uih->xtextpos, xtextposnames);
                        save_nstring(uih, uih->ytextpos, ytextposnames);
                        stop_save(uih);
                        s->xtextpos = uih->xtextpos;
                        s->ytextpos = uih->ytextpos;
                    }
                    save_stringc(uih, "text", uih->text[i]);
                    s->nonfractalscreen = 1;
                }
            }
            save_noparam(uih, "textsleep");
            uih->displaytext = 0;
        }
        if (s->autorotate && changed &&
            tl_lookup_timer(uih->savec->synctimer) > 500000)
            save_float2c(uih, "angle", uih->fcontext->angle, 5), resetsync = 1;
        if (s->mode == UIH_SAVEPOS)
            savepos(uih);
        else {
            if (uih->viewchanged)
                savepos(uih), uih->viewchanged = 0;
            else if (uih->moved)
                savepos3(uih), uih->moved = 0;
            else if (((changed && uih->step) || last) &&
                     tl_lookup_timer(uih->savec->synctimer) > 500000)
                resetsync = 1, savepos2(uih);
        }
        if (uih->savec->firsttime)
            uih->savec->firsttime = 0;
        if (s->writefailed)
            uih_save_disable(uih);
        if (resetsync)
            tl_reset_timer(uih->savec->synctimer);
    } /*if uih->save */
}

int uih_save_enable(struct uih_context *uih, xio_file f, int mode)
{
    struct uih_savedcontext *s;
    int i;
    last = 0;
    if (uih->save) {
        uih_error(uih, TR("Error", "Recording is already enabled"));
        return 0;
    }
    s = (struct uih_savedcontext *)calloc(1, sizeof(*s));
    if (f == NULL || s == NULL) {
        uih_error(uih, TR("Error", "File could not be opened or out of memory"));
        return 0;
    }
    uih->savec = s;
    s->fcontext = make_fractalc(1, uih->image->pixelwidth * uih->image->width,
                                uih->image->pixelheight * uih->image->height);
    if (s->fcontext == NULL) {
        uih_error(uih, TR("Error", "File could not be opened or out of memory"));
        return 0;
    }
    s->mode = mode;
    /*Invalidate context to force save everything first */
    s->speedup = STEP;
    s->maxstep = MAXSTEP;
    s->xcenter = INT_MAX;
    s->fastmode = 2;
    s->juliamode = 0;
    s->cycling = 0;
    for (i = 0; i < uih_nfilters; i++)
        s->filter[i] = 0;
    s->pressed = 0;
    s->firsttime = 1;
    uih->palettechanged = 1;
    s->cyclingspeed = 30;
    s->fcontext->pre = s->fcontext->pim = 0;
    s->fcontext->bre = s->fcontext->bim = 0;
    s->fcontext->currentformula = NULL;
    s->fcontext->periodicity = 1;
    s->fcontext->maxiter = 170;
    s->fcontext->bailout = 4;
    s->fcontext->coloringmode = OutColormodeType::ColOut_iter;
    s->fcontext->incoloringmode = 0;
    s->fcontext->outtcolor = 0;
    s->fcontext->intcolor = 0;
    s->fcontext->mandelbrot = 1;
    s->fcontext->incolorfun = 0;
    s->fcontext->incolorspeed = 1.0f;
    s->fcontext->incolorshift = 0;
    s->fcontext->incolorfrequency = 1;
    s->fcontext->incolorperiodicfun = 0;
    s->fcontext->indecompsectors = s->fcontext->outdecompsectors = 0;
    s->fcontext->indecomprotation = s->fcontext->outdecomprotation = 0;
    s->fcontext->indecompof = 0;
    s->fcontext->indecompradial = s->fcontext->outdecompradial = 0;
    s->fcontext->indecomprepeats = s->fcontext->outdecomprepeats = 1;
    s->fcontext->indecompway = s->fcontext->outdecompway = 0;
    s->fcontext->indecompspeed = s->fcontext->outdecompspeed = 1;
    s->fcontext->outcolorfun = 0;
    s->fcontext->outcolorspeed = 1.0f;
    s->fcontext->outcolorshift = 0;
    s->fcontext->incolorfuninitial = s->fcontext->outcolorfuninitial = 1;
    s->fcontext->incolorfuninitialspeed = s->fcontext->outcolorfuninitialspeed = 1;
    s->fcontext->incolorfunspeed = s->fcontext->outcolorfunspeed = 1;
    s->fcontext->incolorfun2speed = s->fcontext->outcolorfun2speed = 1;
    s->fcontext->incolorfuniter = s->fcontext->outcolorfuniter = COLORFUN_ITER;
    s->fcontext->incolorfunexp = s->fcontext->outcolorfunexp = COLORFUN_EXP;
    s->fcontext->incolorfunmode = s->fcontext->outcolorfunmode = 0;
    s->fcontext->incolorfunradius = s->fcontext->outcolorfunradius = COLORFUN_RADIUS;
    s->fcontext->incolorfun2iter = s->fcontext->outcolorfun2iter = COLORFUN2_ITER;
    s->fcontext->incolorfun2exp = s->fcontext->outcolorfun2exp = COLORFUN2_EXP;
    s->fcontext->incolorfun2mode = s->fcontext->outcolorfun2mode = 0;
    s->fcontext->incolorfun2radius = s->fcontext->outcolorfun2radius = COLORFUN2_RADIUS;
    s->fcontext->infbmintensity = s->fcontext->outfbmintensity = 4;
    s->fcontext->infbmfrequency = s->fcontext->outfbmfrequency = 8;
    s->fcontext->infbmroughness = s->fcontext->outfbmroughness = (number_t)1 / 2;
    s->fcontext->infbmoctaves = s->fcontext->outfbmoctaves = 4;
    s->fcontext->infbmseed = s->fcontext->outfbmseed = 1;
    s->fcontext->outfbmscaling = FBM_SCALE_FIXED;
    s->fcontext->outfbmfactor = 1;
    s->fcontext->intessellation = 3;
    s->fcontext->intessellcolor = TESSELL_ANGLE;
    s->fcontext->intessellfill = TESSELL_TURNED;
    s->fcontext->pndefault = 0;
    s->fcontext->newtonmodesffe = 0;
    s->fcontext->newtonconvergence = 1E-6;
    s->fcontext->plane = 0;
    s->fcontext->range = 3;
    s->fcontext->angle = 0;
    s->rotatepressed = 0;
    s->autorotate = 0;
    s->fastrotate = 0;
    s->rotationspeed = 10;
    s->clearscreen = 0;
    s->color = 0;
    s->xtextpos = 1;
    s->ytextpos = 1;
    s->file = f;
    s->timer = tl_create_timer();
    s->synctimer = tl_create_timer();
    uih->viewchanged = 1;
    uih->palettechanged = 1;
    uih->save = 1;
    uih_emulatetimers(uih);
    tl_reset_timer(s->timer);
    uih->moved = 0;
    if (mode == UIH_SAVEANIMATION)
        myputs(";Animation file automatically generated by " XaoS_NAME " " XaoS_VERSION
               "\n"
               ";  - a realtime interactive fractal zoomer\n"
               ";Use xaos -play <filename> to replay it\n");
    else if (mode == UIH_SAVEPOS)
        myputs(";Position file automatically generated by " XaoS_NAME " " XaoS_VERSION
               "\n"
               ";  - a realtime interactive fractal zoomer\n"
               ";Use xaos -loadpos <filename> to display it\n");
#ifdef USE_FLOAT128
    /* The quad build records what it computed at, before anything else. It is
     * a command and not a comment so that a build which does not know it stops
     * rather than opening the file and drawing a different picture: 113 bits
     * of mantissa is not something a narrower build can reproduce, and saying
     * so is worth losing the ability to open it there. */
    save_intc(uih, "precision", NUMBER_MANTISSA_BITS);
#endif
    uih_saveframe(uih);
    uih_updatemenus(uih, "save");
    xio_putc('\n', f);
    return 1;
}

void uih_save_disable(struct uih_context *uih)
{
    if (uih->save) {
        last = 1;
        if (uih->savec->mode >= UIH_SAVEANIMATION)
            uih_saveframe(uih);
        if (xio_close(uih->savec->file))
            outputerror(uih);
        uih->save = 0;
        free_fractalc(uih->savec->fcontext);
        tl_free_timer(uih->savec->timer);
        tl_free_timer(uih->savec->synctimer);
        free(uih->savec);
        uih_updatemenus(uih, "save");
    }
}

void uih_save_position(struct uih_context *uih, xio_file f, int mode)
{
    struct uih_savedcontext *c = uih->savec;
    int save = uih->save;
    int vc = uih->viewchanged;
    int pc = uih->palettechanged;
    uih->moved = 0;
    uih->save = 0;
    uih->savec = NULL;
    uih_save_enable(uih, f, mode);
    uih_save_disable(uih);
    uih->savec = c;
    uih->save = save;
    uih->viewchanged = vc;
    uih->palettechanged = pc;
}
