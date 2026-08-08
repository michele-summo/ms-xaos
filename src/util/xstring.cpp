#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "config.h"
#include "xio.h"
#include "misc-f.h"

#ifdef USE_FLOAT128
#include <quadmath.h>
#endif

struct fr {
    char *string;
    int pos;
    int allocedsize;
    int stringsize;
};

number_t xstrtonum(const char *s, char **sp)
{
#ifdef USE_FLOAT128
    return strtoflt128(s, sp);
#else
#ifdef USE_LONG_DOUBLE
    return strtold(s, sp);
#else
    return strtod(s, sp);
#endif
#endif
}

/* The counterpart of xstrtonum: writes a number_t into buf with the requested
 * number of significant digits, and returns buf so it can be used inline.
 *
 * %G rather than a fixed number of decimal places, because a view that has
 * been zoomed into is far smaller than any fixed format can show -- printed as
 * %1.12f, a size of 1e-13 reads as 0.000000000000, which is what the status
 * bar used to display once past the first dozen digits of the zoom. */
const char *xnumtostr(char *buf, int size, number_t number, int digits)
{
    char fs[16];
    if (digits < 1)
        digits = 1;
    if (digits > NUMBER_DIGITS)
        digits = NUMBER_DIGITS;
#ifdef USE_FLOAT128
    snprintf(fs, sizeof(fs), "%%.%iQG", digits);
    quadmath_snprintf(buf, size, fs, (__float128)number);
#else
#ifdef USE_LONG_DOUBLE
    snprintf(fs, sizeof(fs), "%%.%iLG", digits);
    snprintf(buf, size, fs, (long double)number);
#else
    snprintf(fs, sizeof(fs), "%%.%iG", digits);
    snprintf(buf, size, fs, (double)number);
#endif
#endif
    return buf;
}

char *mystrdup(const char *c)
{
    int l = strlen(c);
    char *d = (char *)malloc(l + 1);
    if (!d)
        return NULL;
    memcpy(d, c, l + 1);
    return d;
}

static int sputc(int c, xio_file s)
{
    struct fr *f = (struct fr *)s->data;
    if (f->pos >= f->allocedsize - 1) {
        char *c = (char *)realloc(f->string, f->allocedsize * 2);
        if (!c)
            return XIO_EOF;
        f->string = c;
        f->allocedsize *= 2;
    }
    f->string[f->pos++] = c;
    if (f->pos >= f->stringsize)
        f->string[f->pos] = 0, f->stringsize = f->pos;
    return 0;
}

static int sputs(const char *c, xio_file s)
{
    int l = strlen(c);
    struct fr *f = (struct fr *)s->data;
    while (f->pos + l >= f->allocedsize - 1) {
        char *c = (char *)realloc(f->string, f->allocedsize * 2);
        if (!c)
            return XIO_EOF;
        f->string = c;
        f->allocedsize *= 2;
    }
    memcpy(f->string + f->pos, c, l);
    f->pos += l;
    if (f->pos >= f->stringsize)
        f->string[f->pos] = 0, f->stringsize = f->pos;
    return 0;
}

static int sungetc(int /*c*/, xio_file s)
{
    struct fr *f = (struct fr *)s->data;
    f->pos--;
    /*f->string[f->pos]=c; */
    return 0;
}

static int sgetc(xio_file s)
{
    struct fr *f = (struct fr *)s->data;
    if (f->pos == f->stringsize)
        return XIO_EOF;
    return f->string[f->pos++];
}

static int sfeof(xio_file s)
{
    struct fr *f = (struct fr *)s->data;
    return (f->pos == f->stringsize);
}

static int srclose(xio_file s)
{
    struct fr *f = (struct fr *)s->data;
    free(f->string);
    free(f);
    free(s);
    return 0;
}

static int swclose(xio_file s)
{
    struct fr *f = (struct fr *)s->data;
    f->string = (char *)realloc(f->string, f->stringsize + 1);
    /*free(s);
       free(f); */
    return 0;
}

char *xio_getstring(xio_file s)
{
    struct fr *f = (struct fr *)s->data;
    char *c = f->string;
    free(f);
    free(s);
    return c;
}

xio_file xio_strropen(const char *string)
{
    xio_file s = (xio_file)calloc(1, sizeof(*s));
    struct fr *f = (struct fr *)calloc(1, sizeof(*f));
    s->data = f;
    f->pos = 0;
    f->string = (char *)string;
    f->stringsize = strlen(string);
    s->fclose = srclose;
    s->xeof = sfeof;
    s->fgetc = sgetc;
    s->fungetc = sungetc;
    return s;
}

#define PAGE 4096
xio_file xio_strwopen(void)
{
    xio_file s = (xio_file)calloc(1, sizeof(*s));
    struct fr *f = (struct fr *)calloc(1, sizeof(*f));
    s->data = f;
    f->pos = 0;
    f->string = (char *)malloc(PAGE);
    f->allocedsize = PAGE;
    f->stringsize = 0;
    s->fputc = sputc;
    s->fputs = sputs;
    s->fclose = swclose;
    s->flush = NULL;
    return s;
}
